#include "../dictionary.hpp"
#include "../gz/zip_stream.hpp"
#include "../../external/pthash/external/essentials/include/essentials.hpp"

#include "util.hpp"
#include "builder.hpp"
#include "build_index.hpp"
#include "build_skew_index.hpp"

#include <numeric>  // for std::accumulate

namespace sshash {

void parse_file(std::istream& is, builder& b, build_configuration const& build_config) {
    for (std::string sequence; !is.eof();) {
        std::getline(is, sequence);  // header sequence
        if (build_config.weighted) b.parse_header(sequence.data(), sequence.length());
        std::getline(is, sequence);  // DNA sequence
        if (build_config.weighted and sequence.length() != b.expected_sequence_length()) {
            std::cout << "ERROR: expected a sequence of length " << b.expected_sequence_length()
                      << " but got one of length " << sequence.length() << std::endl;
            throw std::runtime_error("file is malformed");
        }
        b.add_sequence(sequence.data(), sequence.length());
    }
}

void parse_file(std::string const& filename, builder& b, build_configuration const& build_config) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    std::cout << "reading file '" << filename << "'..." << std::endl;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        parse_file(zis, b, build_config);
    } else {
        parse_file(is, b, build_config);
    }
    is.close();
}

void dictionary::build_from(builder& b, build_configuration const& build_config) {
    build_config.validate();

    m_k = build_config.k;
    m_m = build_config.m;
    m_seed = build_config.seed;
    m_canonical_parsing = build_config.canonical_parsing;
    m_skew_index.min_log2 = build_config.l;
    m_size = b.data.num_kmers;

    b.finalize();
    build_from(b.data, build_config);
}

void dictionary::build_from(std::string const& filename, build_configuration const& build_config) {
    build_config.validate();

    m_k = build_config.k;
    m_m = build_config.m;
    m_seed = build_config.seed;
    m_canonical_parsing = build_config.canonical_parsing;
    m_skew_index.min_log2 = build_config.l;

    /* step 1: parse the input file and build compact string pool ***/
    essentials::timer_type timer;
    timer.start();
    builder b(build_config);
    parse_file(filename, b, build_config);
    b.finalize();
    m_size = b.data.num_kmers;
    timer.stop();
    b.data.parse_time = timer.elapsed();
    /******/

    build_from(b.data, build_config);
}

void dictionary::build_from(parse_data& data, build_configuration const& build_config) {
    std::vector<double> timings;
    timings.reserve(5);
    essentials::timer_type timer;

    if (data.parse_time) {
        timings.push_back(data.parse_time);
        print_time(timings.back(), data.num_kmers, "step 1: 'parse_file'");
    }

    if (build_config.weighted) {
        /* step 1.1: compress weights ***/
        timer.start();
        data.weights_builder.build(m_weights);
        timer.stop();
        timings.push_back(timer.elapsed());
        print_time(timings.back(), data.num_kmers, "step 1.1.: 'build_weights'");
        timer.reset();
        /******/
        if (build_config.verbose) {
            double entropy_weights = data.weights_builder.print_info(data.num_kmers);
            double avg_bits_per_weight = static_cast<double>(m_weights.num_bits()) / data.num_kmers;
            std::cout << "weights: " << avg_bits_per_weight << " [bits/kmer]" << std::endl;
            std::cout << "  (" << entropy_weights / avg_bits_per_weight
                      << "x smaller than the empirical entropy)" << std::endl;
        }
    }

    /* step 2: merge minimizers and build MPHF ***/
    timer.start();
    data.minimizers.merge();
    {
        mm::file_source<minimizer_tuple> input(data.minimizers.get_minimizers_filename(),
                                               mm::advice::sequential);
        minimizers_tuples_iterator iterator(input.data(), input.data() + input.size());
        m_minimizers.build(iterator, data.minimizers.num_minimizers(), build_config);
        input.close();
    }
    timer.stop();
    timings.push_back(timer.elapsed());
    print_time(timings.back(), data.num_kmers, "step 2: 'build_minimizers'");
    timer.reset();
    /******/

    /* step 3: build index ***/
    timer.start();
    auto buckets_stats = build_index(data, m_minimizers, m_buckets, build_config);
    timer.stop();
    timings.push_back(timer.elapsed());
    print_time(timings.back(), data.num_kmers, "step 3: 'build_index'");
    timer.reset();
    /******/

    /* step 4: build skew index ***/
    timer.start();
    build_skew_index(m_skew_index, data, m_buckets, build_config, buckets_stats);
    timer.stop();
    timings.push_back(timer.elapsed());
    print_time(timings.back(), data.num_kmers, "step 4: 'build_skew_index'");
    timer.reset();
    /******/

    if (data.parse_time) {
        double total_time = std::accumulate(timings.begin(), timings.end(), 0.0);
        print_time(total_time, data.num_kmers, "total_time");
    }

    print_space_breakdown();

    if (build_config.verbose) buckets_stats.print_less();

    data.minimizers.remove_tmp_file();
}

}  // namespace sshash
