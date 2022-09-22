#include "../gz/zip_stream.cpp"

namespace sshash {

struct parse_data {
    parse_data(std::string const& tmp_dirname) : num_kmers(0), minimizers(tmp_dirname) {}
    uint64_t num_kmers;
    minimizers_tuples minimizers;
    compact_string_pool strings;
    weights::builder weights_builder;
};

void parse_file(std::istream& is, parse_data& data, build_configuration const& build_config) {
    uint64_t k = build_config.k;
    uint64_t m = build_config.m;
    uint64_t seed = build_config.seed;
    uint64_t max_num_kmers_in_super_kmer = k - m + 1;
    uint64_t block_size = 2 * k - m;  // max_num_kmers_in_super_kmer + k - 1

    if (max_num_kmers_in_super_kmer >= (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8))) {
        throw std::runtime_error(
            "max_num_kmers_in_super_kmer " + std::to_string(max_num_kmers_in_super_kmer) +
            " does not fit into " + std::to_string(sizeof(num_kmers_in_super_kmer_uint_type) * 8) +
            " bits");
    }

    /* fit into the wanted number of bits */
    assert(max_num_kmers_in_super_kmer < (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8)));

    compact_string_pool::builder builder(k);

    std::string sequence;
    uint64_t prev_minimizer = constants::invalid;

    uint64_t begin = 0;  // begin of parsed super_kmer in sequence
    uint64_t end = 0;    // end of parsed super_kmer in sequence
    uint64_t num_sequences = 0;
    uint64_t num_bases = 0;
    bool glue = false;

    std::vector<uint64_t> positions_distr1(k - m + 1, 0);
    std::vector<uint64_t> positions_distr2(k - m + 1, 0);
    std::vector<uint64_t> cardinalities_distr(k - m + 1 + 1, 0);

    /* left-right-max */
    uint64_t kmers_in_lr_skms = 0;
    uint64_t num_lr_skms = 0;

    /* right-max */
    uint64_t kmers_in_r_skms = 0;
    uint64_t num_r_skms = 0;

    /* left-max */
    uint64_t kmers_in_l_skms = 0;
    uint64_t num_l_skms = 0;

    /* non-max (all others) */
    uint64_t kmers_in_all_other_skms = 0;
    uint64_t num_all_other_skms = 0;

    /* pathological cases where minimizer repeats in skm */
    uint64_t num_skms_with_repeated_minimizer = 0;

    uint64_t num_skms_with_non_repeated_minimizer = 0;
    uint64_t num_skms = 0;  // this is <= than data.strings.num_super_kmers() because of the
                            // splitting of skms into blocks

    auto append_super_kmer = [&]() {
        if (sequence.empty() or prev_minimizer == constants::invalid or begin == end) return;

        assert(end > begin);
        char const* super_kmer = sequence.data() + begin;
        uint64_t size = (end - begin) + k - 1;
        assert(util::is_valid(super_kmer, size));

        // std::cout << std::string(super_kmer, size) << '\n';
        // for (uint64_t i = 0; i != size - k + 1; ++i) {
        //     char const* kmer = super_kmer + i;
        //     uint64_t uint64_kmer = util::string_to_uint64_no_reverse(kmer, k);
        //     auto [minimizer, pos] = util::compute_minimizer_pos(uint64_kmer, k, m, seed);
        //     std::cout << util::uint64_to_string_no_reverse(minimizer, m) << ' ';
        //     std::cout << pos << '\n';
        // }
        // std::cout << '\n';

        /* if num_kmers_in_super_kmer > k - m + 1, then split the super_kmer into blocks */
        uint64_t num_kmers_in_super_kmer = end - begin;
        num_skms += 1;

        uint64_t uint64_first_kmer = util::string_to_uint64_no_reverse(super_kmer, k);
        auto [minimizer_first_kmer, pos_of_min_in_first_kmer_of_skm] =
            util::compute_minimizer_pos(uint64_first_kmer, k, m, seed);
        (void)minimizer_first_kmer;
        assert(pos_of_min_in_first_kmer_of_skm <= k - m);

        uint64_t uint64_last_kmer =
            util::string_to_uint64_no_reverse(super_kmer + num_kmers_in_super_kmer - 1, k);
        auto [minimizer_last_kmer, pos_of_min_in_last_kmer_of_skm] =
            util::compute_minimizer_pos(uint64_last_kmer, k, m, seed);
        (void)minimizer_last_kmer;

        // positions_distr1[pos_of_min_in_first_kmer_of_skm] += 1;
        // positions_distr2[pos_of_min_in_last_kmer_of_skm] += 1;

        if (num_kmers_in_super_kmer <= k - m + 1) {
            num_skms_with_non_repeated_minimizer += 1;

            // cardinalities_distr[num_kmers_in_super_kmer] += 1;

            if (pos_of_min_in_first_kmer_of_skm == k - m and pos_of_min_in_last_kmer_of_skm == 0) {
                /* left-right-max */
                assert(num_kmers_in_super_kmer == k - m + 1);
                kmers_in_lr_skms += num_kmers_in_super_kmer;
                num_lr_skms += 1;
            } else if (pos_of_min_in_first_kmer_of_skm == k - m) {
                /* right-max */
                assert(pos_of_min_in_last_kmer_of_skm != 0);
                kmers_in_r_skms += num_kmers_in_super_kmer;
                num_r_skms += 1;
            } else if (pos_of_min_in_last_kmer_of_skm == 0) {
                /* left-max */
                assert(pos_of_min_in_first_kmer_of_skm != k - m);
                kmers_in_l_skms += num_kmers_in_super_kmer;
                num_l_skms += 1;
            } else {
                /* non-max */
                assert(pos_of_min_in_first_kmer_of_skm != k - m);
                assert(pos_of_min_in_last_kmer_of_skm != 0);
                kmers_in_all_other_skms += num_kmers_in_super_kmer;
                num_all_other_skms += 1;
            }
        } else {
            // pathological cases where num_kmers_in_super_kmer > k - m + 1 should be rare if m is
            // sufficiently large...

            num_skms_with_repeated_minimizer += 1;

            // if (pos_of_min_in_first_kmer_of_skm == k - m and pos_of_min_in_last_kmer_of_skm == 0)
            // {
            //     /* left-right-max */
            //     kmers_in_lr_skms += num_kmers_in_super_kmer;
            //     num_lr_skms += 1;
            // } else if (pos_of_min_in_first_kmer_of_skm == k - m) {
            //     /* right-max */
            //     assert(pos_of_min_in_last_kmer_of_skm != 0);
            //     kmers_in_r_skms += num_kmers_in_super_kmer;
            //     num_r_skms += 1;
            // } else if (pos_of_min_in_last_kmer_of_skm == 0) {
            //     /* left-max */
            //     assert(pos_of_min_in_first_kmer_of_skm != k - m);
            //     kmers_in_l_skms += num_kmers_in_super_kmer;
            //     num_l_skms += 1;
            // } else {
            //     /* non-max */
            //     assert(pos_of_min_in_first_kmer_of_skm != k - m);
            //     assert(pos_of_min_in_last_kmer_of_skm != 0);
            //     kmers_in_all_other_skms += num_kmers_in_super_kmer;
            //     num_all_other_skms += 1;
            // }
        }

        uint64_t num_blocks = num_kmers_in_super_kmer / max_num_kmers_in_super_kmer +
                              (num_kmers_in_super_kmer % max_num_kmers_in_super_kmer != 0);
        assert(num_blocks > 0);
        for (uint64_t i = 0; i != num_blocks; ++i) {
            uint64_t n = block_size;
            if (i == num_blocks - 1) n = size;
            uint64_t num_kmers_in_block = n - k + 1;
            assert(num_kmers_in_block <= max_num_kmers_in_super_kmer);
            data.minimizers.emplace_back(prev_minimizer, builder.offset, num_kmers_in_block);
            builder.append(super_kmer + i * max_num_kmers_in_super_kmer, n, glue);
            if (glue) {
                assert(data.minimizers.back().offset > k - 1);
                data.minimizers.back().offset -= k - 1;
            }
            size -= max_num_kmers_in_super_kmer;
            glue = true;
        }
    };

    uint64_t seq_len = 0;
    uint64_t sum_of_weights = 0;
    data.weights_builder.init();

    /* intervals of weights */
    uint64_t weight_value = constants::invalid;
    uint64_t weight_length = 0;

    auto parse_header = [&]() {
        if (sequence.empty()) return;

        /*
            Heder format:
            >[id] LN:i:[seq_len] ab:Z:[weight_seq]
            where [weight_seq] is a space-separated sequence of integer counters (the weights),
            whose length is equal to [seq_len]-k+1
        */

        // example header: '>12 LN:i:41 ab:Z:2 2 2 2 2 2 2 2 2 2 2'

        expect(sequence[0], '>');
        uint64_t i = 0;
        i = sequence.find_first_of(' ', i);
        if (i == std::string::npos) throw parse_runtime_error();

        i += 1;
        expect(sequence[i + 0], 'L');
        expect(sequence[i + 1], 'N');
        expect(sequence[i + 2], ':');
        expect(sequence[i + 3], 'i');
        expect(sequence[i + 4], ':');
        i += 5;
        uint64_t j = sequence.find_first_of(' ', i);
        if (j == std::string::npos) throw parse_runtime_error();

        seq_len = std::strtoull(sequence.data() + i, nullptr, 10);
        i = j + 1;
        expect(sequence[i + 0], 'a');
        expect(sequence[i + 1], 'b');
        expect(sequence[i + 2], ':');
        expect(sequence[i + 3], 'Z');
        expect(sequence[i + 4], ':');
        i += 5;

        for (uint64_t j = 0; j != seq_len - k + 1; ++j) {
            uint64_t weight = std::strtoull(sequence.data() + i, nullptr, 10);
            i = sequence.find_first_of(' ', i) + 1;

            data.weights_builder.eat(weight);
            sum_of_weights += weight;

            if (weight == weight_value) {
                weight_length += 1;
            } else {
                if (weight_value != constants::invalid) {
                    data.weights_builder.push_weight_interval(weight_value, weight_length);
                }
                weight_value = weight;
                weight_length = 1;
            }
        }
    };

    while (!is.eof()) {
        std::getline(is, sequence);  // header sequence
        if (build_config.weighted) parse_header();

        std::getline(is, sequence);  // DNA sequence
        if (sequence.size() < k) continue;

        if (++num_sequences % 100000 == 0) {
            std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
                      << data.num_kmers << " kmers" << std::endl;
        }

        begin = 0;
        end = 0;
        glue = false;  // start a new piece
        prev_minimizer = constants::invalid;
        num_bases += sequence.size();

        if (build_config.weighted and seq_len != sequence.size()) {
            std::cout << "ERROR: expected a sequence of length " << seq_len
                      << " but got one of length " << sequence.size() << std::endl;
            throw std::runtime_error("file is malformed");
        }

        while (end != sequence.size() - k + 1) {
            char const* kmer = sequence.data() + end;
            assert(util::is_valid(kmer, k));
            uint64_t uint64_kmer = util::string_to_uint64_no_reverse(kmer, k);
            uint64_t minimizer = util::compute_minimizer(uint64_kmer, k, m, seed);

            if (build_config.canonical_parsing) {
                uint64_t uint64_kmer_rc = util::compute_reverse_complement(uint64_kmer, k);
                uint64_t minimizer_rc = util::compute_minimizer(uint64_kmer_rc, k, m, seed);
                minimizer = std::min<uint64_t>(minimizer, minimizer_rc);
            }

            if (prev_minimizer == constants::invalid) prev_minimizer = minimizer;
            if (minimizer != prev_minimizer) {
                append_super_kmer();
                begin = end;
                prev_minimizer = minimizer;
                glue = true;
            }

            ++data.num_kmers;
            ++end;
        }

        append_super_kmer();
    }

    data.minimizers.finalize();
    builder.finalize();
    builder.build(data.strings);

    std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
              << data.num_kmers << " kmers" << std::endl;
    std::cout << "num_kmers " << data.num_kmers << std::endl;
    std::cout << "num_super_kmers " << num_skms << std::endl;
    std::cout << "num_super_kmers (after splitting into blocks) " << data.strings.num_super_kmers()
              << std::endl;
    std::cout << "num_pieces " << data.strings.pieces.size() << " (+"
              << (2.0 * data.strings.pieces.size() * (k - 1)) / data.num_kmers << " [bits/kmer])"
              << std::endl;
    assert(data.strings.pieces.size() == num_sequences + 1);

    if (build_config.weighted) {
        std::cout << "sum_of_weights " << sum_of_weights << std::endl;
        data.weights_builder.push_weight_interval(weight_value, weight_length);
        data.weights_builder.finalize(data.num_kmers);
    }

    assert(data.strings.num_super_kmers() >= num_skms);
    assert(num_skms == num_skms_with_repeated_minimizer + num_skms_with_non_repeated_minimizer);

    std::cout << "num_skms_with_repeated_minimizer " << num_skms_with_repeated_minimizer << " ("
              << (num_skms_with_repeated_minimizer * 100.0) / num_skms << "%)" << std::endl;
    std::cout << "num_skms_with_non_repeated_minimizer " << num_skms_with_non_repeated_minimizer
              << " (" << (num_skms_with_non_repeated_minimizer * 100.0) / num_skms << "%)"
              << std::endl;

    // std::cout << "positions_distr1:\n";
    // for (uint64_t i = 0; i != k - m + 1; ++i) {
    //     std::cout << "% of first km in skms whose min appears at pos " << i << ": "
    //               << (positions_distr1[i] * 100.0) / num_skms_with_non_repeated_minimizer <<
    //               '\n';
    // }
    // std::cout << "positions_distr2:\n";
    // for (uint64_t i = 0; i != k - m + 1; ++i) {
    //     std::cout << "% of last km in skms whose min appears at pos " << i << ": "
    //               << (positions_distr2[i] * 100.0) / num_skms_with_non_repeated_minimizer <<
    //               '\n';
    // }
    // std::cout << "cardinalities_distr:\n";
    // for (uint64_t i = 0; i != k - m + 1 + 1; ++i) {
    //     if (i == 0) continue;
    //     std::cout << "% of skms that contain " << i << " kmers: "
    //               << (cardinalities_distr[i] * 100.0) / num_skms_with_non_repeated_minimizer << "
    //               ("
    //               << (cardinalities_distr[i] * i * 100.0) / data.num_kmers << "% of total kmers)"
    //               << '\n';
    // }

    std::cout << "kmers_in_lr_skms = " << kmers_in_lr_skms << "/" << data.num_kmers << "("
              << (kmers_in_lr_skms * 100.0) / data.num_kmers << "%)" << std::endl;
    std::cout << "num_lr_skms = " << num_lr_skms << "/" << num_skms_with_non_repeated_minimizer
              << "(" << (num_lr_skms * 100.0) / num_skms_with_non_repeated_minimizer << "%)"
              << std::endl;

    std::cout << "kmers_in_r_skms = " << kmers_in_r_skms << "/" << data.num_kmers << "("
              << (kmers_in_r_skms * 100.0) / data.num_kmers << "%)" << std::endl;
    std::cout << "num_r_skms = " << num_r_skms << "/" << num_skms_with_non_repeated_minimizer << "("
              << (num_r_skms * 100.0) / num_skms_with_non_repeated_minimizer << "%)" << std::endl;

    std::cout << "kmers_in_l_skms = " << kmers_in_l_skms << "/" << data.num_kmers << "("
              << (kmers_in_l_skms * 100.0) / data.num_kmers << "%)" << std::endl;
    std::cout << "num_l_skms = " << num_l_skms << "/" << num_skms_with_non_repeated_minimizer << "("
              << (num_l_skms * 100.0) / num_skms_with_non_repeated_minimizer << "%)" << std::endl;

    std::cout << "kmers_in_all_other_skms = " << kmers_in_all_other_skms << "/" << data.num_kmers
              << "(" << (kmers_in_all_other_skms * 100.0) / data.num_kmers << "%)" << std::endl;
    std::cout << "num_all_other_skms = " << num_all_other_skms << "/"
              << num_skms_with_non_repeated_minimizer << "("
              << (num_all_other_skms * 100.0) / num_skms_with_non_repeated_minimizer << "%)"
              << std::endl;
}

parse_data parse_file(std::string const& filename, build_configuration const& build_config) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    std::cout << "reading file '" << filename << "'..." << std::endl;
    parse_data data(build_config.tmp_dirname);
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        parse_file(zis, data, build_config);
    } else {
        parse_file(is, data, build_config);
    }
    is.close();
    return data;
}

}  // namespace sshash