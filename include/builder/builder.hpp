#pragma once

#include "util.hpp"

namespace sshash {

struct parse_data {
    parse_data(std::string const& tmp_dirname)
        : num_kmers(0), minimizers(tmp_dirname), parse_time(0.0) {}
    uint64_t num_kmers;
    minimizers_tuples minimizers;
    compact_string_pool strings;
    weights::builder weights_builder;
    double parse_time;
};

struct builder {
    builder(build_configuration const& build_config)
        : data(build_config.tmp_dirname)

        , m_k(build_config.k)
        , m_m(build_config.m)
        , m_seed(build_config.seed)
        , m_max_num_kmers_in_super_kmer(m_k - m_m + 1)
        , m_block_size(2 * m_k - m_m)  // m_max_num_kmers_in_super_kmer + m_k - 1

        , m_prev_minimizer(constants::invalid_uint64)
        , m_begin(0)  // begin of parsed super_kmer in sequence
        , m_end(0)    // end of parsed super_kmer in sequence
        , m_num_sequences(0)
        , m_num_bases(0)
        , m_glue(false)

        , m_expected_sequence_length(0)
        , m_sum_of_weights(0)
        , m_weight_value(constants::invalid_uint64)
        , m_weight_length(0)

        , m_compact_string_pool_builder(m_k)
        , m_build_config(build_config)

    {
        data.weights_builder.init();

        if (m_max_num_kmers_in_super_kmer >=
            (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8))) {
            throw std::runtime_error(
                "max_num_kmers_in_super_kmer " + std::to_string(m_max_num_kmers_in_super_kmer) +
                " does not fit into " +
                std::to_string(sizeof(num_kmers_in_super_kmer_uint_type) * 8) + " bits");
        }

        /* fit into the wanted number of bits */
        assert(m_max_num_kmers_in_super_kmer <
               (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8)));
    }

    void add_sequence(char const* sequence, uint64_t length) {
        if (length < m_k) return;

        if (++m_num_sequences % 100000 == 0) {
            std::cout << "processed " << m_num_sequences << " sequences, " << m_num_bases
                      << " bases, " << data.num_kmers << " kmers" << std::endl;
        }

        m_begin = 0;
        m_end = 0;
        m_glue = false;  // start a new piece
        m_prev_minimizer = constants::invalid_uint64;
        m_num_bases += length;

        while (m_end != length - m_k + 1) {
            char const* kmer = sequence + m_end;
            assert(util::is_valid(kmer, m_k));
            kmer_t uint_kmer = util::string_to_uint_kmer_no_reverse(kmer, m_k);
            uint64_t minimizer = util::compute_minimizer(uint_kmer, m_k, m_m, m_seed);

            if (m_build_config.canonical_parsing) {
                kmer_t uint_kmer_rc = util::compute_reverse_complement(uint_kmer, m_k);
                uint64_t minimizer_rc = util::compute_minimizer(uint_kmer_rc, m_k, m_m, m_seed);
                minimizer = std::min<uint64_t>(minimizer, minimizer_rc);
            }

            if (m_prev_minimizer == constants::invalid_uint64) m_prev_minimizer = minimizer;
            if (minimizer != m_prev_minimizer) {
                append_super_kmer(sequence, length);
                m_begin = m_end;
                m_prev_minimizer = minimizer;
                m_glue = true;
            }

            ++data.num_kmers;
            ++m_end;
        }

        append_super_kmer(sequence, length);
    }

    void parse_header(char const* sequence, uint64_t length) {
        if (length == 0) {
            m_expected_sequence_length = 0;
            return;
        }

        /*
            Heder format:
            >[id] LN:i:[seq_len] ab:Z:[weight_seq]
            where [weight_seq] is a space-separated sequence of integer counters (the weights),
            whose length is equal to [seq_len]-k+1
        */

        // example header: '>12 LN:i:41 ab:Z:2 2 2 2 2 2 2 2 2 2 2'

        std::string str(sequence, length);

        expect(str[0], '>');
        uint64_t i = 0;
        i = str.find_first_of(' ', i);
        if (i == std::string::npos) throw parse_runtime_error();

        i += 1;
        expect(str[i + 0], 'L');
        expect(str[i + 1], 'N');
        expect(str[i + 2], ':');
        expect(str[i + 3], 'i');
        expect(str[i + 4], ':');
        i += 5;
        uint64_t j = str.find_first_of(' ', i);
        if (j == std::string::npos) throw parse_runtime_error();

        uint64_t seq_len = std::strtoull(str.data() + i, nullptr, 10);
        i = j + 1;
        expect(str[i + 0], 'a');
        expect(str[i + 1], 'b');
        expect(str[i + 2], ':');
        expect(str[i + 3], 'Z');
        expect(str[i + 4], ':');
        i += 5;

        for (uint64_t j = 0; j != seq_len - m_k + 1; ++j) {
            uint64_t weight = std::strtoull(str.data() + i, nullptr, 10);
            i = str.find_first_of(' ', i) + 1;

            data.weights_builder.eat(weight);
            m_sum_of_weights += weight;

            if (weight == m_weight_value) {
                m_weight_length += 1;
            } else {
                if (m_weight_value != constants::invalid_uint64) {
                    data.weights_builder.push_weight_interval(m_weight_value, m_weight_length);
                }
                m_weight_value = weight;
                m_weight_length = 1;
            }
        }

        m_expected_sequence_length = seq_len;
    }

    uint64_t expected_sequence_length() const { return m_expected_sequence_length; }

    void finalize() {
        data.minimizers.finalize();
        m_compact_string_pool_builder.finalize();
        m_compact_string_pool_builder.build(data.strings);

        std::cout << "processed " << m_num_sequences << " sequences, " << m_num_bases << " bases, "
                  << data.num_kmers << " kmers" << std::endl;
        std::cout << "num_kmers " << data.num_kmers << std::endl;
        std::cout << "num_super_kmers " << data.strings.num_super_kmers() << std::endl;
        std::cout << "num_pieces " << data.strings.pieces.size() << " (+"
                  << (2.0 * data.strings.pieces.size() * (m_k - 1)) / data.num_kmers
                  << " [bits/kmer])" << std::endl;
        assert(data.strings.pieces.size() == m_num_sequences + 1);

        if (m_build_config.weighted) {
            std::cout << "sum_of_weights " << m_sum_of_weights << std::endl;
            data.weights_builder.push_weight_interval(m_weight_value, m_weight_length);
            data.weights_builder.finalize(data.num_kmers);
        }
    }

    parse_data data;

private:
    /* general parameters */
    uint64_t m_k, m_m, m_seed, m_max_num_kmers_in_super_kmer, m_block_size;

    /* streaming variables */
    uint64_t m_prev_minimizer, m_begin, m_end, m_num_sequences, m_num_bases;
    bool m_glue;

    /* weights */
    uint64_t m_expected_sequence_length, m_sum_of_weights;
    uint64_t m_weight_value, m_weight_length;

    compact_string_pool::builder m_compact_string_pool_builder;
    build_configuration const& m_build_config;

    void append_super_kmer(char const* sequence, uint64_t length) {
        if (length == 0 or m_prev_minimizer == constants::invalid_uint64 or m_begin == m_end) {
            return;
        }

        assert(m_end > m_begin);
        char const* super_kmer = sequence + m_begin;
        uint64_t size = (m_end - m_begin) + m_k - 1;
        assert(util::is_valid(super_kmer, size));

        /* if num_kmers_in_super_kmer > m_k - m_m + 1, then split the super_kmer into blocks */
        uint64_t num_kmers_in_super_kmer = m_end - m_begin;
        uint64_t num_blocks = num_kmers_in_super_kmer / m_max_num_kmers_in_super_kmer +
                              (num_kmers_in_super_kmer % m_max_num_kmers_in_super_kmer != 0);
        assert(num_blocks > 0);
        for (uint64_t i = 0; i != num_blocks; ++i) {
            uint64_t n = m_block_size;
            if (i == num_blocks - 1) n = size;
            uint64_t num_kmers_in_block = n - m_k + 1;
            assert(num_kmers_in_block <= m_max_num_kmers_in_super_kmer);
            data.minimizers.emplace_back(m_prev_minimizer, m_compact_string_pool_builder.offset,
                                         num_kmers_in_block);
            m_compact_string_pool_builder.append(super_kmer + i * m_max_num_kmers_in_super_kmer, n,
                                                 m_glue);
            if (m_glue) {
                assert(data.minimizers.back().offset > m_k - 1);
                data.minimizers.back().offset -= m_k - 1;
            }
            size -= m_max_num_kmers_in_super_kmer;
            m_glue = true;
        }
    };
};

}  // namespace sshash