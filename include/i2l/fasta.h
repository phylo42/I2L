#ifndef I2L_FASTA_H
#define I2L_FASTA_H


//#include <boost/iostreams/device/mapped_file.hpp>
//#include <boost/iostreams/stream.hpp>
//#include <boost/algorithm/string/predicate.hpp>
#include <fstream>
#include <string>
#include <vector>
#include <i2l/seq_record.h>

namespace i2l::io
{
    class batch_fasta;

    namespace impl
    {
        //namespace bio = boost::iostreams;

        /// An efficient iterator over fasta files. Implements batched reading.
        /// Returns references to sequences stored in the iterator.
        class fasta_iterator
        {
            friend class i2l::io::batch_fasta;
        public:
            using value_type = const i2l::seq_record&;

            fasta_iterator(const std::string& filename, size_t batch_size, bool clean_sequences=true);
            fasta_iterator(const fasta_iterator&) = delete;
            fasta_iterator(fasta_iterator&&) = default;
            fasta_iterator& operator=(const fasta_iterator&) = delete;
            fasta_iterator& operator=(fasta_iterator&&) = delete;
            ~fasta_iterator() noexcept = default;

            fasta_iterator& operator++();

            bool operator==(const fasta_iterator& rhs) const noexcept;
            bool operator!=(const fasta_iterator& rhs) const noexcept;

            value_type operator*() const noexcept;
        private:
            void _read_batch();

            //bio::mapped_file_source _mmap;
            //bio::stream<bio::mapped_file_source> _is;
            std::ifstream _is;

            /// current batch of sequences
            std::vector<seq_record> _seqs;

            size_t _batch_size;
            /// the index of the current sequence in the batch
            size_t _seq_id;
            /// the global index of the current sequence
            size_t _global_seq_id;
            bool _last_batch;

            std::string _header;

            /// Indicates if we need to clean sequences on-the-fly with i2l::clean_sequence
            bool _clean_sequences;
        };
    }

    /// Efficiently iterates over sequences of a fasta file. Unaligns sequences if needed
    /// Usage:
    ///     for (const auto& seq: i2l::io::read_fasta(filename)) { ... }
    class read_fasta
    {
        using const_iterator = impl::fasta_iterator;
    public:
        explicit read_fasta(std::string filename, bool clean_sequences=true, size_t batch_size=1024);

        [[nodiscard]]
        const_iterator begin() const;
        [[nodiscard]]
        const_iterator end() const;

    private:
        std::string _filename;
        size_t _batch_size;
        bool _clean_sequences;
    };

    /// Efficiently reads batches of fasta sequences. Unaligns sequences if needed
    class batch_fasta
    {
        using const_iterator = impl::fasta_iterator;
    public:
        batch_fasta(const std::string& filename, size_t batch_size = 1024, bool clean_sequences = true);

        std::vector<seq_record> next_batch();
    private:
        impl::fasta_iterator _it;
    };

    /// Clean an input sequence from gaps. Characters '*', '!', '.' are also skipped
    /// and ignored
    std::string unalign(std::string sequence);

    std::string remove_whitespaces(std::string sequence);
}


/// \brief Outputs a collection of fasta records
std::ostream& operator<<(std::ostream& out, const std::vector<i2l::seq_record>& sequences);

#endif
