#include <iostream>
#include <boost/algorithm/string/predicate.hpp>
#include <i2l/fasta.h>

using namespace i2l::io;
using std::vector;
using std::string, std::string_view;
using std::cout, std::endl;
using std::move;


//------------------------------------------------------------------------------------

impl::fasta_iterator::fasta_iterator(const std::string& filename, size_t batch_size, bool clean_sequences)
    : _is{ filename, std::ios::in }
    , _batch_size{ batch_size }
    , _seq_id{ 0 }
    , _global_seq_id{ 0 }
    , _last_batch{ false }
    , _clean_sequences{ clean_sequences }
    , _bytes_read(0)
{
    /// end() iterator
    if (_batch_size == 0)
    {
        _global_seq_id = 0;
    }
    else
    {
        _read_batch();
    }
}

impl::fasta_iterator& impl::fasta_iterator::operator++()
{
    ++_seq_id;
    ++_global_seq_id;

    /// if end of file, make it an end iterator
    if (_last_batch && _seq_id == _batch_size)
    {
        _global_seq_id = 0;
        _batch_size = 0;
    }
    else
    {
        if (_seq_id == _batch_size)
        {
            _read_batch();
        }
    }
    return *this;
}

bool impl::fasta_iterator::operator==(const fasta_iterator& rhs) const noexcept
{
    return _global_seq_id == rhs._global_seq_id && _batch_size == rhs._batch_size;
}

bool impl::fasta_iterator::operator!=(const fasta_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

impl::fasta_iterator::value_type impl::fasta_iterator::operator*() const noexcept
{
    return _seqs[_seq_id];
}

void impl::fasta_iterator::_read_batch()
{
    if (_last_batch)
    {
        _seqs.clear();
        return;
    }

    /// start iterating from the first sequence of the new batch
    _seq_id = 0;

    /// clean the old batch
    _seqs.clear();
    _seqs.reserve(_batch_size);

    /// read sequences until eof or batch limit
    string line, sequence;
    size_t i = 0;
    while ((i < _batch_size) && std::getline(_is, line))
    {
        _bytes_read += line.size();
        if (boost::starts_with(line, ">"))
        {
            if (!sequence.empty())
            {
                ++i;
                if (_clean_sequences)
                {
                    _seqs.emplace_back(_header, unalign(sequence));
                }
                else
                {
                    _seqs.emplace_back(_header, remove_whitespaces(sequence));
                }

                sequence = "";
                sequence.reserve(1024);
            }

            _header = {line.c_str() + 1, line.size() - 1};
        }
        else if (!line.empty())
        {
            sequence.append(line);
        }
    }

    /// if eof, it was the last batch
    if (i < _batch_size)
    {

        /// do not forget the last sequence
        if (_clean_sequences)
        {
            _seqs.emplace_back(_header, unalign(sequence));
        }
        else
        {
            _seqs.emplace_back(_header, remove_whitespaces(sequence));
        }
        _last_batch = true;
        _batch_size = i + 1;
    }
}

read_fasta::read_fasta(std::string filename, bool clean_sequences, size_t batch_size)
    : _filename{ std::move(filename) }
    , _batch_size{ batch_size }
    , _clean_sequences{ clean_sequences }
{}

read_fasta::const_iterator read_fasta::begin() const
{
    return { _filename, _batch_size, _clean_sequences };
}

read_fasta::const_iterator read_fasta::end() const
{
    return { _filename, 0, false };
}

batch_fasta::batch_fasta(const std::string& filename, size_t batch_size, bool clean_sequences)
    : _it(filename, batch_size, clean_sequences)
{}

std::vector<i2l::seq_record> batch_fasta::next_batch()
{
    auto batch = std::move(_it._seqs);
    _it._read_batch();
    return batch;
}

size_t batch_fasta::bytes_read() const
{
    return _it._bytes_read;
}

string i2l::io::remove_whitespaces(string sequence)
{
    sequence.erase(
        std::remove_if(sequence.begin(), sequence.end(), [](unsigned char x){ return x == ' '; }), sequence.end());
    return sequence;
}

string i2l::io::unalign(string sequence)
{
    sequence.erase(
        std::remove_if(sequence.begin(), sequence.end(),
                       [](unsigned char x){ return x == '-' || x == '.' || x == '!' || x == '*' || x == ' '; }),
        sequence.end());
    return sequence;
}

std::ostream& operator<<(std::ostream& out, const std::vector<i2l::seq_record>& sequences)
{
    for (const auto& seq : sequences)
    {
        out << seq.header() << '\n';
        out << seq.sequence() << '\n';
    }
    return out;
}