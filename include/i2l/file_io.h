#ifndef I2L_FILE_IO_H
#define I2L_FILE_IO_H

#include <string>
#include <memory>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>

namespace i2l::io
{
    /// \brief A buffered reader class for memory mapped files.
    class buffered_reader
    {
    public:
        explicit buffered_reader(const std::string& file_name);
        buffered_reader(const buffered_reader&) = delete;
        buffered_reader(buffered_reader&&) = delete;
        ~buffered_reader();

        buffered_reader& operator=(const buffered_reader&) = delete;
        buffered_reader& operator=(buffered_reader&&) = delete;

        std::string_view read_next_chunk();
        bool empty() const;
        bool good() const;

        static constexpr size_t buffer_size = 4096;

    private:
        std::fpos<mbstate_t> _get_file_legth();
        void _start_reading();
        void _read_next_chunk();

    private:
        using mf_source_t = boost::iostreams::mapped_file_source;
        using stream_t = boost::iostreams::stream<mf_source_t>;

        mf_source_t _msource;
        stream_t _stream;

        bool _started;
        std::fpos<mbstate_t> _file_length;
        std::fpos<mbstate_t> _read = 0;


        char _buffer[buffer_size];
    };

    /// \brief Reads the whole file in a string
    std::string read_as_string(const std::string& filename);
}

#endif