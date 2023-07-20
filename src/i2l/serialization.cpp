#include <iostream>
#include <fstream>
#include <boost/iostreams/filter/zlib.hpp>
#include <i2l/serialization.h>


void print_boost_version()
{
    std::cout << "Boost version: " << BOOST_VERSION / 100000 << "."
              << BOOST_VERSION / 100 % 1000 << "."
              << BOOST_VERSION % 100 << std::endl;
}

namespace i2l::legacy
{
    i2l::phylo_kmer_db load_compressed(const std::string& filename)
    {
        std::ifstream ifs(filename);
        boost::iostreams::filtering_istream in;
        boost::iostreams::zlib_params zp(boost::iostreams::zlib::best_speed);
        in.push(boost::iostreams::zlib_decompressor(zp));
        in.push(ifs);

        boost::archive::binary_iarchive ar(in);

        i2l::phylo_kmer_db db { 0, 0.0, "", "" };
        ar & db;
        return db;
    }

    i2l::phylo_kmer_db load_uncompressed(const std::string& filename)
    {
        std::ifstream ifs(filename);
        boost::archive::binary_iarchive ar(ifs);

        i2l::phylo_kmer_db db { 0, 0.0, "", "" };
        ar & db;
        return db;
    }

    i2l::phylo_kmer_db load_pre_v4_1(const std::string& filename)
    {
        if (!fs::exists(filename))
        {
            throw std::runtime_error("No such file: " + filename);
        }

        try
        {
            return load_compressed(filename);
        }
        catch (const boost::iostreams::zlib_error& error)
        {
            return load_uncompressed(filename);
        }
    }
}



i2l::phylo_kmer_db i2l::load_compressed(const std::string& filename, float user_mu,
                                        float user_omega, size_t max_entries)
{
    std::ifstream ifs(filename);
    boost::iostreams::filtering_istream in;
    boost::iostreams::zlib_params zp(boost::iostreams::zlib::best_speed);
    in.push(boost::iostreams::zlib_decompressor(zp));
    in.push(ifs);

    boost::archive::binary_iarchive ar(in);

    i2l::phylo_kmer_db db { 0, user_omega, "", "", user_mu, max_entries };
    unsigned int version = i2l::protocol::ERROR;
    ar & version;
    db.set_version(version);
    load_db(ar, db, version);
    return db;
}

i2l::phylo_kmer_db i2l::load_uncompressed(const std::string& filename, float user_mu,
                                          float user_omega, size_t max_entries)
{
    std::ifstream ifs(filename);
    boost::archive::binary_iarchive ar(ifs);

    i2l::phylo_kmer_db db { 0, user_omega, "", "", user_mu, max_entries };
    unsigned int version = i2l::protocol::ERROR;
    ar & version;
    db.set_version(version);
    load_db(ar, db, version);
    return db;
}

i2l::phylo_kmer_db i2l::load(const std::string& filename, float mu, float user_epsilon, size_t max_entries)
{
    print_boost_version();

    if (!fs::exists(filename))
    {
        throw std::runtime_error("No such file: " + filename);
    }

    /// Versions earlier than v0.2.1 were not using zlib compression.
    /// There is no good way to figure out if the file is compressed or not
    /// than to just try to decompress first.
    try
    {
        return i2l::load_compressed(filename, mu, user_epsilon, max_entries);
    }
    catch (const boost::iostreams::zlib_error& error)
    {
        return i2l::load_uncompressed(filename, mu, user_epsilon, max_entries);
    }
    /// std::bad_alloc is thrown as format error if we try to
    /// load old databases. Use legacy deserialization in this case
    catch (const std::bad_alloc&)
    {
        return legacy::load_pre_v4_1(filename);
    }
}

/// Archive-agnostic serialization of the database: write the header and the content.
/// Works both for compressed and uncompressed formats
template<class Archive>
void save_db(Archive& ar, const i2l::phylo_kmer_db& db)
{
    /// Serialize the protocol header
    const i2l::ipk_header header = {
        db.sequence_type(),
        db.tree_index(),
        db.tree(),
        db.kmer_size(),
        db.omega(),
        db.size(),
        get_num_entries(db)
    };
    save_header(ar, header);

    /// Serialize content
    if (!db.kmer_order.empty())
    {
        save_filter_ordered_phylo_kmers(ar, db);
    }
    else
    {
        save_phylo_kmers(ar, db);
    }
}

void i2l::save_compressed(const i2l::phylo_kmer_db& db, const std::string& filename)
{
    std::ofstream ofs(filename);

    boost::iostreams::filtering_ostream out;
    boost::iostreams::zlib_params zp(boost::iostreams::zlib::best_speed);
    out.push(boost::iostreams::zlib_compressor(zp));
    out.push(ofs);

    ::boost::archive::binary_oarchive ar(out);
    save_db(ar, db);
}

void i2l::save_uncompressed(const i2l::phylo_kmer_db& db, const std::string& filename)
{
    std::ofstream out(filename);
    ::boost::archive::binary_oarchive ar(out);
    save_db(ar, db);
}

void i2l::save(const i2l::phylo_kmer_db& db, const std::string& filename, bool uncompressed)
{
    print_boost_version();

    if (uncompressed)
    {
        i2l::save_uncompressed(db, filename);
    }
    else
    {
        i2l::save_compressed(db, filename);
    }
}
