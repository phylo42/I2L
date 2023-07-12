#include <iostream>
#include <fstream>
#include <boost/iostreams/filter/zlib.hpp>
#include <i2l/serialization.h>

i2l::phylo_kmer_db i2l::load_compressed(const std::string& filename, float user_mu, float user_omega)
{
    std::ifstream ifs(filename);
    boost::iostreams::filtering_istream in;
    boost::iostreams::zlib_params zp(boost::iostreams::zlib::best_speed);
    in.push(boost::iostreams::zlib_decompressor(zp));
    in.push(ifs);

    boost::archive::binary_iarchive ia(in);

    i2l::phylo_kmer_db db { 0, user_omega, "", "" };
    db.set_mu(user_mu);

    ia & db;
    return db;
}

i2l::phylo_kmer_db i2l::load_uncompressed(const std::string& filename, float user_mu, float user_omega)
{
    std::ifstream ifs(filename);
    boost::archive::binary_iarchive ia(ifs);

    i2l::phylo_kmer_db db { 0, user_omega, "", "" };
    db.set_mu(user_mu);

    ia & db;
    return db;
}


i2l::phylo_kmer_db i2l::load(const std::string& filename, float mu, float user_epsilon)
{
    std::cout << "Boost version: " << BOOST_VERSION / 100000 << "."
              << BOOST_VERSION / 100 % 1000 << "."
              << BOOST_VERSION % 100 << std::endl;

    if (!fs::exists(filename))
    {
        throw std::runtime_error("No such file: " + filename);
    }

    /// Versions earlier than v0.2.1 were not using zlib compression.
    /// There is no good way to figure out if the file is compressed or not
    /// than to just try to decompress first.
    try
    {
        return i2l::load_compressed(filename, mu, user_epsilon);
    }
    catch (const boost::iostreams::zlib_error& error)
    {
        return i2l::load_uncompressed(filename, mu, user_epsilon);
    }
}

/// Archive-agnostic serialization of the database: write the header and the content.
/// Works both for compressed and uncompressed formats
template<class Archive>
void save_db(Archive& ar, const i2l::phylo_kmer_db& db)
{
    i2l::ipk_header header = {
        db.sequence_type(),
        db.tree_index(),
        db.tree(),
        db.kmer_size(),
        db.omega(),
        db.size(),
        get_num_entries(db)
    };
    ar & header;

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
    std::cout << "Boost version: " << BOOST_VERSION / 100000 << "."
              << BOOST_VERSION / 100 % 1000 << "."
              << BOOST_VERSION % 100 << std::endl;

    if (uncompressed)
    {
        i2l::save_uncompressed(db, filename);
    }
    else
    {
        i2l::save_compressed(db, filename);
    }
}
