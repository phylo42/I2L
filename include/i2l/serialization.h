#ifndef I2L_SERIALIZATION_H
#define I2L_SERIALIZATION_H

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include "phylo_kmer_db.h"
#include "phylo_node.h"
#include "version.h"

namespace fs = boost::filesystem;

namespace i2l
{
    /// A struct with definitions for different version of the serialization protocol
    struct protocol
    {
        /// xpas v0.1.x
        static const unsigned int v0_1_x = 2;

        /// xpas v0.2.x + v0.3.1
        static const unsigned int v0_2_WITHOUT_POSITIONS = 3;
        static const unsigned int v0_2_WITH_POSITIONS = 4;

        /// xpas v0.3.2, added tree indexing
        static const unsigned int v0_3_2_WITHOUT_POSITIONS = 5;
        static const unsigned int v0_3_2_WITH_POSITIONS = 6;

#ifdef KEEP_POSITIONS
        static const unsigned int EARLIEST_INDEX = v0_3_2_WITH_POSITIONS;
        static const unsigned int CURRENT = v0_3_2_WITH_POSITIONS;
#else
        static const unsigned int EARLIEST_INDEX = v0_3_2_WITHOUT_POSITIONS;
        static const unsigned int CURRENT = v0_3_2_WITHOUT_POSITIONS;
#endif
    };

    i2l::phylo_kmer_db load_compressed(const std::string& filename)
    {
        std::ifstream ifs(filename);
        boost::iostreams::filtering_istream in;
        boost::iostreams::zlib_params zp(boost::iostreams::zlib::best_speed);
        in.push(boost::iostreams::zlib_decompressor(zp));
        in.push(ifs);

        boost::archive::binary_iarchive ia(in);

        i2l::phylo_kmer_db db { 0, 0.0, "", "" };
        ia & db;
        return db;
    }

    i2l::phylo_kmer_db load_uncompressed(const std::string& filename)
    {
        std::ifstream ifs(filename);
        boost::archive::binary_iarchive ia(ifs);

        i2l::phylo_kmer_db db { 0, 0.0, "", "" };
        ia & db;
        return db;
    }


    i2l::phylo_kmer_db load(const std::string& filename)
    {
        if (!fs::exists(filename))
        {
            throw std::runtime_error("No such file: " + filename);
        }

        /// Versions earlier than v0.2.1 were not using zlib compression.
        /// There is no good way to figure out if the file is compressed or not
        /// than to just try to decompress first.
        try
        {
            return load_compressed(filename);
        }
        catch (const boost::iostreams::zlib_error& error)
        {
            return load_uncompressed(filename);
        }
    }

    void save_compressed(const i2l::phylo_kmer_db& db, const std::string& filename)
    {
        std::ofstream ofs(filename);

        boost::iostreams::filtering_ostream out;
        boost::iostreams::zlib_params zp(boost::iostreams::zlib::best_speed);
        out.push(boost::iostreams::zlib_compressor(zp));
        out.push(ofs);

        ::boost::archive::binary_oarchive oa(out);
        oa & db;
    }

    void save_uncompressed(const i2l::phylo_kmer_db& db, const std::string& filename)
    {
        std::ofstream out(filename);
        ::boost::archive::binary_oarchive oa(out);
        oa & db;
    }

    void save(const i2l::phylo_kmer_db& db, const std::string& filename, bool uncompressed=false)
    {
        if (uncompressed)
        {
            save_uncompressed(db, filename);
        }
        else
        {
            save_compressed(db, filename);
        }
    }
}

namespace boost::serialization
{
    template<class Archive>
    inline void save(Archive& ar,
                     const i2l::_phylo_kmer_db<i2l::positioned_phylo_kmer>& db,
                     unsigned int /*version*/)
    {
        ar & std::string(db.sequence_type());

        /// Tree index
        const auto& tree_index = db.tree_index();
        ar & tree_index.size();
        for (const auto& x : tree_index)
        {
            ar & x.subtree_num_nodes & x.subtree_total_length;
        }

        const auto original_tree_view = std::string{ db.tree() };
        ar & original_tree_view;

        const auto kmer_size = db.kmer_size();
        ar & kmer_size;

        i2l::phylo_kmer::score_type omega = db.omega();
        ar & omega;

        size_t table_size = db.size();
        ar & table_size;

        for (const auto& [key, entries] : db)
        {
            size_t entries_size = entries.size();
            ar & key & entries_size;

            for (const auto& [branch, score, position] : entries)
            {
                ar & branch & score & position;
            }
        }
    }

    template<class Archive>
    inline void save(Archive& ar,
                     const i2l::_phylo_kmer_db<i2l::unpositioned_phylo_kmer>& db,
                     unsigned int /*version*/)
    {
        ar & std::string(db.sequence_type());

        /// Tree index
        const auto& tree_index = db.tree_index();
        ar & tree_index.size();
        for (const auto& x : tree_index)
        {
            ar & x.subtree_num_nodes & x.subtree_total_length;
        }

        /// Tree newick
        const auto original_tree_view = std::string{ db.tree() };
        ar & original_tree_view;

        /// DB construction parameters
        size_t kmer_size = db.kmer_size();
        ar & kmer_size;

        i2l::phylo_kmer::score_type omega = db.omega();
        ar & omega;

        /// The number of different k-mers
        size_t table_size = db.size();
        ar & table_size;

        /// Phylo-k-mers
        for (const auto& [key, entries] : db)
        {
            size_t entries_size = entries.size();
            ar & key & entries_size;

            for (const auto& [branch, score] : entries)
            {
                ar & branch & score;
            }
        }
    }

    /// Loads the tree index: the number of the nodes in the subtree for every node,
    ///     the total branch length in the subtree
    template<class Archive>
    std::vector<i2l::phylo_node::node_index> load_tree_index(Archive& ar)
    {
        size_t num_nodes;
        ar & num_nodes;

        std::vector<i2l::phylo_node::node_index> index(num_nodes);
        for (size_t i = 0; i < num_nodes; ++i)
        {
            ar & index[i].subtree_num_nodes;
            ar & index[i].subtree_total_length;
        }
        return index;
    }


    template<class Archive>
    inline void load(Archive& ar,
                     i2l::_phylo_kmer_db<i2l::positioned_phylo_kmer>& db,
                     const unsigned int version)
    {
        db.set_version(version);

        /// Early versions are not supported
        if (version < i2l::protocol::v0_2_WITH_POSITIONS)
        {
            throw std::runtime_error("Failed to load database: the database does not have positional information.");
        }

        /// if we did not throw an exception, positions will be loaded.
        db.set_positions_loaded(true);

        /// Deserialization of the content added in versions v0.2.x
        if (version > i2l::protocol::v0_1_x)
        {
            std::string sequence_type;
            ar & sequence_type;
            db.set_sequence_type(std::move(sequence_type));
        }

        /// Deserialization of the tree index, v0.3.2 and later
        if (version >= i2l::protocol::v0_3_2_WITHOUT_POSITIONS)
        {
            std::vector<i2l::phylo_node::node_index> tree_index = load_tree_index(ar);
            db.set_tree_index(std::move(tree_index));
        }

        /// Deserialization of the main content
        {
            std::string tree;
            ar & tree;
            db.set_tree(std::move(tree));

            size_t kmer_size = 0;
            ar & kmer_size;
            db.set_kmer_size(kmer_size);

            i2l::phylo_kmer::score_type omega = 0;
            ar & omega;
            db.set_omega(omega);

            size_t table_size = 0;
            ar & table_size;
            for (size_t i = 0; i < table_size; ++i)
            {
                auto key = i2l::phylo_kmer::na_key;
                size_t entries_size = 0;
                ar & key;
                ar & entries_size;
                for (size_t j = 0; j < entries_size; ++j)
                {
                    auto branch = i2l::phylo_kmer::na_branch;
                    auto score = i2l::phylo_kmer::na_score;
                    auto position = i2l::phylo_kmer::na_pos;
                    ar & branch & score & position;
                    db.unsafe_insert(key, { branch, score, position });
                }
            }
        }
    }

    template<class Archive>
    inline void load(Archive& ar,
                     i2l::_phylo_kmer_db<i2l::unpositioned_phylo_kmer>& db,
                     const unsigned int version)
    {
        db.set_version(version);

        /// Early versions are not supported
        if (version < i2l::protocol::v0_1_x)
        {
            throw std::runtime_error("Failed to load database: this database was built with older version of xpas.");
        }

        db.set_positions_loaded(false);

        /// Deserialization of the content added in versions v0.2.x
        if (version > i2l::protocol::v0_1_x)
        {
            std::string sequence_type;
            ar & sequence_type;
            db.set_sequence_type(std::move(sequence_type));
        }

        /// Deserialization of the tree index, v0.3.2 and later
        if (version >= i2l::protocol::v0_3_2_WITHOUT_POSITIONS)
        {
            auto tree_index = load_tree_index(ar);
            db.set_tree_index(std::move(tree_index));
        }

        std::string tree;
        ar & tree;
        db.set_tree(std::move(tree));

        size_t kmer_size = 0;
        ar & kmer_size;
        db.set_kmer_size(kmer_size);

        i2l::phylo_kmer::score_type omega = 0;
        ar & omega;
        db.set_omega(omega);

        size_t table_size = 0;
        ar & table_size;
        for (size_t i = 0; i < table_size; ++i)
        {
            auto key = i2l::phylo_kmer::na_key;
            size_t entries_size = 0;
            ar & key;
            ar & entries_size;
            for (size_t j = 0; j < entries_size; ++j)
            {
                /// classic deserialization of non-positioned phylo k-mers
                auto branch = i2l::phylo_kmer::na_branch;
                auto score = i2l::phylo_kmer::na_score;
                ar & branch & score;

                /// if the database has positions, read and ignore them
                if (version == i2l::protocol::v0_2_WITH_POSITIONS)
                {
                    auto position = i2l::phylo_kmer::na_pos;
                    ar & position;
                }
                db.unsafe_insert(key, { branch, score });
            }
        }
    }

    // split non-intrusive serialization function member into separate
    // non intrusive save/load member functions
    template<class Archive>
    inline void serialize(Archive& ar, i2l::phylo_kmer_db& db, unsigned int file_version)
    {
        boost::serialization::split_free(ar, db, file_version);
    }
}

BOOST_CLASS_VERSION(i2l::phylo_kmer_db, i2l::protocol::CURRENT)

#endif
