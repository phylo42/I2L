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
        /// IPK v0.1.x
        static const unsigned int v0_1_x = 2;

        /// IPK v0.2.x + v0.3.1
        static const unsigned int v0_2_WITHOUT_POSITIONS = 3;
        static const unsigned int v0_2_WITH_POSITIONS = v0_2_WITHOUT_POSITIONS + 1;

        /// IPK v0.4.0, added tree indexing
        static const unsigned int v0_4_0_WITHOUT_POSITIONS = 5;
        static const unsigned int v0_4_0_WITH_POSITIONS = v0_4_0_WITHOUT_POSITIONS + 1;

        /// IPK v0.4.1, added kmer ordering and filter values
        static const unsigned int v0_4_1_WITHOUT_POSITIONS = 7;
        static const unsigned int v0_4_1_WITH_POSITIONS = v0_4_1_WITHOUT_POSITIONS + 1;

#ifdef KEEP_POSITIONS
        static const unsigned int EARLIEST_INDEX = v0_4_0_WITH_POSITIONS;
        static const unsigned int CURRENT = v0_4_1_WITH_POSITIONS;
#else
        static const unsigned int EARLIEST_INDEX = v0_4_0_WITHOUT_POSITIONS;
        static const unsigned int CURRENT = v0_4_1_WITHOUT_POSITIONS;
        static const unsigned int CURRENT_WITH_POSITIONS = CURRENT + 1;

#endif
    };

    i2l::phylo_kmer_db load_compressed(const std::string& filename, float user_mu, float user_omega)
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

    i2l::phylo_kmer_db load_uncompressed(const std::string& filename, float user_mu, float user_omega)
    {
        std::ifstream ifs(filename);
        boost::archive::binary_iarchive ia(ifs);

        i2l::phylo_kmer_db db { 0, user_omega, "", "" };
        db.set_mu(user_mu);

        ia & db;
        return db;
    }


    i2l::phylo_kmer_db load(const std::string& filename, float mu, float user_epsilon)
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
            return load_compressed(filename, mu, user_epsilon);
        }
        catch (const boost::iostreams::zlib_error& error)
        {
            return load_uncompressed(filename, mu, user_epsilon);
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
    template<class T, class K>
    void found_or_die(const T& optional, const K& key)
    {
        if (!optional)
        {
            throw std::runtime_error("Error while serializing " + std::to_string(key) + ". K-mer not found");
        }
    }

    template<class Archive>
    void save_entry(Archive& ar, const i2l::positioned_pkdb_value& entry)
    {
        const auto& [branch, score, position] = entry;
        ar & branch & score & position;
    }

    template<class Archive>
    void save_entry(Archive& ar, const i2l::unpositioned_pkdb_value& entry)
    {
        const auto& [branch, score] = entry;
        ar & branch & score;
    }

    /// Saves k-mers with their filter values followed by vectors of scores.
    /// Introduced in v0.4.1
    template<class Archive, class Database>
    inline void save_filter_ordered_phylo_kmers(Archive& ar, const Database& db)
    {
        /// Phylo-k-mers order by Mutual Information
        for (const auto& [ key, filter_value ] : db.kmer_order)
        {
            const auto entries = db.search(key);
            found_or_die(entries, key);

            size_t entries_size = entries->size();
            ar & key & entries_size & filter_value;

            for (const auto& entry : *entries)
            {
                save_entry(ar, entry);
            }
        }
    }

    /// Saves k-mers with their vectors of scores. Pre-v0.4.1
    template<class Archive, class Database>
    inline void save_phylo_kmers(Archive& ar, const Database& db)
    {
        /// Phylo-k-mers order by Mutual Information
        for (const auto& [ key, filter_value ] : db.kmer_order)
        {
            const auto entries = db.search(key);
            found_or_die(entries, key);

            size_t entries_size = entries->size();
            ar & key & entries_size & filter_value;

            for (const auto& entry : *entries)
            {
                save_entry(ar, entry);
            }
        }
    }

    /// Type-agnostic save. Used for both positioned and unpositioned databases
    template<class Archive, class Database>
    inline void save(Archive& ar, const Database& db, unsigned int /*version*/)
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

        size_t num_kmers = db.size();
        ar & num_kmers;

        size_t num_entries = 0;
        ar & num_entries;

        if (!db.kmer_order.empty())
        {
            save_filter_ordered_phylo_kmers(ar, db);
        }
        else
        {
            save_phylo_kmers(ar, db);
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
    void load_entry(Archive& ar, i2l::positioned_pkdb_value& entry)
    {
        ar & entry.branch & entry.score & entry.position;
    }

    template<class Archive>
    void load_entry(Archive& ar, i2l::unpositioned_pkdb_value& entry)
    {
        ar & entry.branch & entry.score;
    }

    template <class T>
    T na_phylo_kmer();

    template<>
    i2l::positioned_pkdb_value na_phylo_kmer()
    {
        return { i2l::phylo_kmer::na_branch, i2l::phylo_kmer::na_score, i2l::phylo_kmer::na_pos };
    }

    template<>
    i2l::unpositioned_pkdb_value na_phylo_kmer()
    {
        return { i2l::phylo_kmer::na_branch, i2l::phylo_kmer::na_score };
    }

    template<class Archive, class Database>
    inline void load_db(Archive& ar, Database& db, const unsigned int version)
    {
        db.set_version(version);

        /// Early versions are not supported
        if (version < i2l::protocol::v0_1_x)
        {
            throw std::runtime_error("Failed to load database: this database was built with older version of IPK.");
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
        if (version >= i2l::protocol::v0_4_0_WITHOUT_POSITIONS)
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

        /// Ignore omega if we want to dynamic load
        if (version < i2l::protocol::v0_4_1_WITHOUT_POSITIONS)
        {
            db.set_omega(omega);
        }

        size_t num_kmers = 0;
        ar & num_kmers;

        size_t num_entries = 0;
        ar & num_entries;

        if (version < i2l::protocol::v0_4_1_WITHOUT_POSITIONS)
        {
            /// v0.4.0 and earlier: load all unordered k-mers
            for (size_t i = 0; i < num_kmers; ++i)
            {
                auto key = i2l::phylo_kmer::na_key;
                size_t entries_size = 0;
                float filter_value = 0.0f;

                ar & key;
                ar & entries_size;
                ar & filter_value;

                for (size_t j = 0; j < entries_size; ++j)
                {
                    /// classic deserialization of non-positioned phylo k-mers
                    auto entry = na_phylo_kmer<typename Database::pkdb_value_type>();
                    load_entry(ar, entry);
                    db.unsafe_insert(key, entry);
                }
            }
        }
        else
        {
            const auto user_threshold = std::log10(i2l::score_threshold(db.omega(), db.kmer_size()));
            const auto user_mu = db.get_mu();
            double fv_sum = 0.0f;
            const auto max_entries_allowed = user_mu * num_entries;

            /// Load ordered k-mers with dynamic mu and omega
            for (size_t i = 0; i < num_kmers; ++i)
            {
                auto key = i2l::phylo_kmer::na_key;
                size_t entries_size = 0;
                float filter_value;
                ar & key;
                ar & entries_size;
                ar & filter_value;

                size_t entries_added = 0;
                for (size_t j = 0; j < entries_size; ++j)
                {
                    /// classic deserialization of non-positioned phylo k-mers
                    auto entry = na_phylo_kmer<typename Database::pkdb_value_type>();
                    load_entry(ar, entry);

                    if (entry.score >= user_threshold)
                    {
                        db.unsafe_insert(key, entry);

                        if (entries_added == 0)
                        {
                            fv_sum += filter_value;
                        }
                        entries_added++;
                    }
                }

                // stop reading the database
                if (entries_added >= max_entries_allowed)
                {
                    break;
                }
            }
        }
    }

    template<class Archive>
    inline void load(Archive& ar,
                     i2l::_phylo_kmer_db<i2l::unpositioned_phylo_kmer>& db,
                     const unsigned int version)
    {
        load_db(ar, db, version);
    }

    template<class Archive>
    inline void load(Archive& ar,
                     i2l::_phylo_kmer_db<i2l::positioned_phylo_kmer>& db,
                     const unsigned int version)
    {
        load_db(ar, db, version);
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
