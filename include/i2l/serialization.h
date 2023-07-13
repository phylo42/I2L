#ifndef I2L_SERIALIZATION_H
#define I2L_SERIALIZATION_H

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/filesystem.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <utility>
#include <fstream>
#include <memory>
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
        static const unsigned int ERROR = 42;
    };

    i2l::phylo_kmer_db load(const std::string& filename, float mu=1.0, float user_epsilon=0.0f);
    i2l::phylo_kmer_db load_compressed(const std::string& filename, float user_mu, float user_omega);
    i2l::phylo_kmer_db load_uncompressed(const std::string& filename, float user_mu, float user_omega);

    void save(const i2l::phylo_kmer_db& db, const std::string& filename, bool uncompressed = false);
    void save_compressed(const i2l::phylo_kmer_db& db, const std::string& filename);
    void save_uncompressed(const i2l::phylo_kmer_db& db, const std::string& filename);

    using oarchive = boost::archive::binary_oarchive;

    /// Tree index: numbers of the nodes in the subtree for every node,
    //  the total branch lengths of subtrees
    using tree_index = std::vector<i2l::phylo_node::node_index>;

    /// IPK binary format header
    struct ipk_header
    {
        /// DNA or Proteins
        std::string sequence_type;

        /// A number of tree nodes and a vector of every branch length
        tree_index index;

        /// Newick-formatted tree
        std::string tree;

        size_t kmer_size;

        i2l::phylo_kmer::score_type omega;

        /// Total number of distinct k-mers
        size_t num_kmers;

        /// Total number of branch-score pairs
        size_t num_entries;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int /*version*/)
        {
            ar & sequence_type;
            ar & index;
            ar & tree;
            ar & kmer_size;
            ar & omega;
            ar & num_kmers;
            ar & num_entries;
        }
    };

    /// One k-mer with its filter value and scores
    template<class Database>
    struct serialization_unit
    {
        i2l::phylo_kmer::key_type key;
        float filter_value;
        std::vector<typename Database::pkdb_value_type> entries;

        serialization_unit()
            : key(i2l::phylo_kmer::na_key), filter_value(0.0), entries()
        {}

        serialization_unit(i2l::phylo_kmer::key_type _key, float _filter_value,
                           std::vector<typename Database::pkdb_value_type>&& _entries)
            : key(_key), filter_value(_filter_value), entries(std::move(_entries))
        {}

        serialization_unit(const serialization_unit&) = delete;
        serialization_unit(serialization_unit&&) noexcept = default;
        serialization_unit& operator=(const serialization_unit&) = delete;
        serialization_unit& operator=(serialization_unit&&) noexcept = default;
        ~serialization_unit() = default;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int /*version*/)
        {
            ar & key;
            ar & filter_value;
            ar & entries;
        }

        [[nodiscard]]
        bool is_valid() const
        {
            return key != i2l::phylo_kmer::na_key;
        }
    };

    template<class T, class K>
    void found_or_die(const T& optional, const K& key)
    {
        if (!optional)
        {
            throw std::runtime_error("Error while serializing " + std::to_string(key) + ". K-mer not found");
        }
    }

    /// Deserializes one phylo-k-mer while the on-demand deserialization.
    /// Matches the protocol of save_filter_ordered_phylo_kmers (see below).
    /// Both are needed (one for on-demand deserialization, another for
    /// batched merge)
    template<class Archive>
    inline void save_phylo_kmer(Archive& ar, i2l::phylo_kmer::key_type key,
                                i2l::phylo_kmer::score_type filter_value,
                                const std::vector<pkdb_value>& entries)
    {
        size_t entries_size = entries.size();
        ar & key & filter_value & entries_size;
        for (const auto& entry: entries)
        {
            ar & entry;
        }
    }

    /// Saves k-mers with their filter values followed by vectors of scores.
    /// Introduced in v0.4.1
    template<class Archive, class Database>
    inline void save_filter_ordered_phylo_kmers(Archive& ar, const Database& db)
    {
        /// Phylo-k-mers ordered by Mutual Information
        for (const auto& [key, filter_value]: db.kmer_order)
        {
            const auto entries = db.search(key);
            found_or_die(entries, key);

            size_t entries_size = entries->size();
            ar & key & filter_value & entries_size;
            for (const auto& entry: *entries)
            {
                ar & entry;
            }
        }
    }

    template<class Archive>
    inline void save_header(Archive& ar, const i2l::ipk_header& header)
    {
        /// Serialize the protocol version. It is not saved by Boost because
        /// our serialization of i2l::phylo_kmer_db is not canonically defined
        unsigned int version = i2l::protocol::CURRENT;
        ar & version;

        /// Serialize the header
        ar & header;
    }

    /// Saves k-mers with their vectors of scores in a random order. Pre-v0.4.1
    template<class Archive, class Database>
    inline void save_phylo_kmers(Archive& ar, const Database& db)
    {
        for (const auto& [key, entries]: db)
        {
            size_t entries_size = entries.size();
            ar & key;
            ar & entries_size;
            for (const auto& entry: entries)
            {
                ar & entry;
            }
        }
    }

    template<class Database>
    inline size_t get_num_entries(const Database& db)
    {
        size_t num_entries = 0;
        for (const auto& [kmer, entries]: db)
        {
            (void) kmer;
            num_entries += entries.size();
        }
        return num_entries;
    }

    /// Type-agnostic save. Used for both positioned and unpositioned databases
    template<class Archive, class Database>
    inline void save(Archive& ar, const Database& db, unsigned int /*version*/)
    {
        ipk_header header = {
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

    namespace legacy
    {
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
        inline void load_content(Archive& ar,
                                 i2l::_phylo_kmer_db<i2l::unpositioned_phylo_kmer>& db,
                                 const unsigned int version, size_t table_size)
        {
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

        template<class Archive>
        inline void load_content(Archive& ar,
                                 i2l::_phylo_kmer_db<i2l::positioned_phylo_kmer>& db,
                                 const unsigned int version, size_t table_size)
        {
            (void)version;
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

        template<class Archive, class Database>
        inline void load(Archive& ar,
                         Database& db,
                         const unsigned int version)
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
            db.set_omega(omega);

            size_t table_size = 0;
            ar & table_size;

            load_content(ar, db, version, table_size);
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
    }

    template<class Archive, class Database>
    inline void load_db(Archive& ar, Database& db, unsigned int version)
    {
        /// Version is not written in the latest protocol. See save_db()
        if (version == 0)
        {
            version = protocol::CURRENT;
        }

        /// Prior to v0.4.1, phylo-k-mers are deserialized in random order
        if (version < protocol::v0_4_1_WITHOUT_POSITIONS)
        {
            legacy::load(ar, db, version);
        }
        else
        {
            ipk_header header;
            ar & header;

            db.set_sequence_type(header.sequence_type);
            db.set_tree_index(std::move(header.index));
            db.set_tree(std::move(header.tree));
            db.set_kmer_size(header.kmer_size);

            if (db.omega() == 0.0f)
            {
                db.set_omega(header.omega);
            }

            /// The total number of phylo-k-mer pairs loaded from disk
            size_t num_pk_loaded = 0;
            const auto user_threshold = std::log10(i2l::score_threshold(db.omega(), db.kmer_size()));
            const auto user_mu = db.get_mu();
            const auto max_pk_allowed = user_mu * header.num_entries;

            /// Load ordered k-mers with dynamic mu and omega
            for (size_t i = 0; i < header.num_kmers; ++i)
            {
                auto key = i2l::phylo_kmer::na_key;
                size_t entries_size = 0;
                float filter_value;
                ar & key & filter_value & entries_size;

                for (size_t j = 0; j < entries_size; ++j)
                {
                    /// classic deserialization of non-positioned phylo k-mers
                    auto entry = i2l::na_phylo_kmer<typename Database::pkdb_value_type>();
                    ar & entry;

                    if (entry.score >= user_threshold)
                    {
                        db.unsafe_insert(key, entry);
                        num_pk_loaded++;
                    }
                }

                // stop reading the database
                if (num_pk_loaded >= max_pk_allowed)
                {
                    break;
                }
            }
            db.set_num_entries_loaded(num_pk_loaded);
        }
    }

    /// Lazy deserializion of phylo-k-mer databases. Retrieves the next k-mer and the vector of
    /// its scores on demand. Keeps track of the last k-mer read.
    template<class Database>
    class lazy_load
    {
    public:
        explicit lazy_load(const std::string& filename)
            : _ifs(std::make_unique<std::ifstream>(filename))
            , _ar(std::make_unique<boost::archive::binary_iarchive>(*_ifs))
            , _db(0, 0.0, "", "")
            , _num_kmers_read(0)
            , _current{ i2l::phylo_kmer::na_key, 0.0, {} }
        {
            if (!_ifs->good())
            {
                throw std::runtime_error("Failed to open file: " + filename);
            }

            (*_ar) & _header;
        }
        lazy_load(const lazy_load& other) = delete;
        lazy_load operator=(const lazy_load& other) = delete;

        lazy_load(lazy_load&& other) noexcept
            : _ifs(std::move(other._ifs))
            , _ar(std::move(other._ar))
            , _db(std::move(other._db))
            , _header(std::move(other._header))
            , _num_kmers_read(std::move(other._num_kmers_read))
            , _current(std::move(other._current))
        {
        }

        lazy_load& operator=(lazy_load&& other) noexcept
        {
            _ifs = std::move(other._ifs);
            _ar = std::move(other._ar);
            _db = std::move(other._db);
            _header = std::move(other._header);
            _num_kmers_read = std::move(other._num_kmers_read);
            _current = std::move(other._current);
            return *this;
        }

        ~lazy_load() noexcept = default;

        [[nodiscard]]
        bool has_next() const
        {
            return _num_kmers_read < _header.num_kmers;
        }

        [[nodiscard]]
        serialization_unit<Database>& current()
        {
            return _current;
        }

        [[nodiscard]]
        const serialization_unit<Database>& current() const
        {
            return _current;
        }

        serialization_unit<Database>& next()
        {
            if (_ifs->eof() || !has_next())
            {
                _current = {};
            }
            else
            {
                auto key = i2l::phylo_kmer::na_key;
                size_t entries_size = 0;
                float filter_value;

                (*_ar) & key & filter_value & entries_size;

                std::vector<pkdb_value> entries;
                entries.reserve(entries_size);
                for (size_t j = 0; j < entries_size; ++j)
                {
                    /// classic deserialization of non-positioned phylo k-mers
                    auto entry = i2l::na_phylo_kmer<typename Database::pkdb_value_type>();
                    (*_ar) & entry;
                    entries.push_back(entry);
                }

                _current = serialization_unit<Database>(key, filter_value, std::move(entries));
                _num_kmers_read++;
            }

            return _current;
        }

    private:
        std::unique_ptr<std::ifstream> _ifs;
        std::unique_ptr<boost::archive::binary_iarchive> _ar;
        Database _db;
        ipk_header _header;
        size_t _num_kmers_read;
        serialization_unit<Database> _current;
    };

    using batch_loader = lazy_load<phylo_kmer_db>;

    /// "Less" comparator for batch loaders. Compares its top elements
    /// if they have any. If one element is empty, it is "greater"
    /// to push it away if used in the min-heap
    struct batch_loader_compare
    {
        bool operator()(const batch_loader* a, const batch_loader* b) const
        {
            const auto& a_current = a->current();
            const auto& b_current = b->current();

            if (!a_current.is_valid() && !b_current.is_valid())
            {
                return true;
            }
            // a goes before b as it has no current element
            if (!a_current.is_valid())
            {
                return true;
            }
            // b goes before a as it has no current element
            if (!b_current.is_valid())
            {
                return false;
            }
            return a_current.filter_value > b_current.filter_value;
        }
    };
}


namespace boost::serialization
{
    /// Non-intrusive (de-)serialization functions of other i2l types:

    /// A branch-score pair
    template<class Archive>
    void save(Archive & ar, const i2l::unpositioned_pkdb_value& t, const unsigned int /*version*/)
    {
        ar & t.branch;
        ar & t.score;
    }

    template<class Archive>
    void load(Archive & ar, i2l::unpositioned_pkdb_value& t, const unsigned int /*version*/)
    {
        ar & t.branch;
        ar & t.score;
    }

    /// A branch-score-position triple
    template<class Archive>
    void save(Archive & ar, const i2l::positioned_pkdb_value& t, const unsigned int /*version*/)
    {
        ar & t.branch;
        ar & t.score;
        ar & t.position;
    }

    template<class Archive>
    void load(Archive & ar, i2l::positioned_pkdb_value& t, const unsigned int /*version*/)
    {
        ar & t.branch;
        ar & t.score;
        ar & t.position;
    }

    /// Node index struct (num of nodes in the subtree, branch length)
    template<class Archive>
    void save(Archive & ar, const i2l::phylo_node::node_index& t, const unsigned int /*version*/)
    {
        ar & t.subtree_num_nodes;
        ar & t.subtree_total_length;
    }

    template<class Archive>
    void load(Archive & ar, i2l::phylo_node::node_index& t, const unsigned int /*version*/)
    {
        ar & t.subtree_num_nodes;
        ar & t.subtree_total_length;
    }

    template<class Archive>
    inline void load(Archive& ar,
                     i2l::_phylo_kmer_db<i2l::unpositioned_phylo_kmer>& db,
                     const unsigned int version)
    {
        i2l::load_db(ar, db, version);
    }

    template<class Archive>
    inline void load(Archive& ar,
                     i2l::_phylo_kmer_db<i2l::positioned_phylo_kmer>& db,
                     const unsigned int version)
    {
        i2l::load_db(ar, db, version);
    }

    // split non-intrusive serialization function member into separate
    // non intrusive save/load member functions
    template<class Archive>
    inline void serialize(Archive& ar, i2l::phylo_kmer_db& db, unsigned int file_version)
    {
        boost::serialization::split_free(ar, db, file_version);
    }
}

BOOST_SERIALIZATION_SPLIT_FREE(i2l::unpositioned_pkdb_value)
BOOST_SERIALIZATION_SPLIT_FREE(i2l::positioned_pkdb_value)
BOOST_SERIALIZATION_SPLIT_FREE(i2l::phylo_node::node_index)
BOOST_CLASS_VERSION(i2l::phylo_kmer_db, i2l::protocol::CURRENT)

#endif
