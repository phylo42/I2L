#ifndef I2L_PHYLO_KMER_DB_H
#define I2L_PHYLO_KMER_DB_H

#include "hash_map.h"
#include "phylo_kmer.h"
#include "optional.h"
#include "phylo_node.h"


namespace i2l
{
    namespace impl {
        template <typename PhyloKmer>
        class search_result;

        template <>
        class search_result<unpositioned_phylo_kmer>;
    }

    struct kmer_fv
    {
        i2l::phylo_kmer::key_type key;
        float filter_value;

        kmer_fv(i2l::phylo_kmer::key_type k, float fv)
            : key(k), filter_value(fv)
        {}
    };

    /// \brief A phylo-kmer database that stores all phylo-kmers.
    template <typename PhyloKmer>
    class _phylo_kmer_db
    {
    public:
        /// Member types
        using key_type = phylo_kmer::key_type;
        using pkdb_value_type = typename get_pkdb_value_type<PhyloKmer>::type;
        using value_type = std::vector<pkdb_value_type>;

        /// \brief A storage of a phylo-kmer information.
        /// \details Note that phylo-kmers are not stored as objects of phylo_kmer,
        /// which is just a temporary storage for a phylo-kmer information.
        using storage = hash_map<key_type, value_type>;
        using const_iterator = typename storage::const_iterator;

        /// Ctors, dtor and operator=
        _phylo_kmer_db(size_t kmer_size, phylo_kmer::score_type omega, std::string seq_type,
                       std::string tree)
            : _kmer_size{ kmer_size }
            , _omega{ omega }
            , _sequence_type(std::move(seq_type))
            , _tree(std::move(tree))
            , _version(0)
            , _mu(1)
        {}
        _phylo_kmer_db(const _phylo_kmer_db&) noexcept = delete;
        _phylo_kmer_db(_phylo_kmer_db&&) noexcept = default;
        _phylo_kmer_db& operator=(const _phylo_kmer_db&) = delete;
        _phylo_kmer_db& operator=(_phylo_kmer_db&&) noexcept = default;
        ~_phylo_kmer_db() noexcept = default;

        /// \brief Returns the sequence type: DNA or Protein
        [[nodiscard]]
        std::string_view sequence_type() const noexcept
        {
            return _sequence_type;
        }

        void set_sequence_type(std::string sequence_type)
        {
            _sequence_type = std::move(sequence_type);
        }

        /// \brief Returns the k-mer size.
        [[nodiscard]]
        size_t kmer_size() const noexcept
        {
            return _kmer_size;
        }

        void set_kmer_size(size_t kmer_size) noexcept
        {
            _kmer_size = kmer_size;
        }

        /// \brief Returns omega (the parameter of i2l::score_threshold)
        [[nodiscard]]
        phylo_kmer::score_type omega() const noexcept
        {
            return _omega;
        }

        void set_omega(phylo_kmer::score_type omega) noexcept
        {
            _omega = omega;
        }

        /// \brief Returns if the positional information was loaded during deserialization.
        [[nodiscard]]
        bool positions_loaded() const noexcept
        {
            return _positions_loaded;
        }

        void set_positions_loaded(bool value) noexcept
        {
            _positions_loaded = value;
        }

        /// \brief Returns a view to the newick formatted phylogenetic tree
        [[nodiscard]]
        std::string_view tree() const noexcept
        {
            return _tree;
        }

        void set_tree(std::string tree)
        {
            _tree = std::move(tree);
        }

        [[nodiscard]]
        std::vector<phylo_node::node_index> tree_index() const noexcept
        {
            return _tree_index;
        }

        [[nodiscard]]
        std::vector<phylo_node::node_index>& tree_index()
        {
            return _tree_index;
        }

        void set_tree_index(std::vector<phylo_node::node_index> index)
        {
            _tree_index = std::move(index);
        }


        /// Access
        /// \brief Searches for a key against the database.
        /// \details WARNING: This method does not know how the key was calculated. It is required
        /// to provide keys of substrings of size _kmer_size to get correct results.
        /// \sa _kmer_size

        [[nodiscard]]
        optional<impl::search_result<value_type>> search(key_type key) const noexcept
        {
            if (auto it = _map.find(key); it != _map.end())
            {
                return impl::search_result<value_type>{ it->second.begin(), it->second.end() };
            }
            else
            {
                return nullopt;
            }
        }
        /// Iterators
        /// \brief Returns an iterator to the beginning
        [[nodiscard]]
        const_iterator begin() const noexcept
        {
            return std::begin(_map);
        }

        /// \brief Returns an iterator to the end
        [[nodiscard]]
        const_iterator end() const noexcept
        {
            return std::end(_map);
        }

        /// Capacity
        /// \brief Returns the number of keys
        [[nodiscard]]
        size_t size() const noexcept
        {
            return _map.size();
        }

        /// \brief Returns a hash function used to hash kmers
        [[nodiscard]]
        typename storage::hasher hash_function() const noexcept
        {
            return _map.hash_function();
        }

        /// Modifiers
        /// \brief Puts a phylo-kmer information in the database.
        /// \details This method is unsafe, which means it does not control if the value has
        /// a correct branch id, the score is maximal etc. All of this checks must be done before calling
        /// this method. It just puts the value in a hash map.
        /// WARNING: This method does not know how the key was calculated. Here we assume it represents
        /// a string of size _kmer_size.
        /// \sa _kmer_size
        void unsafe_insert(key_type key, const pkdb_value_type& value)
        {
            _map[key].push_back(value);
        }

        /// \brief Replace a phylo-kmer in the database.
        /// \details This method will remove all values associated to the given key,
        /// and will add a new value.
        void replace(key_type key, const pkdb_value_type& value)
        {
            _map[key].clear();
            _map[key].push_back(value);
        }

        void sort()
        {
            for (auto& [key, scores] : _map)
            {
                std::sort(_map[key].begin(), _map[key].end(),
                          [](auto pk1, auto pk2) { return pk1.branch < pk2.branch; });
            }
        }

        /// Returns the serialation protocol version
        [[nodiscard]]
        unsigned int version() const noexcept
        {
            return _version;
        }

        void set_version(unsigned int version)
        {
            _version = version;
        }

        void set_mu(float mu)
        {
            _mu = mu;
        }

        float get_mu()
        {
            return _mu;
        }

        /// An array of pairs (kmer, filter value). Should be sorted by filter value.
        /// Determines the order in which k-mers are serialized.
        std::vector<kmer_fv> kmer_order;

    private:
        storage _map;

        /// \brief K-mer size.
        /// \details This number is given by user to the constructor. We can not guarantee
        /// that the keys stored in hash tables actually correspond to substrings of the size _kmer_size.
        /// Example: DNA ('A', 'C', 'G', 'T")
        ///     key('AAA') == key('AA') == 0
        /// e.g. putting 0 in hashtable, we assume it corresponds to 'AAA' having _kmer_size == 3,
        /// but we can not guarantee that it was not calculated for another k-mer size by mistake.
        size_t _kmer_size;

        /// \brief Score threshold paramenter.
        /// \sa core::score_threshold
        phylo_kmer::score_type _omega;

        /// \brief Sequence type: DNA or Proteins.
        /// \details Serialized within the database, needed to ensure
        /// we work with the correct sequence type.
        std::string _sequence_type;

        /// KEEP_POSITIONS is a compile-time constant, and the serialization
        /// protocol depends on it. This variable is true if during the deserialization
        /// of a database with positions we loaded them.
        bool _positions_loaded;

        /// \brief Newick formatted phylogenetic tree
        std::string _tree;

        /// The tree index (# nodes in the subtree, total subtree branch length
        /// for every node in a plain array indexed by postorder IDs
        std::vector<phylo_node::node_index> _tree_index;

        /// Serialization protocol version
        unsigned int _version;

        /// The proportion of phylo-k-mers that have been loaded
        float _mu;
    };

    using phylo_kmer_db = _phylo_kmer_db<phylo_kmer>;

    namespace impl
    {
        /// \brief A search result wrapper around a collection of pairs [branch, score]
        /// to iterate over.

        template <typename ValueType>
        class search_result
        {
        public:
            using const_iterator = typename ValueType::const_iterator;

            search_result() noexcept = default;
            search_result(const_iterator begin, const_iterator end) noexcept
                : _begin{ begin }, _end{ end }
            {}

            search_result(const search_result&) noexcept = default;
            ~search_result() noexcept = default;

            [[nodiscard]]
            const_iterator begin() const noexcept
            {
                return _begin;
            }

            [[nodiscard]]
            const_iterator end() const noexcept
            {
                return _end;
            }

            [[nodiscard]]
            size_t size() const noexcept
            {
                return std::distance(_begin, _end);
            }
        private:
            const_iterator _begin;
            const_iterator _end;
        };
    }

}

#endif
