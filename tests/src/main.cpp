#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()

#include <boost/filesystem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <i2l/phylo_kmer_db.h>
#include <i2l/serialization.h>
#include <i2l/kmer_iterator.h>
#include <i2l/newick.h>
#include <i2l/phylo_tree.h>
#include <i2l/fasta.h>

namespace fs = boost::filesystem;
using Catch::Approx;

auto create_test_map()
{
    std::unordered_map<i2l::phylo_kmer::key_type,
        std::unordered_map<i2l::phylo_kmer::branch_type, i2l::phylo_kmer::score_type>> values =
        {
            {
                0, { { 0, 0.00f } }
            },
            {
                1, { { 0, 0.10f }, { 1, 0.11f } }
            },
            {
                2, { { 0, 0.20f }, { 1, 0.21f }, { 2, 0.22f } }
            },
            {
                3, { { 1, 0.31f }, { 2, 0.32f } }
            },
            {
                4, { { 2, 0.42f } }
            }
        };
    return values;
}

template<typename MapType>
i2l::phylo_kmer_db create_db_from_map(const MapType& values, size_t kmer_size, i2l::phylo_kmer::score_type omega)
{
    i2l::phylo_kmer_db db { kmer_size, omega, i2l::seq_type::name, "" };
    for (const auto& [key, entries] : values)
    {
        for (const auto& [branch, score] : entries)
        {
            db.unsafe_insert(key, { branch, score });
        }
    }
    return db;
}

TEST_CASE("Database size", "[database]")
{
    {
        const auto values = create_test_map();
        const auto db = create_db_from_map(values, 3, 1.0);
        REQUIRE(db.size() == values.size());
    }

    {
        const i2l::phylo_kmer_db db { 3, 1.0, i2l::seq_type::name, "" };
        REQUIRE(db.size() == 0);
    }
}

TEST_CASE("K-mer size and omega", "[database]")
{
    const size_t kmer_size = 5;
    const i2l::phylo_kmer::score_type omega = 1.0;
    const i2l::phylo_kmer_db db { kmer_size, omega, i2l::seq_type::name, "" };

    REQUIRE(db.kmer_size() == kmer_size);
}

template<typename MapType>
void compare_db(const MapType& values, const i2l::phylo_kmer_db& db)
{
    for (const auto& [key, entries] : values)
    {
        auto db_entries = db.search(key);
        REQUIRE((bool)db_entries);
        for (const auto&[branch, score] : *db_entries)
        {
            REQUIRE(entries.find(branch) != entries.end());
            REQUIRE(entries.find(branch)->second == Approx(score));
        }
    }
}

TEST_CASE("Database search", "[database]")
{
    const auto filename = fs::unique_path().string();
    const auto values = create_test_map();
    const auto db = create_db_from_map(values, 3, 1.0 );
    compare_db(values, db);
}


TEST_CASE("(De-)serialization", "[database]")
{
    const auto filename = fs::unique_path().string();
    const auto values = create_test_map();
    const size_t kmer_size = 3;
    const i2l::phylo_kmer::score_type omega = 1.0;

    {
        const auto db = create_db_from_map(values, kmer_size, omega);
        i2l::save(db, filename);
    }

    {
        const auto db = i2l::load(filename);
        REQUIRE(db.size() == values.size());
        REQUIRE(db.kmer_size() == kmer_size);
        REQUIRE(db.omega() == omega);
    }
}

/// Iterate over all the combinations with repetition
/// http://shoaib-ahmed.com/2018/for-each-combination-with-repetetion-c++/
template<typename V, typename Callable>
void for_each_combination(V &v, size_t gp_sz, Callable f) {
    V gp(gp_sz);
    auto total_n = std::pow(v.size(), gp.size());
    for (auto i = 0; i < total_n; ++i) {
        auto n = i;
        for (auto j = 0ul; j < gp.size(); ++j) {
            gp[gp.size() - j - 1] = v[n % v.size()];
            n /= v.size();
        }
        f(gp);
    }
}

TEST_CASE("Encoding and decoding k-mers", "[kmers]")
{
    auto alphabet = std::vector<char>{ 'A', 'C', 'G', 'T', '-', 'N' };
    const size_t kmer_size = 3;
    size_t count = 0;
    for_each_combination(alphabet, kmer_size,
                         [&count](std::vector<char>& bases) {
                             const auto kmer = std::string{ bases.begin(), bases.end() };
                             if (const auto key = i2l::encode_kmer<i2l::no_ambiguity_policy>(kmer); key)
                             {
                                 REQUIRE(kmer == i2l::decode_kmer(*key, kmer.size()));
                                 REQUIRE(*key == count);
                                 ++count;
                             }
                         });
}

TEST_CASE("i2l::to_kmers iteration", "[kmers]")
{
    /// Simple iteration
    const auto long_read = std::string{ "--TTTAT-AAATGNNNN-CAAAN.NNTTTT---" };
    const size_t kmer_size = 4;

    size_t count = 0;
    for (const auto& [kmer, code] : i2l::to_kmers<i2l::no_ambiguity_policy>(long_read, kmer_size))
    {
        REQUIRE(i2l::encode_kmer<i2l::no_ambiguity_policy>(kmer) == code);
        ++count;
    }
    REQUIRE(count == 6);
}

TEST_CASE("i2l::to_kmers empty set", "[kmers]")
{
    const size_t kmer_size = 3;

    for (const auto& [kmer, code] : i2l::to_kmers<i2l::no_ambiguity_policy>("", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : i2l::to_kmers<i2l::no_ambiguity_policy>("-", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : i2l::to_kmers<i2l::no_ambiguity_policy>("---", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : i2l::to_kmers<i2l::no_ambiguity_policy>("----", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : i2l::to_kmers<i2l::no_ambiguity_policy>("NNN", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : i2l::to_kmers<i2l::no_ambiguity_policy>("AAN", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : i2l::to_kmers<i2l::no_ambiguity_policy>("ANA", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : i2l::to_kmers<i2l::no_ambiguity_policy>("ANA", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : i2l::to_kmers<i2l::no_ambiguity_policy>("NAA", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : i2l::to_kmers<i2l::no_ambiguity_policy>("-AA", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : i2l::to_kmers<i2l::no_ambiguity_policy>("A-A", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : i2l::to_kmers<i2l::no_ambiguity_policy>("AA-", kmer_size))
    {
        REQUIRE_FALSE(true);
    }
}

TEST_CASE("i2l::io::parse_newick", "[tree]")
{
    std::string newick = "((A:0.05,B:0.1):0.15,(C:0.2,D:0.25):0.3):0.35;";
    /// labels in the DFS post-order
    const std::vector<std::string> labels = { "A", "B", "", "C", "D", "", ""};
    /// branch lengths in the DFS post-order
    const std::vector<double> lengths = { 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35 };

    auto tree = i2l::io::parse_newick(newick);

    /// iterate over nodes and check if the labels and branch lengths are correct
    size_t i = 0;
    for (auto& node : tree)
    {
        REQUIRE(node.get_label() == labels[i]);
        REQUIRE(Approx(node.get_branch_length()) == lengths[i]);

        // ensure we can not create trees from non-root nodes
        if (node.get_parent())
        {
            REQUIRE_THROWS(i2l::phylo_tree{ &node });
        }

        ++i;
    }

    REQUIRE(tree.get_node_count() == 7);
}

TEST_CASE("i2l::io::parse_newick complex labels", "[tree]")
{
    std::string newick = "('A':0.1,'B is complex':0.2,('C, (C is complex, really)':0.3,'D does not make sense,,,((((':0.4)E:0.5)F;";
    /// labels in the DFS post-order
    const std::vector<std::string> labels = { "A", "B is complex",
                                              "C, (C is complex, really)",
                                              "D does not make sense,,,((((", "E", "F"};
    /// branch lengths in the DFS post-order
    const std::vector<double> lengths = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.0 };

    auto tree = i2l::io::parse_newick(newick);

    /// iterate over nodes and check if the labels and branch lengths are correct
    size_t i = 0;
    for (auto& node : tree)
    {
        REQUIRE(node.get_label() == labels[i]);

        REQUIRE(Approx(node.get_branch_length()) == lengths[i]);

        ++i;
    }

    REQUIRE(tree.get_node_count() == 6);
}


TEST_CASE("i2l::impl::next_by_postorder", "[tree]")
{
    std::string newick = "((A:0.05,B:0.1):0.15,(C:0.2,D:0.25):0.3):0.35;";
    const auto tree = i2l::io::parse_newick(newick);

    for (const auto& node : tree)
    {
        auto next = i2l::impl::next_by_postorder(&node);

        /// next can be null only for the root node
        if (!next)
        {
            REQUIRE(node.get_parent() == nullptr);
        }
        else
        {
            REQUIRE(node.get_postorder_id() + 1 == next->get_postorder_id());
        }
    }
}

TEST_CASE("i2l::phylo_tree::get_by_postorder_id", "[tree]")
{
    std::string newick = "((A:0.05,B:0.1):0.15,(C:0.2,D:0.25):0.3):0.35;";

    const auto tree = i2l::io::parse_newick(newick);

    int iteration_count = 0;
    /// test phylo_tree::get_post_order
    for (const auto& node : tree)
    {
        // make sure this node has got the right post-order id
        REQUIRE(node.get_postorder_id() == iteration_count);

        /// find the same node in the tree by their post-order id
        const auto postorder_id = node.get_postorder_id();
        const auto found = tree.get_by_postorder_id(postorder_id);

        // make sure we found it
        REQUIRE(found);

        const auto& node_found = *found;
        // make sure it's not a null pointer
        REQUIRE(node_found);

        // compare fields
        REQUIRE(node_found->get_label() == node.get_label());
        REQUIRE(node_found->get_postorder_id() == node.get_postorder_id());
        REQUIRE(node_found->get_preorder_id() == node.get_preorder_id());
        REQUIRE(node_found->get_children() == node.get_children());
        REQUIRE(node_found->get_branch_length() == node.get_branch_length());

        ++iteration_count;
    }
}

TEST_CASE("i2l::visit_tree", "[tree]")
{
    std::string newick = "((A:0.05,B:0.1):0.15,(C:0.2,D:0.25):0.3):0.35;";
    auto tree = i2l::io::parse_newick(newick);

    std::vector<double> total_lengths = { 0.05, 0.1, 0.3, 0.2, 0.25, 0.75, 1.4 };
    size_t i = 0;

    /// can't start visiting nullptr
    REQUIRE_THROWS(i2l::visit_subtree(nullptr));

    /// Here we also test non-const iteration
    for (auto& node : tree)
    {
        /// Test if we can start visiting the subtree
        using iterator = i2l::postorder_tree_iterator<false>;
        REQUIRE_NOTHROW(i2l::visit_subtree<iterator>(&node));

        /// run DFS from a node, calculating the total subtree branch length
        double total_length = 0.0;
        for (const auto& subtree_node : i2l::visit_subtree(&node))
        {
            total_length += subtree_node.get_branch_length();
        }

        REQUIRE(Approx(total_length) == total_lengths[i]);

        ++i;
    }
}

TEST_CASE("i2l::io::read_fasta", "[utils]")
{
    const auto seq_records = std::vector<i2l::seq_record>{
        {"1", "AAAAA"},
        {"2", "AAAAC"},
        {"3", "AAAAG"},
        {"4", "AAAAT"},
        {"5", "AAACA"},
        {"6", "AAACC"},
        {"7", "AAACG"},
        {"8", "AAACT"},
        {"9", "AAAGA"},
        {"10", "AAAGC"},
        {"11", "AAAGG"},
        {"12", "AAAGT"},
    };

    /// write sequences to a temporary file
    const auto filename = fs::unique_path().string();
    std::ofstream out(filename);
    for (const auto& seq : seq_records)
    {
        out << ">" << seq.header() << std::endl << seq.sequence() << std::endl;
    }

    /// read with default batch size
    size_t i = 0;
    for (const auto& seq : i2l::io::read_fasta(filename))
    {
        REQUIRE(seq == seq_records[i]);
        ++i;
    }

    /// batch_size < seq_records.size()
    /// seq_records.size() % batch_size == 0
    size_t batch_size = 4;
    i = 0;
    for (const auto& seq : i2l::io::read_fasta(filename, batch_size))
    {
        REQUIRE(seq == seq_records[i]);
        ++i;
    }

    /// batch_size < seq_records.size()
    /// seq_records.size() % batch_size != 0
    batch_size = 5;
    i = 0;
    for (const auto& seq : i2l::io::read_fasta(filename, batch_size))
    {
        REQUIRE(seq == seq_records[i]);
        ++i;
    }

    /// batch_size == seq_records.size()
    batch_size = seq_records.size();
    i = 0;
    for (const auto& seq : i2l::io::read_fasta(filename, batch_size))
    {
        REQUIRE(seq == seq_records[i]);
        ++i;
    }

    /// batch_size > seq_records.size()
    batch_size = 2 * seq_records.size();
    i = 0;
    for (const auto& seq : i2l::io::read_fasta(filename, batch_size))
    {
        REQUIRE(seq == seq_records[i]);
        ++i;
    }
}