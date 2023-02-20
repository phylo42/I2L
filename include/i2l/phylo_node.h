#ifndef XPAS_BUILD_PHYLO_NODE_H
#define XPAS_BUILD_PHYLO_NODE_H

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <type_traits>

namespace i2l
{
    class phylo_tree;
    class phylo_node;
}

namespace i2l
{
    //////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief A node of a phylogenetic tree.
    class phylo_node
    {
        friend i2l::phylo_tree;
    public:
        /// Member types

        /// \brief Node post-/pre-order id type
        using id_type = int;

        /// \brief Branch length type
        using branch_length_type = double;

        /// Indexing information for every node.
        /// This struct is auxilary, used in serialization.h
        struct node_index
        {
            size_t subtree_num_nodes;
            branch_length_type subtree_total_length;
        };

        phylo_node();
        phylo_node(std::string label, branch_length_type branch_length, phylo_node* parent);
        phylo_node(const phylo_node& other) = delete;
        phylo_node& operator=(const phylo_node&) = delete;
        ~phylo_node() noexcept;

        /// WARNING: this operator only checks for the id and label fields
        bool operator==(const phylo_node& rhs) const noexcept;
        bool operator!=(const phylo_node& rhs) const noexcept;

        [[nodiscard]]
        std::string get_label() const noexcept;

        void set_label(const std::string& label);

        [[nodiscard]]
        phylo_node* get_parent() const noexcept;

        void set_parent(phylo_node* parent);

        [[nodiscard]]
        id_type get_preorder_id() const noexcept;

        [[nodiscard]]
        id_type get_postorder_id() const noexcept;

        [[nodiscard]]
        branch_length_type get_branch_length() const noexcept;

        void set_branch_length(branch_length_type length);

        [[nodiscard]]
        branch_length_type get_subtree_branch_length() const noexcept;

        void set_subtree_branch_length(branch_length_type length);

        /// Returns the total number of nodes in the subtree
        [[nodiscard]]
        size_t get_num_nodes() const noexcept;

        void set_num_nodes(size_t num_nodes);

        /// Returns the total number of leaves in the subtree
        [[nodiscard]]
        size_t get_num_leaves() const noexcept;

        void set_num_leaves(size_t num_leaves);

        [[nodiscard]]
        const std::vector<phylo_node*>& get_children() const;

        /// Clean node and fill with the default values. Used in the default constructor
        void clean();

        void add_child(phylo_node* node);

        void remove_child(phylo_node* node);

        [[nodiscard]]
        bool is_leaf() const noexcept;

        [[nodiscard]]
        bool is_root() const noexcept;

        /// Creates a deep copy of this node. We prefer to have this method, not the
        /// copy constructor to make sure we never copy nodes by mistake
        [[nodiscard]]
        phylo_node* copy();

    private:
        id_type _preorder_id;
        id_type _postorder_id;

        std::string _label;

        /// The length of the branch to the parent
        branch_length_type _branch_length;

        /// The total branch length in the subtree
        branch_length_type _subtree_branch_length;

        /// The total number of nodes in the subtree excluding this node
        size_t _num_nodes;

        /// The total number of leaves in the subtree including this node
        size_t _num_leaves;

        std::vector<phylo_node*> _children;
        phylo_node* _parent;
    };

}
#endif //XPAS_BUILD_PHYLO_NODE_H
