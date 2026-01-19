#pragma once

#include <algorithm>
#include <cstddef>
#include <memory>
#include <stack>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#include <sdsl/vectors.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/bp_support.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <maki/utils.hpp>

namespace maki_tree {

    /* Node
    */
    struct Node {
        size_t id,      // id is pre-order depth-first rank (note rank starts at 1)
               edge;    // edge to its parent
        double dist;    // length of edge
        std::unordered_map<std::string,std::string> metadata;

        Node(size_t id_);

        /* Write data to Newick format
        */
        void toNewick(std::ostream &os, const std::string &name) const;

        /* Check/guess if a custom edge label was assigned.
        */
        bool hasEdgeLabel() const noexcept;
    };

    class NodeNameMap {
    public:
        using Pair = std::pair<std::string,size_t>;
        using ConstIterator = std::vector<Pair>::const_iterator;

    protected:
        std::vector<Pair>   _ids;
        std::vector<size_t> _names;

        friend class Tree;

        struct PairLessThan {
            bool operator()(const Pair &lhs, const Pair &rhs) const;
            bool operator()(const Pair &ref, const std::string &q) const;
        };

    protected:
        void sortIds();
        void mapNames();
        void finaliseIds();

    public:
        // constructors

        NodeNameMap() =default;

        NodeNameMap(NodeNameMap&&) =default;
        NodeNameMap& operator=(NodeNameMap&&) =default;

        NodeNameMap(const NodeNameMap &) =delete;
        NodeNameMap& operator=(const NodeNameMap &) =delete;

        // insertion

        void insert(std::string &&name, size_t id);

        // lookup

        ConstIterator name(size_t id) const;
        ConstIterator id(const std::string &name) const;

        // iteration

        ConstIterator begin() const;
        ConstIterator end() const;
    };

    namespace newick
    {

        struct NewickParsingError : public std::runtime_error {
            explicit NewickParsingError(std::string &&message) : std::runtime_error(std::move(message)) {}
        };
        
        namespace detail
        {
            
            /* Auxilliary function. Parse segment of newick string
            corresponding to branch name and metadata and assign to
            `Node` in pool.
            */
            size_t parseNodeData(const std::string &newick, NodeNameMap &names, size_t pos, Node *node);

        } // namespace detail

    } // namespace newick

    /* Tree class
    * 
    * Balanced parentheses representation of tree.
    * 
    * Tree nodes are assigned internal IDs accroding to a pre-order
    * depth-first enumeration of each node, starting from 0.
    */
    class Tree {
        std::vector<Node>           pool;
        NodeNameMap                 names;
        sdsl::bit_vector            bpr;     // balanced parenthesis repr.
        sdsl::bp_support_sada<>     bps;     // rs-support for bpr
        sdsl::rank_support_v5<10,2> lrs;     // leaf rank support

    public:
        static constexpr size_t RootId = 1u;

        // constructors

        Tree() =default;
        Tree(Tree&&);
        Tree(const char *);
        Tree(const std::string&);

        /* Load a newick-formatted file into a `Tree` object.
        */
        void load(const char*);

        /* Parse a newick-formatted string into a `Tree` object.
        */
        void parse(const std::string &data);

        /* Write structure to Newick
        */
        void toNewick(std::ostream &os) const;

        // access

        const Node& operator[](size_t id) const;
        static constexpr size_t root() { return Tree::RootId; }
        size_t countNodes() const;
        size_t countEdges() const;
        size_t parent(size_t id) const;
        const void* bprAddress() const { return static_cast<const sdsl::bit_vector*>(&bpr); }

        // name lookup

        bool hasLabel(size_t id) const;
        NodeNameMap::ConstIterator name(size_t id) const;
        NodeNameMap::ConstIterator id(const std::string &name) const;
        NodeNameMap::ConstIterator end() const;

        void printNames() const;

        // modifiers

        Tree reduce(const std::unordered_set<std::string> &nodes) const;

        // structural and traversal queries

        /* Find pool position of (n+1)-order child (indexing starts at 0).
        Does not perform checking of validity of request.
        */
        size_t childId(size_t id, size_t n) const;

        /* Count the number of immediate children of a given node.
        */
        size_t countChildren(size_t id) const;

        /* Count the number of nodes below a given node.
        */
        size_t countSubnodes(size_t id) const;

        /* Count the number of leaves below a given node.
        */
       size_t countSubleaves(size_t id) const;

    protected:
        struct cursor {
            size_t nid, cx, length;
            cursor(const Tree &tree, size_t nid_);
            bool hasChildren() const;
        };

        /** cluster
         * 
         * Represent a cluster of nodes.
         */
        struct cluster  {
            std::vector<size_t> nodes;      // values refer to node IDs
            std::vector<size_t> children;   // values refer to other clusters
            size_t              parent;     // values refers to other clusters

            constexpr size_t size() const noexcept { return nodes.size(); }
            constexpr size_t weight() const noexcept { return children.size(); }
            constexpr size_t empty() const noexcept { return nodes.empty(); }
            constexpr size_t head() const noexcept { return nodes[0]; }
            constexpr bool singular() const noexcept { return size() == 1; }
            
            void removeChild(size_t cx);
            void merge(cluster &rhs);
        };

        /* Export current structure as clusters
        */
        std::vector<cluster> currentStructure() const;

        /* Collect set of node IDs corresponding to the query node names,
        throwing if no match is found for any query.
        */
        std::unordered_set<size_t> retrieveIds(const std::unordered_set<std::string> &nodes) const;

        /* Collapse clusters using targets.
        */
        void reduceStructures(std::vector<cluster> &clusters, const std::unordered_set<size_t> &targets) const;

        /* Construct tree from node clustering.
        */
        Tree(const Tree &t, const std::vector<cluster> &clusters);
        
        /* Push cluster as new node. Return `true` if it the new node has children
        (i.e. it is a branch).
        */
        bool pushCluster(const Tree &src, const cluster &cstr, Node &newNode);
    };

    namespace iter
    {
        
        struct Hptr {
            const Node &n;
            size_t height;
        };

        struct NodeRange {
            size_t nid,     // pointer to node
                   cx,      // pointer to child of node
                   length,  // number of children of this node
                   depth;   // node path depth

            NodeRange(size_t nid_, size_t cx_, size_t length_, size_t depth_);

            bool hasChildren() const;
        };

        struct PreOrderTourSentinel{};

        class PreOrderIterator {
            const Tree            &tree;
            std::stack<NodeRange> stack;

            bool deadEnd() const;

        public:
            PreOrderIterator(const Tree &t_, size_t nid_, size_t cx_ = 0);

            Hptr operator*() const;
            void operator++();
            void nextNew();
            bool operator!=(const PreOrderTourSentinel &rhs) const;
        };

        /* pre-order tour of tree nodes

        Perform pre-order, depth-first traversal of tree nodes. Returns
        an iterator with `const` references to nodes of `tree`.
        */
        class PreOrderTour {
            const Tree& tree;
            size_t      root;

        public:            
            PreOrderTour(const Tree& tree, size_t subtree);
            PreOrderTour(const Tree& tree);

            PreOrderIterator begin() const;
            PreOrderTourSentinel end() const;
        };

    } // namespace iter
    
    namespace lca
    {

        /* Lowest Common Ancestor Table
        */
        class LcaTable {
            size_t N, TL;
            std::vector<size_t> ranks, tour;
            std::unique_ptr<sdsl::rmq_succinct_sct<>> rmq;
        public:
            LcaTable(const Tree &t_);
            size_t lca(size_t a, size_t b) const;
            size_t lca(const std::vector<size_t> &ids) const;
        };

    } // namespace lca

} // namespace maki_tree
