#include <fstream>
#include <sstream>
#include <oneapi/tbb/parallel_for.h>
#include <maki/tree.hpp>

namespace maki_tree {

    Node::Node(size_t id_) : id(id_), edge(), dist(), metadata() {}

    bool Node::hasEdgeLabel() const noexcept { return edge != id; }

    bool NodeNameMap::PairLessThan::operator()(const Pair &lhs, const Pair &rhs) const
    { return lhs.first < rhs.first; }

    bool NodeNameMap::PairLessThan::operator()(const Pair &ref, const std::string &q) const
    { return ref.first < q; }

    void NodeNameMap::sortIds()
    { std::sort(_ids.begin(), _ids.end(), PairLessThan{}); }

    void NodeNameMap::mapNames() {
        auto mx = std::max_element(_ids.cbegin(), _ids.cend(), [](const Pair &lhs, const Pair &rhs)->bool
            { return lhs.second < rhs.second; });
        _names.resize(mx->second + 1, 0);

        for (size_t i = 0; i < _ids.size(); ++i)
            _names[_ids[i].second] = i;
    }

    void NodeNameMap::finaliseIds() {
        sortIds();
        mapNames();
    }

    void NodeNameMap::insert(std::string &&name, size_t id)
    { _ids.emplace_back(std::move(name), id); }

    NodeNameMap::ConstIterator NodeNameMap::name(size_t id) const
    { return _ids.cbegin() + _names[id]; }

    NodeNameMap::ConstIterator NodeNameMap::id(const std::string &name) const {
        auto it = std::lower_bound(_ids.cbegin(), _ids.cend(), name, PairLessThan{});
        if (it->first == name)
            return it;
        else return end();
    }

    NodeNameMap::ConstIterator NodeNameMap::begin() const
    { return _ids.begin(); }

    NodeNameMap::ConstIterator NodeNameMap::end() const
    { return _ids.cend(); }

    const Node& Tree::operator[](size_t id) const
    { return pool[id-1]; }

    size_t Tree::countNodes() const
    { return pool.size(); }

    size_t Tree::countEdges() const
    { return countNodes() - 1; }

    size_t Tree::parent(size_t id) const
    { return bps.rank(bps.enclose(bps.select(id))); }

    bool Tree::hasLabel(size_t id) const {
        auto it = name(id);
        return it != end() && !it->first.empty();
    }

    NodeNameMap::ConstIterator Tree::name(size_t id) const
    { return names.name(id); }

    NodeNameMap::ConstIterator Tree::id(const std::string &name) const
    { return names.id(name); }

    NodeNameMap::ConstIterator Tree::end() const
    { return names.end(); }

    void Tree::printNames() const {
        for (const auto &[k, v] : names) {
            std::cout << k << '\n';
        }
        std::cout << std::endl;
    }

    size_t Tree::childId(size_t id, size_t n) const {
        size_t pos = bps.select(id) + 1U;
        for (size_t c = 0; c < n; ++c)
            pos = bps.find_close(pos) + 1u;
        return bps.rank(pos);
    }

    size_t Tree::countChildren(size_t id) const {
        size_t pos = bps.select(id) + 1,
               n = 0;
        while (bpr[pos] == 1) {
            pos = bps.find_close(pos) + 1;
            ++n;
        }
        return n;
    }

    size_t Tree::countSubnodes(size_t id) const {
        size_t pos = bps.select(id);
        return bps.rank(bps.find_close(pos)) - bps.rank(pos);
    }

    size_t Tree::countSubleaves(size_t id) const {
        size_t lpos = bps.select(id);
        return lrs.rank(bps.find_close(lpos)) - lrs.rank(lpos);
    }

    size_t newick::detail::parseNodeData(const std::string &newick, NodeNameMap &names, size_t pos, Node *node) {
        while (pos < newick.size() && std::isspace(newick[pos]))
            ++pos;

        // parse label
        if (pos < newick.size() && newick[pos] == '\'') {
            // quoted label
            size_t labelEnd = newick.find_first_of('\'', pos + 1);
            if (labelEnd == std::string::npos) {
                throw NewickParsingError("Unterminated quoted label at position " + std::to_string(pos));
            } else {
                names.insert(newick.substr(pos + 1, labelEnd - pos - 1), node->id);
                pos = labelEnd + 1;
            }
        } else {
            // unquoted label
            size_t labelEnd = newick.find_first_of(":,);{[", pos);
            names.insert(newick.substr(pos, labelEnd - pos), node->id);
            pos = labelEnd;
        }

        // parse edge id if present
        if (pos < newick.size() && newick[pos] == '{') {
            size_t labelEnd = newick.find_first_of('}', pos + 1);
            if (labelEnd == std::string::npos) {
                throw NewickParsingError("Unterminated curly brace at position " + std::to_string(pos));
            } else {
                size_t cvt;
                try {
                    node->edge = std::stoul(newick.substr(pos + 1, labelEnd - pos - 1), &cvt);
                } catch(const std::exception& e) {
                    goto noparse;
                }
                
                if (cvt != (labelEnd - pos - 1))
                    goto noparse;
                else goto parse;

                noparse: throw NewickParsingError("Could not parse edge label (must be integer edge labels) at position " + std::to_string(pos + 1));
                parse: pos = labelEnd + 1;
            }
        } else {
            node->edge = node->id;
        }

        // parse edge length if present
        if (pos < newick.size() && newick[pos] == ':') {
            size_t distBegin = ++pos;
            while (pos < newick.size() &&
                   (std::isdigit(newick[pos]) || newick[pos] == '.' ||
                    newick[pos] == 'e' || newick[pos] == '-' || newick[pos] == '+'))
            {
                ++pos;
            }

            try {
                node->dist = std::stod(newick.substr(distBegin, pos - distBegin));
            } catch(const std::exception& e) {
                throw NewickParsingError("Invalid edge length at position " + std::to_string(distBegin));
            }
        }
        
        // parse additional metadata if present
        if (pos < newick.size() && newick[pos] == '[') {
            while (newick[pos] != ']') {
                size_t kvpEnd = newick.find_first_of(",]", pos + 1);
                if (kvpEnd == std::string::npos) {
                    throw NewickParsingError("Unterminated metadata section at position " + std::to_string(pos + 1));
                } else {
                    size_t eqPos = newick.find_first_of('=', pos + 1);
                    if (eqPos == std::string::npos || eqPos > kvpEnd) {
                        throw NewickParsingError("Misformatted metadata key-value pair at position " + std::to_string(pos + 1));
                    } else {
                        size_t keyStart = pos + 1;
                        if (keyStart < newick.size() && newick[keyStart] == '&')
                            ++keyStart;
                        node->metadata[newick.substr(keyStart, eqPos - keyStart)] = newick.substr(eqPos + 1, kvpEnd - eqPos - 1);
                    }
                }
                pos = kvpEnd;
            }
            ++pos;
        }

        if (pos < newick.size() && newick[pos] == ',')
            ++pos;

        return pos;
    }

    void Tree::load(const char *fpath) {
        std::ifstream     ifs(fpath);
        std::stringstream buffer{};
        buffer << ifs.rdbuf();

        std::string data = buffer.str();
        data.erase(std::remove(data.begin(), data.end(), '\n'), data.end());

        parse(data);
    }

    Tree::Tree(Tree &&rval)
        : pool(std::move(rval.pool))
        , names(std::move(rval.names))
        , bpr(std::move(rval.bpr))
        , bps(std::move(rval.bps))
        , lrs(std::move(rval.lrs))
    {
        bps.set_vector(&bpr);
        lrs.set_vector(&bpr);
    }

    Tree::Tree(const char *fpath) { load(fpath); }

    void Tree::parse(const std::string &data) {
        std::stack<size_t> lineage{};

        size_t i = 0,
               L = data.size(),
               finished = L + 10,
               id = Tree::RootId;
        while (i < L) {
            const char &x = data[i];
            switch (x)
            {
            case '(':
                pool.emplace_back(id++);
                lineage.push(pool.size() - 1);
                bpr.push_back(1);
                ++i;
                break;
            case ')':
                i = newick::detail::parseNodeData(data, names, i+1, &pool[lineage.top()]);
                lineage.pop();
                bpr.push_back(0);
                break;
            case ';':
                i = finished;
                break;
            default:
                pool.emplace_back(id++);
                i = newick::detail::parseNodeData(data, names, i, &pool.back());
                bpr.push_back(1);
                bpr.push_back(0);
                break;
            }
        }

        if (i != finished)
            throw newick::NewickParsingError("Data ended before ';'");
        else if (!lineage.empty())
            throw newick::NewickParsingError("Uneven number of parentheses in Newick string");

        sdsl::util::init_support(bps, &bpr);
        sdsl::util::init_support(lrs, &bpr);
        names.finaliseIds();
    }

    Tree::Tree(const std::string &data) { parse(data); }

    Tree::cursor::cursor(const Tree &tree, size_t nid_)
        : nid(nid_)
        , cx(0)
        , length(tree.countChildren(nid))
    {}

    bool Tree::cursor::hasChildren() const { return cx < length; }

    void Tree::cluster::removeChild(size_t cx) {
        children.erase(std::remove(children.begin(), children.end(), cx),
                       children.end());
    }

    void Tree::cluster::merge(Tree::cluster &rhs) {
        nodes.reserve(size() + rhs.size());
        nodes.insert(nodes.end(), rhs.nodes.begin(), rhs.nodes.end());
        children = std::move(rhs.children);
    }

    void Node::toNewick(std::ostream &os, const std::string &name) const {
        if (name.contains(' ')) os << '\'' << name << '\'';
        else                    os << name;
        // write edge label if present
        if (edge != id)
            os << '{' << edge << '}';
        // write edge distance
        if (dist > 0.0)
            os << ':' << dist;
        // write edge metadata
        if (!metadata.empty()) {
            os << "[&";
            auto it = metadata.cbegin(), mend = metadata.cend();
            os << it->first << '=' << it->second;
            while (++it != mend)
                os << ',' << it->first << '=' << it->second;
            os << ']';
        }
    }

    void Tree::toNewick(std::ostream &os) const {
        std::stack<cursor> lineage;
        
        // root case
        lineage.emplace(*this, root());
        os << '(';

        // recurse tree structure
        while (!lineage.empty()) {
            while (!lineage.empty() && !lineage.top().hasChildren()) {
                const Node &n = operator[](lineage.top().nid);
                os << ')';
                // write node label
                n.toNewick(os, name(n.id)->first);
                lineage.pop();
            }

            if (!lineage.empty()) {
                cursor &c = lineage.top();
                if (c.cx > 0) os << ',';
                size_t cid = childId(c.nid, c.cx++);
                if (countChildren(cid)) {
                    os << '(';
                    lineage.emplace(*this, cid);
                } else {
                    const Node &n = operator[](cid);
                    n.toNewick(os, name(cid)->first);
                }
            }
        }

        // finish
        os << ';';
    }

    Tree Tree::reduce(const std::unordered_set<std::string> &nodes) const {
        assert(!nodes.empty());
        auto targets = retrieveIds(nodes);
        auto clusters = currentStructure();
        reduceStructures(clusters, targets);
        return Tree(*this, clusters);
    }

    std::unordered_set<size_t> Tree::retrieveIds(const std::unordered_set<std::string> &nodes) const {
        std::unordered_set<size_t> results;
        results.reserve(nodes.size());
        std::vector<std::string> missingNames;
        for (const std::string &qry : nodes) {
            auto it = id(qry);
            if (it == end()) missingNames.push_back(qry);
            else             results.insert(it->second);
        }
        if (!missingNames.empty()) {
            std::string missingNameString = std::move(missingNames[0]);
            for (size_t i = 1; i < std::min(missingNames.size(), (size_t)10); ++i) {
                missingNameString.push_back(',');
                missingNameString += missingNames[i];
            }
            if (missingNameString.size() > 10) missingNameString += ",...";
            throw std::runtime_error("Unrecognised genome names (" + std::to_string(missingNames.size()) + "): " + missingNameString);
        }
        return results;
    }

    std::vector<Tree::cluster> Tree::currentStructure() const {
        std::vector<cluster> clusters(countNodes());
        for (size_t i = 0; i < clusters.size(); ++i) {
            clusters[i].nodes.push_back(i+1);
            size_t nchild = countChildren(i+1);
            for (size_t c = 0; c < nchild; ++c) {
                size_t j = childId(i+1, c);
                clusters[j-1].parent = i;
                clusters[i].children.push_back(j-1);
            }
        }
        return clusters;
    }

    void Tree::reduceStructures(std::vector<cluster> &clusters, const std::unordered_set<size_t> &targets) const {
        for (int i = clusters.size() - 1; i > -1; --i) {
            cluster &x = clusters[i];
            // can we remove this node?
            if (!targets.contains(x.head())/*not requested*/) {
                if (x.weight() == 0 /*no children*/) {
                    clusters[x.parent].removeChild(i);
                }
                // can we merge this edge?
                else if (x.weight() == 1) {
                    cluster &y = clusters[x.children[0]];
                    if (y.size() > 1 || !hasLabel(x.head()) || !targets.contains(y.head()))
                        x.merge(y);
                }
            }
        }
    }

    bool Tree::pushCluster(const Tree &src, const cluster &cstr, Node &newNode) {       
        // set edge label
        if (cstr.singular() && src[cstr.head()].hasEdgeLabel())
            newNode.edge = src[cstr.head()].edge;
        else
            newNode.edge = newNode.id;

        // collect name and node data
        std::string name;
        for (const auto &nid : cstr.nodes) {
            const Node &q = src[nid];
            // update name
            auto nameIter = src.name(nid);
            if (!nameIter->first.empty()) {
                if (name.empty())   name = nameIter->first;
                else                name += "; " + nameIter->first;
            }
            // update dist
            newNode.dist += q.dist;
            // update metadata
            for (const auto &[k_, v_] : q.metadata) newNode.metadata[k_] = v_;
        }
        names.insert(std::move(name), newNode.id);

        // update bpr structure
        if (cstr.weight()) {
            bpr.push_back(1);/*branch*/
            return true;
        } else {
            bpr.push_back(1); bpr.push_back(0);/*leaf*/
            return false;
        }
    }

    Tree::Tree(const Tree &t, const std::vector<cluster> &clusters) {
        assert(!clusters.empty());
        std::stack<std::pair<size_t,size_t>> walk;
        size_t id = Tree::RootId;

        // base case - root
        if (pushCluster(t, clusters[0], pool.emplace_back(id++)))
            walk.emplace(0, 0);

        // recurse
        while (!walk.empty()) {
            auto &[clusterId, childId] = walk.top();
            if (childId == clusters[clusterId].weight()) {
                walk.pop();
                bpr.push_back(0);
                continue;
            } else {
                size_t nextChildId = clusters[clusterId].children[childId++];
                if (pushCluster(t, clusters[nextChildId], pool.emplace_back(id++)))
                    walk.emplace(nextChildId, 0);/*recurse*/
            }
        }

        // initialise rank-select support
        sdsl::util::init_support(bps, &bpr);
        sdsl::util::init_support(lrs, &bpr);
        names.finaliseIds();
    }

    namespace iter
    {
        
        NodeRange::NodeRange(size_t nid_, size_t cx_, size_t length_, size_t depth_)
            : nid(nid_)
            , cx(cx_)
            , length(length_)
            , depth(depth_)
        {}

        bool NodeRange::hasChildren() const { return cx < length; }

        PreOrderIterator::PreOrderIterator(const maki_tree::Tree &t_, size_t nid_, size_t cx_)
            : tree(t_)
            , stack()
        {
            stack.emplace(nid_, cx_, tree.countChildren(nid_), 0);
        }

        Hptr PreOrderIterator::operator*() const {
            const NodeRange &r = stack.top();
            return Hptr{ tree[r.nid], r.depth };
        }

        bool PreOrderIterator::deadEnd() const {
            return !stack.empty() && !stack.top().hasChildren();
        }

        void PreOrderIterator::nextNew() {
            while (deadEnd())
                stack.pop();

            if (!stack.empty()) {
                NodeRange &r = stack.top();
                size_t cid = tree.childId(r.nid, r.cx++);
                stack.emplace(cid, 0, tree.countChildren(cid), r.depth + 1);
            }
        }

        void PreOrderIterator::operator++() {
            if (stack.top().hasChildren()) {
                NodeRange &r = stack.top();
                size_t cid = tree.childId(r.nid, r.cx++);
                stack.emplace(cid, 0, tree.countChildren(cid), r.depth + 1);
            } else stack.pop();
        }

        bool PreOrderIterator::operator!=([[maybe_unused]] const PreOrderTourSentinel &) const {
            return !stack.empty();
        }

        PreOrderTour::PreOrderTour(const maki_tree::Tree &tree, size_t subtree)
            : tree(tree)
            , root(subtree)
        {}

        PreOrderTour::PreOrderTour(const Tree& tree)
            : PreOrderTour(tree, tree.root())
        {}

        PreOrderIterator PreOrderTour::begin() const
        { return PreOrderIterator(tree, root); }

        PreOrderTourSentinel PreOrderTour::end() const
        { return PreOrderTourSentinel{}; }

    } // namespace iter

    namespace lca
    {

        LcaTable::LcaTable(const Tree &t_)
            : N(t_.countNodes())
            , TL(2 * t_.countEdges())
            , ranks(N + 1, (size_t)-1)
            , tour(TL)
            , rmq()
        {
            std::vector<size_t> heights(TL);

            // perform euler tour of tree
            iter::PreOrderTour Tour(t_);
            for (auto [i, it] = std::make_pair((size_t)0u, Tour.begin()); it != Tour.end(); ++i, ++it) {
                auto p = *it;
                if (ranks[p.n.id] == (size_t)-1)
                    ranks[p.n.id] = i;
                tour[i] = p.n.id;
                heights[i] = p.height;
            }

            // initialise RMQ support
            rmq = std::make_unique<sdsl::rmq_succinct_sct<>>(&heights);
        }

        size_t LcaTable::lca(size_t ia, size_t ib) const {
            size_t ra = ranks[ia],
                   rb = ranks[ib];
            if (ra > rb)
                std::swap(ra, rb);
            return tour[(*rmq)(ra, rb)];
        }

        size_t LcaTable::lca(const std::vector<size_t> &ids) const {
            if (ids.empty())
                return 0;
            else {
                size_t min, max, r;
                min = max = ranks[ids[0]];
                for (size_t i = 1; i < ids.size(); ++i) {
                    r = ranks[ids[i]];
                    if (r < min)
                        min = r;
                    else if (r > max)
                        max = r;
                }
                return tour[(*rmq)(min, max)];
            }
        }

    } // namespace lca

} // namespace maki_tree

