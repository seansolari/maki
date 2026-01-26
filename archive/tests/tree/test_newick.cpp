#include <string>
#include <unordered_set>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <maki/tree.hpp>

class NewickParserTest : public ::testing::Test
{
protected:
    maki_tree::Tree autoReduce(const std::string &newick, const std::unordered_set<std::string> &nodes) {
        maki_tree::Tree t;
        t.parse(newick);
        return t.reduce(nodes);
    }
};

// 1. Basic tree without labels or edge lengths
TEST_F(NewickParserTest, ParsesSimpleTree)
{
    std::string newick = "(A,B,C);";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_EQ(tree.countNodes(), 4);
    EXPECT_EQ(tree.countEdges(), 3);
    EXPECT_EQ(tree.countChildren(tree.root()), 3);

    // ids
    EXPECT_EQ(tree.root(), maki_tree::Tree::RootId);
    EXPECT_EQ(tree.childId(tree.root(), 0), 2);
    EXPECT_EQ(tree.childId(tree.root(), 1), 3);
    EXPECT_EQ(tree.childId(tree.root(), 2), 4);

    // traversing
    for (size_t i = 0; i < tree.countChildren(tree.root()); ++i)
    {
        EXPECT_EQ(tree.parent(tree.childId(tree.root(), i)), tree.root());
    }
}

// 2. Tree with node labels
TEST_F(NewickParserTest, ParsesTreeWithNodeLabels)
{
    std::string newick = "(A,B,C)Root;";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_EQ(tree.name(tree.root())->first, "Root");
    EXPECT_EQ(tree.name(tree.childId(tree.root(), 0))->first, "A");
    EXPECT_EQ(tree.name(tree.childId(tree.root(), 1))->first, "B");
    EXPECT_EQ(tree.name(tree.childId(tree.root(), 2))->first, "C");
}

// 3. Tree with edge lengths
TEST_F(NewickParserTest, ParsesTreeWithEdgeLengths)
{
    std::string newick = "(A:0.1,B:0.2,C:0.3);";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_NEAR(tree[2].dist, 0.1, 1e-6);
    EXPECT_NEAR(tree[3].dist, 0.2, 1e-6);
    EXPECT_NEAR(tree[4].dist, 0.3, 1e-6);
}

// 4. Tree with node labels and edge lengths
TEST_F(NewickParserTest, ParsesTreeWithLabelsAndLengths)
{
    std::string newick = "(A:0.1,B:0.2,C:0.3)Root:1.0;";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_EQ(tree.name(tree.root())->first, "Root");
    EXPECT_NEAR(tree[tree.root()].dist, 1.0, 1e-6);
}

// 5. Nested tree structure
TEST_F(NewickParserTest, ParsesNestedTree)
{
    std::string newick = "((A,B),(C,D));";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_EQ(tree.countChildren(tree.root()), 2);

    // ids
    EXPECT_EQ(tree.childId(tree.root(), 0), 2);
    EXPECT_EQ(tree.childId(tree.root(), 1), 5);

    // subgroups
    EXPECT_EQ(tree.countChildren(tree.childId(tree.root(), 0)), 2);
    EXPECT_EQ(tree.countChildren(tree.childId(tree.root(), 1)), 2);
}

// 6. Multifurcation
TEST_F(NewickParserTest, ParsesMultifurcation)
{
    std::string newick = "(A,B,C,D,E);";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_EQ(tree.countChildren(tree.root()), 5);
}

// 7. Malformed input: missing semicolon
TEST_F(NewickParserTest, FailsOnMissingSemicolon)
{
    std::string newick = "(A,B,C)";
    maki_tree::Tree tree;
    EXPECT_THROW(tree.parse(newick), maki_tree::newick::NewickParsingError);
}

// 8. Malformed input: unmatched parentheses
TEST_F(NewickParserTest, FailsOnUnmatchedParentheses)
{
    std::string newick = "((A,B),C;";
    maki_tree::Tree tree;
    EXPECT_THROW(tree.parse(newick), maki_tree::newick::NewickParsingError);
}

// 9. Empty input
TEST_F(NewickParserTest, FailsOnEmptyInput)
{
    std::string newick = "";
    maki_tree::Tree tree;
    EXPECT_THROW(tree.parse(newick), maki_tree::newick::NewickParsingError);
}

// 10. Tree with quoted labels
TEST_F(NewickParserTest, ParsesQuotedLabels)
{
    std::string newick = "('A label with spaces','B_label'):0.5;";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_EQ(tree.name(tree.childId(tree.root(), 0))->first, "A label with spaces");
    EXPECT_EQ(tree.name(tree.childId(tree.root(), 1))->first, "B_label");
    EXPECT_NEAR(tree[tree.root()].dist, 0.5, 1e-6);
}

// 11. Tree with integer edge labels
TEST_F(NewickParserTest, ParsesTreeWithIntegerEdgeLabels)
{
    std::string newick = "(A{1},B{2},C{3});";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_EQ(tree.name(tree.childId(tree.root(), 0))->first, "A");
    EXPECT_EQ(tree[tree.childId(tree.root(), 0)].edge, 1);
    EXPECT_EQ(tree[tree.childId(tree.root(), 1)].edge, 2);
    EXPECT_EQ(tree[tree.childId(tree.root(), 2)].edge, 3);
}

// 12. Tree with integer edge labels and edge lengths
TEST_F(NewickParserTest, ParsesTreeWithIntegerEdgeLabelsAndLengths)
{
    std::string newick = "(A{10}:0.1,B{20}:0.2,C{30}:0.3);";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_EQ(tree[tree.childId(tree.root(), 0)].edge, 10);
    EXPECT_NEAR(tree[tree.childId(tree.root(), 0)].dist, 0.1, 1e-6);
    EXPECT_EQ(tree[tree.childId(tree.root(), 1)].edge, 20);
    EXPECT_NEAR(tree[tree.childId(tree.root(), 1)].dist, 0.2, 1e-6);
}

// 13. Tree with root node having an integer edge label
TEST_F(NewickParserTest, ParsesRootWithIntegerEdgeLabel)
{
    std::string newick = "(A{1},B{2})Root{99};";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_EQ(tree.name(tree.root())->first, "Root");
    EXPECT_EQ(tree[tree.root()].edge, 99);
}

// 14. Tree with quoted node labels and integer edge labels
TEST_F(NewickParserTest, ParsesQuotedLabelsWithIntegerEdgeLabels)
{
    std::string newick = "('Node A'{100}:0.5,'Node B'{200}:0.6);";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_EQ(tree.name(tree.childId(tree.root(), 0))->first, "Node A");
    EXPECT_EQ(tree[tree.childId(tree.root(), 0)].edge, 100);
    EXPECT_NEAR(tree[tree.childId(tree.root(), 0)].dist, 0.5, 1e-6);
}

// 15. Tree with non-integer edge label (should fail)
TEST_F(NewickParserTest, FailsOnNonIntegerEdgeLabel)
{
    std::string newick = "(A{abc},B{2});";
    maki_tree::Tree tree;
    EXPECT_THROW(tree.parse(newick), maki_tree::newick::NewickParsingError);
}

// 16. Tree with missing closing brace in edge label (should fail)
TEST_F(NewickParserTest, FailsOnMalformedIntegerEdgeLabel)
{
    std::string newick = "(A{1,B{2});";
    maki_tree::Tree tree;
    EXPECT_THROW(tree.parse(newick), maki_tree::newick::NewickParsingError);
}

// 17. Tree with metadata on leaf nodes
TEST_F(NewickParserTest, ParsesMetadataOnLeafNodes)
{
    std::string newick = "(A[&color=red],B[&size=large],C[&confidence=0.95]);";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_EQ(tree[tree.childId(tree.root(), 0)].metadata.at("color"), "red");
    EXPECT_EQ(tree[tree.childId(tree.root(), 1)].metadata.at("size"), "large");
    EXPECT_EQ(tree[tree.childId(tree.root(), 2)].metadata.at("confidence"), "0.95");
}

// 18. Tree with multiple metadata fields
TEST_F(NewickParserTest, ParsesMultipleMetadataFields)
{
    std::string newick = "(A[&color=blue,shape=circle],B[&x=1,y=2]);";
    maki_tree::Tree tree;
    tree.parse(newick);

    const auto &metaA = tree[tree.childId(tree.root(), 0)].metadata;
    EXPECT_EQ(metaA.at("color"), "blue");
    EXPECT_EQ(metaA.at("shape"), "circle");

    const auto &metaB = tree[tree.childId(tree.root(), 1)].metadata;
    EXPECT_EQ(metaB.at("x"), "1");
    EXPECT_EQ(metaB.at("y"), "2");
}

// 19. Tree with metadata on internal node
TEST_F(NewickParserTest, ParsesMetadataOnInternalNode)
{
    std::string newick = "(A,B)Root[&support=0.99];";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_EQ(tree.name(tree.root())->first, "Root");
    EXPECT_EQ(tree[tree.root()].metadata.at("support"), "0.99");
}

// 20. Tree with quoted node labels and metadata
TEST_F(NewickParserTest, ParsesQuotedLabelsWithMetadata)
{
    std::string newick = "('Leaf A'[&note=important],'Leaf B'[&note=check]);";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_EQ(tree.name(tree.childId(tree.root(), 0))->first, "Leaf A");
    EXPECT_EQ(tree[tree.childId(tree.root(), 0)].metadata.at("note"), "important");
    EXPECT_EQ(tree[tree.childId(tree.root(), 1)].metadata.at("note"), "check");
}

// 21. Tree with malformed metadata (missing closing bracket)
TEST_F(NewickParserTest, FailsOnMalformedMetadata)
{
    std::string newick = "(A[&color=red,B);";
    maki_tree::Tree tree;
    EXPECT_THROW(tree.parse(newick), maki_tree::newick::NewickParsingError);
}

// 22. Tree with metadata and edge labels
TEST_F(NewickParserTest, ParsesMetadataWithEdgeLabels)
{
    std::string newick = "(A{1}[&tag=leaf],B{2}[&tag=branch])Root{0}[&type=internal];";
    maki_tree::Tree tree;
    tree.parse(newick);

    EXPECT_EQ(tree[tree.root()].metadata.at("type"), "internal");
    EXPECT_EQ(tree[tree.root()].edge, 0);
    EXPECT_EQ(tree[tree.childId(tree.root(), 0)].edge, 1);
    EXPECT_EQ(tree[tree.childId(tree.root(), 0)].metadata.at("tag"), "leaf");
    EXPECT_EQ(tree[tree.childId(tree.root(), 1)].metadata.at("tag"), "branch");
}

// 23. Deep complex tree with all features
TEST_F(NewickParserTest, ParsesDeepComplexTreeWithAllFeatures)
{
    std::string newick =
        "((('Leaf A'{1}:0.1[&color=red],'Leaf B'{2}:0.2[&color=blue])'Subtree 1'{10}:0.5[&support=0.9],"
        "('Leaf C'{3}:0.3[&color=green],'Leaf D'{4}:0.4[&color=yellow])'Subtree 2'{11}:0.6[&support=0.92])"
        "'Mid Node'{20}:0.8[&type=intermediate],"
        "'Leaf E'{5}:0.7[&color=purple])"
        "'Root Node'{99}:1.5[&support=0.99,type=clade];";

    maki_tree::Tree tree;
    tree.parse(newick);

    // Root node
    EXPECT_EQ(tree.name(tree.root())->first, "Root Node");
    EXPECT_EQ(tree[tree.root()].edge, 99);
    EXPECT_NEAR(tree[tree.root()].dist, 1.5, 1e-6);
    EXPECT_EQ(tree[tree.root()].metadata.at("support"), "0.99");
    EXPECT_EQ(tree[tree.root()].metadata.at("type"), "clade");

    // Mid Node
    auto midNode = tree.childId(tree.root(), 0);
    EXPECT_EQ(tree.name(midNode)->first, "Mid Node");
    EXPECT_EQ(tree[midNode].edge, 20);
    EXPECT_NEAR(tree[midNode].dist, 0.8, 1e-6);
    EXPECT_EQ(tree[midNode].metadata.at("type"), "intermediate");

    // Subtree 1
    auto subtree1 = tree.childId(midNode, 0);
    EXPECT_EQ(tree.name(subtree1)->first, "Subtree 1");
    EXPECT_EQ(tree[subtree1].edge, 10);
    EXPECT_NEAR(tree[subtree1].dist, 0.5, 1e-6);
    EXPECT_EQ(tree[subtree1].metadata.at("support"), "0.9");

    auto cc1 = tree.childId(subtree1, 0);
    EXPECT_EQ(tree.name(cc1)->first, "Leaf A");
    EXPECT_EQ(tree[cc1].edge, 1);
    EXPECT_EQ(tree[cc1].metadata.at("color"), "red");

    auto cc2 = tree.childId(subtree1, 1);
    EXPECT_EQ(tree.name(cc2)->first, "Leaf B");
    EXPECT_EQ(tree[cc2].edge, 2);
    EXPECT_EQ(tree[cc2].metadata.at("color"), "blue");

    // Subtree 2
    auto subtree2 = tree.childId(midNode, 1);
    EXPECT_EQ(tree.name(subtree2)->first, "Subtree 2");
    EXPECT_EQ(tree[subtree2].edge, 11);
    EXPECT_NEAR(tree[subtree2].dist, 0.6, 1e-6);
    EXPECT_EQ(tree[subtree2].metadata.at("support"), "0.92");

    auto cc3 = tree.childId(subtree2, 0);
    EXPECT_EQ(tree.name(cc3)->first, "Leaf C");
    EXPECT_EQ(tree[cc3].edge, 3);
    EXPECT_EQ(tree[cc3].metadata.at("color"), "green");

    auto cc4 = tree.childId(subtree2, 1);
    EXPECT_EQ(tree.name(cc4)->first, "Leaf D");
    EXPECT_EQ(tree[cc4].edge, 4);
    EXPECT_EQ(tree[cc4].metadata.at("color"), "yellow");

    // Leaf E
    auto leafE = tree.childId(tree.root(), 1);
    EXPECT_EQ(tree.name(leafE)->first, "Leaf E");
    EXPECT_EQ(tree[leafE].edge, 5);
    EXPECT_NEAR(tree[leafE].dist, 0.7, 1e-6);
    EXPECT_EQ(tree[leafE].metadata.at("color"), "purple");
}

// 30. countSubleaves on a leaf node should return 0
TEST_F(NewickParserTest, CountSubleavesOnLeafReturnsZero)
{
    std::string newick = "(A,B,C)Root;";
    maki_tree::Tree tree;
    tree.parse(newick);
    EXPECT_EQ(tree.countSubleaves(tree.id("A")->second), 0);
    EXPECT_EQ(tree.countSubleaves(tree.id("B")->second), 0);
    EXPECT_EQ(tree.countSubleaves(tree.id("C")->second), 0);
}

// 31. countSubleaves on internal node with 3 leaves
TEST_F(NewickParserTest, CountSubleavesOnRoot)
{
    std::string newick = "(A,B,C)Root;";
    maki_tree::Tree tree;
    tree.parse(newick);
    EXPECT_EQ(tree.countSubleaves(tree.root()), 3);
}

TEST_F(NewickParserTest, CountSubleavesOnInternalNode)
{
    std::string newick = "(A,((B),(C))Subroot)Root;";
    maki_tree::Tree tree;
    tree.parse(newick);
    EXPECT_EQ(tree.countSubleaves(tree.id("Subroot")->second), 2);
}

// 32. reduce: node with 3 children reduced to 2
TEST_F(NewickParserTest, ReduceNodeWithThreeChildrenToTwo)
{
    maki_tree::Tree tree = autoReduce("(A,B,C)Root;", {"A", "B"});
    EXPECT_EQ(tree.countChildren(tree.root()), 2);
    EXPECT_EQ(tree.name(tree.root())->first, "Root");
    EXPECT_EQ(tree.name(tree.childId(tree.root(), 0))->first, "A");
    EXPECT_EQ(tree.name(tree.childId(tree.root(), 1))->first, "B");
}

// 32. reduce: node with 3 children reduced to 2, preserving distances
TEST_F(NewickParserTest, ReduceNodeWithThreeChildrenToTwoDistances)
{
    maki_tree::Tree tree = autoReduce("(A:0.1,B:0.2,C:0.3)Root;", {"A", "B"});
    EXPECT_EQ(tree.countChildren(tree.root()), 2);
    EXPECT_EQ(tree.name(tree.root())->first, "Root");
    EXPECT_NEAR(tree[tree.childId(tree.root(), 0)].dist, 0.1, 1e-6);
    EXPECT_NEAR(tree[tree.childId(tree.root(), 1)].dist, 0.2, 1e-6);
}

// 33. reduce: unnamed node is collapsed
TEST_F(NewickParserTest, ReduceSingleUnnamedChildCollapse)
{
    maki_tree::Tree tree = autoReduce("((A:0.1,B:0.2):0.3,C:0.5)Root;", {"A", "C"});
    EXPECT_EQ(tree.name(tree.root())->first, "Root");
    EXPECT_EQ(tree.countChildren(tree.root()), 2);
    EXPECT_EQ(tree.name(tree.childId(tree.root(), 0))->first, "A");
    EXPECT_NEAR(tree[tree.childId(tree.root(), 0)].dist, 0.4, 1e-6);
    EXPECT_EQ(tree.name(tree.childId(tree.root(), 1))->first, "C");
    EXPECT_NEAR(tree[tree.childId(tree.root(), 1)].dist, 0.5, 1e-6);
}

// 33. reduce: named node is not collapsed
TEST_F(NewickParserTest, ReduceSingleNamedChildNoCollapse)
{
    maki_tree::Tree tree = autoReduce("((A:0.1,B:0.2)X:0.3,C:0.5)Root;", {"A", "C"});
    EXPECT_EQ(tree.name(tree.root())->first, "Root");
    EXPECT_EQ(tree.countChildren(tree.root()), 2);

    auto c1 = tree.childId(tree.root(), 0);
    EXPECT_EQ(tree.name(c1)->first, "X");
    EXPECT_NEAR(tree[c1].dist, 0.3, 1e-6);
    EXPECT_EQ(tree.countChildren(c1), 1);
    EXPECT_EQ(tree.name(tree.childId(c1, 0))->first, "A");
    EXPECT_NEAR(tree[tree.childId(c1, 0)].dist, 0.1, 1e-6);

    EXPECT_EQ(tree.name(tree.childId(tree.root(), 1))->first, "C");
    EXPECT_NEAR(tree[tree.childId(tree.root(), 1)].dist, 0.5, 1e-6);
}

// 33. reduce: keep subtree
TEST_F(NewickParserTest, RequestKeepSubtree)
{
    maki_tree::Tree tree = autoReduce("((A,B)X,C)Root;", {"X"});
    EXPECT_EQ(tree.name(tree.root())->first, "Root");
    EXPECT_EQ(tree.countChildren(tree.root()), 1);
    EXPECT_EQ(tree.name(tree.childId(tree.root(), 0))->first, "X");
    EXPECT_EQ(tree.countChildren(tree.childId(tree.root(), 0)), 0);
}

// 33. reduce: keep subtree with distances
TEST_F(NewickParserTest, RequestKeepSubtreeWithDistances)
{
    maki_tree::Tree tree = autoReduce("((A:0.2,B:0.3)X:0.1,C:0.4)Root;", {"X"});
    EXPECT_NEAR(tree[tree.childId(tree.root(), 0)].dist, 0.1, 1e-6);
}

// 34. reduce: entire subtree is pruned
TEST_F(NewickParserTest, PruneSubtree)
{
    maki_tree::Tree tree = autoReduce("((A,B)X,(C,D)Y)Root;", {"C", "D"});
    EXPECT_EQ(tree.countChildren(tree.root()), 2);
    EXPECT_EQ(tree.name(tree.root())->first, "Root; Y");
    EXPECT_EQ(tree.name(tree.childId(tree.root(), 0))->first, "C");
    EXPECT_EQ(tree.name(tree.childId(tree.root(), 1))->first, "D");
}

// 34. reduce: entire subtree with distances is pruned
TEST_F(NewickParserTest, PruneSubtreeFixDistances)
{
    maki_tree::Tree tree = autoReduce("((A,B)X,(C:0.3,D:0.4)Y:0.2)Root:0.1;", {"C", "D"});
    EXPECT_EQ(tree.countChildren(tree.root()), 2);
    EXPECT_EQ(tree.name(tree.root())->first, "Root; Y");
    EXPECT_NEAR(tree[tree.root()].dist, 0.3, 1e-6);
    EXPECT_NEAR(tree[tree.childId(tree.root(), 0)].dist, 0.3, 1e-6);
    EXPECT_NEAR(tree[tree.childId(tree.root(), 1)].dist, 0.4, 1e-6);
}

class NewickPreOrderTraversalTest : public ::testing::Test
{
protected:
    NewickPreOrderTraversalTest() : tree() {}

    maki_tree::Tree tree;
};

// 1. Fast pre-order, depth-first traversal for just leaves
TEST_F(NewickPreOrderTraversalTest, ComplexPreOrderTraversalLeaves)
{
    std::string nwk = "(((A,B),(C,(D,E))),((F,G,H),((I,J),(K,(L,M,N)))),O);";
    tree.parse(nwk);

    std::vector<std::string> names;
    std::vector<size_t> depths;

    maki_tree::iter::PreOrderTour tour(tree);
    for (auto it = tour.begin(); it != tour.end(); it.nextNew())
    {
        auto p = *it;
        if (tree.countChildren(p.n.id) == 0)
        {
            names.push_back(tree.name(p.n.id)->first);
            depths.push_back(p.height);
        }
    }

    using _S = std::string;
    EXPECT_THAT(names, testing::ElementsAre(_S("A"), _S("B"), _S("C"), _S("D"), _S("E"), _S("F"), _S("G"), _S("H"), _S("I"), _S("J"), _S("K"), _S("L"), _S("M"), _S("N"), _S("O")));
    EXPECT_THAT(depths, testing::ElementsAre(3, 3, 3, 4, 4, 3, 3, 3, 4, 4, 4, 5, 5, 5, 1));
}

// 2. Euler tour of tree
TEST_F(NewickPreOrderTraversalTest, ComplexPreOrderTraversal)
{
    std::string nwk = "((5,6)2,3,(7)4)1;";
    tree.parse(nwk);

    std::vector<std::string> names;
    std::vector<size_t> depths;

    maki_tree::iter::PreOrderTour tour(tree);
    for (auto it = tour.begin(); it != tour.end(); ++it)
    {
        auto p = *it;
        names.push_back(tree.name(p.n.id)->first);
        depths.push_back(p.height);
    }

    using _S = std::string;
    EXPECT_THAT(names, testing::ElementsAre(_S("1"), _S("2"), _S("5"), _S("2"), _S("6"), _S("2"), _S("1"), _S("3"), _S("1"), _S("4"), _S("7"), _S("4"), _S("1")));
    EXPECT_THAT(depths, testing::ElementsAre(0, 1, 2, 1, 2, 1, 0, 1, 0, 1, 2, 1, 0));
}

class LcaTest : public ::testing::Test
{
protected:
    maki_tree::Tree tree;
    std::unique_ptr<maki_tree::lca::LcaTable> lcaTable;

    void SetUpTree(const std::string &newick)
    {
        tree.parse(newick);
        lcaTable = std::make_unique<maki_tree::lca::LcaTable>(tree);
    }

    size_t GetNodeId(const std::string &name)
    {
        return tree.id(name)->second;
    }

    void ExpectLca(const std::string &node1, const std::string &node2, const std::string &expectedLca)
    {
        size_t id1 = GetNodeId(node1),
               id2 = GetNodeId(node2),
               expectedId = GetNodeId(expectedLca);
        EXPECT_EQ(lcaTable->lca(id1, id2), expectedId);
    }
};

// Test 1: Simple binary tree
TEST_F(LcaTest, SimpleBinaryTree)
{
    SetUpTree("((A,B)E,(C,D)F)root;");
    ExpectLca("A", "B", "E");
    ExpectLca("A", "C", "root");
    ExpectLca("C", "D", "F");
}

// Test 2: Deep nested tree
TEST_F(LcaTest, DeepNestedTree)
{
    SetUpTree("(((A,B)1,(C,(D,E)2)3)4,((F,G,H)5,((I,J)6,(K,(L,M,N)7)8)9)10,O)root;");
    ExpectLca("D", "E", "2");
    ExpectLca("A", "E", "4");
    ExpectLca("L", "N", "7");
    ExpectLca("A", "L", "root"); // Replace with actual root name
}

// Test 3: Unbalanced tree
TEST_F(LcaTest, UnbalancedTree)
{
    SetUpTree("(A,(B,(C,(D,E)1)2)3)root;");
    ExpectLca("D", "E", "1");
    ExpectLca("B", "E", "3");
    ExpectLca("A", "E", "root");
}

// Test 4: Star topology
TEST_F(LcaTest, StarTopology)
{
    SetUpTree("(A,B,C,D,E)root;");
    ExpectLca("A", "B", "root");
    ExpectLca("C", "E", "root");
}

// Test 5: Single node tree
TEST_F(LcaTest, SingleNodeTree)
{
    SetUpTree("(A)root;");
    ExpectLca("A", "A", "A");
}

// Test 6: Complex mixed topology
TEST_F(LcaTest, MixedTopologyTree)
{
    SetUpTree("((A,B)1,((C,D)2,E)3,(F,(G,H)4)5)root;");
    ExpectLca("A", "B", "1");
    ExpectLca("C", "E", "3");
    ExpectLca("G", "H", "4");
    ExpectLca("A", "H", "root");
}

class MultiLcaTest : public LcaTest
{
protected:
    void ExpectLca(const std::vector<std::string> &nodeNames, const std::string &expectedLcaName)
    {
        std::vector<size_t> nodeIds;
        for (const auto &name : nodeNames)
        {
            nodeIds.push_back(GetNodeId(name));
        }
        size_t expectedId = GetNodeId(expectedLcaName);

        std::string qry;
        for (const auto &n : nodeNames)
        {
            qry += n;
            qry += ' ';
        }
        qry += "? " + expectedLcaName;
        qry += " != " + tree.name(lcaTable->lca(nodeIds))->first;
        EXPECT_EQ(lcaTable->lca(nodeIds), expectedId) << qry;
    }
};

// Test 1: Tree with distinct edge labels
TEST_F(MultiLcaTest, EdgeLabeledTreeDistinctNames)
{
    SetUpTree("((leafA,leafB)internalX,(leafC,leafD)internalY)internalZ;");
    ExpectLca({"leafA", "leafB"}, "internalX");
    ExpectLca({"leafC", "leafD"}, "internalY");
    ExpectLca({"leafA", "leafD"}, "internalZ");
    ExpectLca({"leafA", "leafB", "leafC"}, "internalZ");
}

// Test 2: Deep nested tree with distinct labels
TEST_F(MultiLcaTest, DeepNestedEdgeLabeledTreeDistinctNames)
{
    SetUpTree("(((leafA,leafB)internalX,(leafC,(leafD,leafE)internalY)internalW)internalM,"
              "((leafF,leafG,leafH)internalN,((leafI,leafJ)internalO,(leafK,(leafL,leafM,leafN)internalP)internalQ)internalR)internalS)Root;");
    ExpectLca({"leafD", "leafE"}, "internalY");
    ExpectLca({"leafC", "leafE"}, "internalW");
    ExpectLca({"leafA", "leafE"}, "internalM");
    ExpectLca({"leafL", "leafN"}, "internalP");
    ExpectLca({"leafI", "leafJ", "leafK"}, "internalR");
    ExpectLca({"leafA", "leafL"}, "Root");
}

// Test 3: Star topology with internal label
TEST_F(MultiLcaTest, StarTopologyWithDistinctLabel)
{
    SetUpTree("(leafA,leafB,leafC,leafD,leafE)StarRoot;");
    ExpectLca({"leafA", "leafB", "leafC"}, "StarRoot");
    ExpectLca({"leafD", "leafE"}, "StarRoot");
    ExpectLca({"leafA", "leafE"}, "StarRoot");
}

// Test 4: Unbalanced tree with distinct labels
TEST_F(MultiLcaTest, UnbalancedTreeWithDistinctLabels)
{
    SetUpTree("(leafA,(leafB,(leafC,(leafD,leafE)internalX)internalY)internalZ)Root;");
    ExpectLca({"leafD", "leafE"}, "internalX");
    ExpectLca({"leafC", "leafE"}, "internalY");
    ExpectLca({"leafB", "leafE"}, "internalZ");
    ExpectLca({"leafA", "leafE"}, "Root");
}

// Test 5: Single node tree
TEST_F(MultiLcaTest, SingleNodeTree)
{
    SetUpTree("(leafA)root;");
    ExpectLca({"leafA"}, "leafA");
}
