#include "maki/core/graph/cdbg.hpp"
#include "test_common.hpp"

#include "maki/core/seq/random.hpp"
#include "maki/core/utils/tempfile.hpp"

// ==== Tests
// ========================================================================

class TraversableGraph_Kmer : public testing::Test {
protected:
  TraversableGraph_Kmer() : bufferPath(tempio::create_temporary_directory()) {}

  ~TraversableGraph_Kmer() { fs::remove_all(bufferPath); }

  fs::path bufferPath;
};

TEST_F(TraversableGraph_Kmer, SimplePath_NoRepeats) {
  std::string fasta = BuildFasta({"ACGTACGT"});
  size_t k = 4;
  BuildAndCheckGraph(fasta, k, 1, bufferPath);
}

TEST_F(TraversableGraph_Kmer, Homopolymer_Repeats) {
  std::string fasta = BuildFasta({"AAAAAAAAAA"}); // 10 As
  for (size_t k = 4; k <= 9; ++k) {
    BuildAndCheckGraph(fasta, k, 1, bufferPath);
  }
}

TEST_F(TraversableGraph_Kmer, Branching) {
  std::string fasta = BuildFasta({"ACAAAA", "ACAGAA"});
  size_t k = 4;
  BuildAndCheckGraph(fasta, k, 1, bufferPath);
}

TEST_F(TraversableGraph_Kmer, Palindromic) {
  std::string fasta = BuildFasta({"ATGCAT"});
  size_t k = 4;
  BuildAndCheckGraph(fasta, k, 1, bufferPath);
}

TEST_F(TraversableGraph_Kmer, KEqualsSequenceLength) {
  std::string fasta = BuildFasta({"ACGT"});
  size_t k = 4;
  BuildAndCheckGraph(fasta, k, 1, bufferPath);
}

TEST_F(TraversableGraph_Kmer, KGreaterThanSequenceLength_ShouldYieldNoEdges) {
  std::string fasta = BuildFasta({"ACGT"});
  size_t k = 5;
  MakeGraph(fasta, k, 1, bufferPath);
  ColouredGraph g;
  ColouredGraph::FromDisk(g, bufferPath);
  ASSERT_EQ(edge_count(g), 0u)
      << "Graph should have 0 edges when k > sequence length.";
}

TEST_F(TraversableGraph_Kmer, RandomizedFuzz_Small) {
  std::mt19937 rng(123456); // reproducible
  std::uniform_int_distribution<int> numSeqDist(1, 3);
  std::uniform_int_distribution<int> lenDist(10, 100);

  for (int iter = 0; iter < 50; ++iter) {
    int numSeqs = numSeqDist(rng);
    auto [fasta, minLen] = RandomFasta(numSeqs, rng, lenDist);

    std::uniform_int_distribution<size_t> kDist(4,
                                                std::min<size_t>(10, minLen));
    size_t k = kDist(rng);

    BuildAndCheckGraph(fasta, k, 1, bufferPath);
  }
}

TEST_F(TraversableGraph_Kmer, Stress_LongSequence_ModerateK) {
  std::mt19937 rng(424242);
  std::string seq = RandomDNA(2000, rng);
  std::string fasta = BuildFasta({seq});
  size_t k = 21;

  BuildAndCheckGraph(fasta, k, 1, bufferPath);
}

/*

   Node    Edge    W(-)    W(+)
 0 NNN     G       1       0
 1 NNN     T       1       1
 2 TCA     N       1       1
 3 NGA     C       1       1
 4 TGA     G       1       1
 5 GAC     T       1       1
 6 CTC     A       1       1
 7 GTC     N       1       1
 8 NNG     A       1       1
 9 GAG     T       1       1
10 NTG     A       1       1
11 NNT     G       1       1
12 ACT     C       1       1
13 AGT     C       1       1

*/
class GraphTraverseTests : public testing::Test {
protected:
  GraphTraverseTests()
      : k(4), bufferPath(tempio::create_temporary_directory()), g() {
    MakeGraph(">contig_one\nGACTCA\n>contig_two\nTGAGTC\n", k - 1, 1,
              bufferPath);
    ColouredGraph::FromDisk(g, bufferPath);
  }

  ~GraphTraverseTests() { fs::remove_all(bufferPath); }

  uint8_t k;
  fs::path bufferPath;
  ColouredGraph g;
};

TEST_F(GraphTraverseTests, Fwd0) {
  auto i = g.fwd(0);
  ASSERT_EQ(*i, 8);
}

TEST_F(GraphTraverseTests, Fwd1) {
  auto i = g.fwd(1);
  ASSERT_EQ(*i, 11);
}

TEST_F(GraphTraverseTests, Fwd4) {
  auto i = g.fwd(4);
  ASSERT_EQ(*i, 9);
}

TEST_F(GraphTraverseTests, NoFwd) {
  auto i = g.fwd(2);
  ASSERT_FALSE(i.has_value());
}

TEST_F(GraphTraverseTests, Bwd9) {
  auto i = g.bwd(9);
  ASSERT_EQ(*i, 4);
}

TEST_F(GraphTraverseTests, BwdRoot) {
  auto i = g.bwd(8);
  ASSERT_EQ(*i, 0);
}

TEST_F(GraphTraverseTests, BwdNearRoot) {
  auto i = g.bwd(11);
  ASSERT_EQ(*i, 1);
}

TEST_F(GraphTraverseTests, Bwd0) {
  auto i = g.bwd(0);
  ASSERT_FALSE(i.has_value());
}

TEST_F(GraphTraverseTests, Kmer6) {
  auto s = g.kmer(6);
  ASSERT_EQ(s, Dna4Sequence("CTCA"_dna4));
}

TEST_F(GraphTraverseTests, Kmer7) {
  auto s = g.kmer(7);
  ASSERT_EQ(s, Dna4Sequence("GTC"_dna4));
}

TEST_F(GraphTraverseTests, Kmer8) {
  auto s = g.kmer(8);
  ASSERT_EQ(s, Dna4Sequence("GA"_dna4));
}

TEST_F(GraphTraverseTests, Kmer0) {
  auto s = g.kmer(0);
  ASSERT_EQ(s, Dna4Sequence("G"_dna4));
}

/*

   Node    Edge    W(-)    W(+)
 0 NNN     A       1       0
 1 NNN     C       1       1
 2 NNA     C       1       1
 3 CGA     C       1       1
 4 NNC     C       1       1
 5 NAC     G       1       1
 6 GAC     C       1       1
 7 NCC     G       1       1
 8 ACC     G       0       1
 9 GTC     G       1       1
10 ACG     G       1       1
11 CCG     A       1       0
12 CCG     T       1       1
13 TCG     G       0       1
14 CGG     T       1       1
15 GGT     C       1       1

*/
class GraphTraverseRepeatTests : public testing::Test {
protected:
  GraphTraverseRepeatTests()
      : k(4), bufferPath(tempio::create_temporary_directory()), g() {
    MakeGraph(">contig_one\nACGGTCGG\n", k - 1, 1, bufferPath);
    ColouredGraph::FromDisk(g, bufferPath);
  }

  ~GraphTraverseRepeatTests() { fs::remove_all(bufferPath); }

  uint8_t k;
  fs::path bufferPath;
  ColouredGraph g;
};

TEST_F(GraphTraverseRepeatTests, Fwd7) {
  auto i = g.fwd(7);
  ASSERT_EQ(*i, 12);
}

TEST_F(GraphTraverseRepeatTests, Fwd8) {
  auto i = g.fwd(8);
  ASSERT_EQ(*i, 12);
}

TEST_F(GraphTraverseRepeatTests, Bwd12) {
  auto i = g.bwd(12);
  ASSERT_EQ(*i, 7);
}

TEST_F(GraphTraverseRepeatTests, Bwd11) {
  auto i = g.bwd(11);
  ASSERT_EQ(*i, 7);
}

TEST_F(GraphTraverseRepeatTests, Kmer8) {
  auto s = g.kmer(8);
  ASSERT_EQ(s, Dna4Sequence("ACCG"_dna4));
}

TEST_F(GraphTraverseRepeatTests, Kmer8Str) {
  auto s = toString(g.kmer(8));
  ASSERT_EQ(s, std::string("ACCG"));
}

class TraversableGraph_Colour : public testing::Test {
protected:
  TraversableGraph_Colour()
      : bufferPath(tempio::create_temporary_directory()) {}

  ~TraversableGraph_Colour() { fs::remove_all(bufferPath); }

  fs::path bufferPath;
};

/**
 * Collect requested colours from graph
 */
TEST_F(TraversableGraph_Colour, RandomSequences) {
  // Generate random genomes
  std::mt19937 rng(123456);
  std::uniform_int_distribution<int> numSeqDist(5, 20);
  std::uniform_int_distribution<int> lenDist(100u, 10000u);
  int numSeqs = numSeqDist(rng);
  auto [fasta, _] = RandomFasta(numSeqs, rng, lenDist);

  // Build graph
  std::size_t k = 4u;
  BuildAndCheckBufferEdgePositions(fasta, k, 1u, bufferPath);
}
