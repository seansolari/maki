
#include "maki/classify/interleave.hpp"
#include "maki/core/utils/tempfile.hpp"
#include "test_common.hpp"
#include <sdsl/int_vector.hpp>

class LargeMergeSequences : public testing::Test {
protected:
  LargeMergeSequences()
      : k(3), bufferPath(tempio::create_temporary_directory()),
        refBuffer(bufferPath / "ref"), qryBuffer(bufferPath / "qry"), ref(),
        qry() {
    fs::create_directories(refBuffer);
    fs::create_directories(qryBuffer);

    MakeGraph(">contig_one\nTACACTNTACTCG\n", k, 1, refBuffer);
    ColouredGraph::FromDisk(ref, refBuffer);

    MakeGraph(">contig_two\nGACTCA\n", k, 1, qryBuffer);
    ColouredGraph::FromDisk(qry, qryBuffer);
  }

  ~LargeMergeSequences() { fs::remove_all(bufferPath); }

  size_t k;
  fs::path bufferPath, refBuffer, qryBuffer;
  ColouredGraph ref, qry;
};

TEST_F(LargeMergeSequences, InitialRefStructure) {
  ASSERT_EQ(ref.edges(), 23);
  ASSERT_EQ(ref.nodes(), 19);

  sdsl::int_vector<4> W = {0b1001, 0b1010, 0b1100, 0b1011, 0b1010, 0b1011,
                           0b1010, 0b1000, 0b1011, 0b1100, 0b1001, 0b0100,
                           0b1011, 0b1100, 0b0100, 0b1001, 0b1000, 0b1100,
                           0b1001, 0b1010, 0b1001, 0b1011, 0b0001};
  for (size_t i = 0; i < ref.edges(); ++i)
    ASSERT_EQ(ref.W[i], W[i]) << "mismatch at position " << i;

  sdsl::bit_vector last = {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1,
                           1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1};
  ASSERT_THAT(ref.l, ContainerEq(last));

  ASSERT_THAT(ref.F, ElementsAreArray({0, 3, 8, 13, 18}));
  ASSERT_THAT(ref.C, ElementsAreArray({0, 1, 6, 10, 15}));
}

TEST_F(LargeMergeSequences, InitialQryStructure) {
  ASSERT_EQ(qry.edges(), 14);
  ASSERT_EQ(qry.nodes(), 13);

  sdsl::int_vector<4> W = {0b1011, 0b1100, 0b1000, 0b1010, 0b1011,
                           0b1100, 0b1001, 0b1000, 0b1001, 0b1100,
                           0b1001, 0b1011, 0b1010, 0b1010};
  for (size_t i = 0; i < qry.edges(); ++i)
    ASSERT_EQ(qry.W[i], W[i]) << "mismatch at position " << i;

  sdsl::bit_vector last = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  ASSERT_THAT(qry.l, ContainerEq(last));

  ASSERT_THAT(qry.F, ElementsAreArray({0, 2, 5, 8, 11}));
  ASSERT_THAT(qry.C, ElementsAreArray({0, 1, 4, 7, 10}));
}

TEST_F(LargeMergeSequences, H1InterleavingStructure) {
  ELMMergeLarge Ix(&qry, &ref, 1);

  sdsl::bit_vector expectedZ = {0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1,
                                1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1},
                   expectedB = {1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                                0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0};
  ASSERT_THAT(Ix.previous(), ContainerEq(expectedZ));
  ASSERT_THAT(Ix.b.Bp, ContainerEq(expectedB));
}

TEST_F(LargeMergeSequences, H1InterleavingStructureGsize5) {
  ELMMergeLarge Ix(&qry, &ref, 5);

  sdsl::bit_vector expectedZ = {0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1,
                                1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1},
                   expectedB = {1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                                0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0};
  ASSERT_THAT(Ix.previous(), ContainerEq(expectedZ));
  ASSERT_THAT(Ix.b.Bp, ContainerEq(expectedB));
}

TEST_F(LargeMergeSequences, H2InterleavingStructure) {
  ELMMergeLarge Ix(&qry, &ref, 1);
  Ix.DoOneDummyInterleave();

  sdsl::bit_vector expectedZ = {0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1,
                                1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1},
                   expectedB = {1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0,
                                0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0};
  ASSERT_THAT(Ix.previous(), ContainerEq(expectedZ));
  ASSERT_THAT(Ix.b.Bp, ContainerEq(expectedB));
}

TEST_F(LargeMergeSequences, H2InterleavingStructureGsize5) {
  ELMMergeLarge Ix(&qry, &ref, 5);
  Ix.DoOneInterleave();

  sdsl::bit_vector expectedZ = {0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1,
                                1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1},
                   expectedB = {1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0,
                                0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0};
  ASSERT_THAT(Ix.previous(), ContainerEq(expectedZ));
  ASSERT_THAT(Ix.b.Bp, ContainerEq(expectedB));
}

TEST_F(LargeMergeSequences, H3InterleavingStructure) {
  ELMMergeLarge Ix(&qry, &ref, 1);
  Ix.Interleave();

  sdsl::bit_vector expectedZ = {0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1,
                                1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0},
                   expectedB = {1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
                                1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1};
  ASSERT_THAT(Ix.previous(), ContainerEq(expectedZ));
  ASSERT_THAT(Ix.b.Bp, ContainerEq(expectedB));
}

TEST_F(LargeMergeSequences, FullInterleavingBiggerGrainsize) {
  ELMMergeLarge Ix(&qry, &ref, 5);
  Ix.Interleave();

  sdsl::bit_vector expectedZ = {0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1,
                                1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0},
                   expectedB = {1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
                                1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1};
  ASSERT_THAT(Ix.previous(), ContainerEq(expectedZ));
  ASSERT_THAT(Ix.b.Bp, ContainerEq(expectedB));
}

TEST_F(LargeMergeSequences, Classify) {
  std::vector<int64_t> results = ClassifyLarge(&qry, &ref, 5),
                       expected = {-1, 2,  -1, -1, -1, -1, -1,
                                   -1, -1, 14, -1, -1, 19, -1};
  ASSERT_THAT(results, ContainerEq(expected));
}

class IdenticalSequences : public testing::Test {
protected:
  IdenticalSequences()
      : bufferBase(tempio::create_temporary_directory()),
        sequence(">contig\nACGTAGCTGCATGCTAGCTAGTCAGTCGATCTGCAGCTA\n") {}

  ~IdenticalSequences() { fs::remove_all(bufferBase); }

  void doTest(size_t k) {
    MakeGraph(sequence, k, 1, bufferBase);

    ColouredGraph g;
    ColouredGraph::FromDisk(g, bufferBase);
    size_t numQueries = g.nodes();

    ELMMergeLarge Ix(&g, &g, 1);
    Ix.Interleave();

    // expected values

    sdsl::bit_vector expectedZ(2 * numQueries, 0);
    for (std::size_t i = 0; i < numQueries; ++i)
      expectedZ[(2 * i) + 1] = 1;

    sdsl::bit_vector expectedQ(2 * numQueries, 0);
    for (std::size_t i = 0; i < numQueries; ++i)
      expectedQ[2 * i] = 1;

    // check

    ASSERT_THAT(Ix.previous(), ContainerEq(expectedZ));
    ASSERT_THAT(Ix.b.Bp, ContainerEq(expectedQ));
  }

  fs::path bufferBase;
  std::string sequence;
};

TEST_F(IdenticalSequences, k3) { doTest(3); }

TEST_F(IdenticalSequences, k5) { doTest(5); }

TEST_F(IdenticalSequences, k11) { doTest(11); }

class DisjointSequences : public testing::Test {
protected:
  DisjointSequences()
      : bufferPath(tempio::create_temporary_directory()),
        refBuffer(bufferPath / "ref"), qryBuffer(bufferPath / "qry"),
        seqCG(">contig_one\nCGCGCGCGGGGGCGCGCGCGGCGGCGCGCGCGCGGGC\n"),
        seqAT(">contig_two\nATATATTATTTTAAAAAAAAAAATATATTAAAAAAA\n") {
    fs::create_directories(refBuffer);
    fs::create_directories(qryBuffer);
  }

  ~DisjointSequences() { fs::remove_all(bufferPath); }

  void doTest(size_t k) {
    MakeGraph(seqCG, k, 1, refBuffer);
    ColouredGraph::FromDisk(ref, refBuffer);

    MakeGraph(seqAT, k, 1, qryBuffer);
    ColouredGraph::FromDisk(qry, qryBuffer);

    size_t numQueries = qry.nodes(), numReferences = ref.nodes();

    ELMMergeLarge Ix(&qry, &ref, 1);
    Ix.Interleave();

    // expected values
    sdsl::bit_vector expectedQ(numQueries + numReferences, 1);
    expectedQ[1] = 0;

    ASSERT_THAT(Ix.b.Bp, ContainerEq(expectedQ));
  }

  fs::path bufferPath, refBuffer, qryBuffer;
  std::string seqCG, seqAT;
  ColouredGraph ref, qry;
};

TEST_F(DisjointSequences, k3) { doTest(3); }

TEST_F(DisjointSequences, k5) { doTest(5); }

TEST_F(DisjointSequences, k11) { doTest(11); }
