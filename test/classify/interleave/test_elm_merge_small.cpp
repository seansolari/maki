
#include "maki/classify/interleave.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/utils/tempfile.hpp"
#include "test_common.hpp"
#include <filesystem>

/**
Reference
 TACACTNTACTCG

Query
 GACTCA
*/

/*

Graphs
------
kmer edge w+ last

Reference
NNN A 1 0
NNN C 1 0
NNN T 1 1
NNA G 1 1
ACA C 1 1
CGA G 1 1
NTA C 1 1
GTA N 1 1
NNC G 1 1
CAC T 1 1
TAC A 1 0
TAC T 0 1
CTC G 1 1
NAG T 1 1
GAG T 0 1
NCG A 1 1
TCG N 1 1
GTG T 1 1
NNT A 1 1
ACT C 1 1
AGT A 1 0
AGT G 1 1
TGT A 0 1

Query
NNN G 1 0
NNN T 1 1
TCA N 1 1
NGA C 1 1
TGA G 1 1
GAC T 1 1
CTC A 1 1
GTC N 1 1
NNG A 1 1
GAG T 1 1
NTG A 1 1
NNT G 1 1
ACT C 1 1
AGT C 1 1

Nodes
-----

Reference
NNN
NNA
ACA
CGA
NTA
GTA
NNC
CAC
TAC
CTC
NAG
GAG
NCG
TCG
GTG
NNT
ACT
AGT
TGT

Query
NNN
TCA
NGA
TGA
GAC
CTC
GTC
NNG
GAG
NTG
NNT
ACT
AGT

Merged
------
h=1                 R Q     Pos
 0    --N
    NNN             1 0     0
 1    --A
 2    --A
 3    --A
 4    --A
 5    --A
    TCA             1 0     1
    NGA             0 0     0
    TGA             0 0     0
 6    --C
 7    --C
 8    --C
 9    --C
    GAC             1 0     6
    CTC             0 0     0
    GTC             0 0     0
10    --G
11    --G
12    --G
13    --G
14    --G
    NNG             1 0     10
    GAG             0 0     0
    NTG             0 0     0
15    --T
16    --T
17    --T
18    --T
    NNT             1 0     15
    ACT             0 0     0
    AGT             0 0     0
------------------------------

h=2                 R Q     Pos
 0    -NN
    NNN             1 0     0
 1    -NA
 2    -CA
    TCA             1 0     2
 3    -GA
    NGA             1 0     3
    TGA             0 0     0
 4    -TA
 5    -TA
 6    -NC
 7    -AC
 8    -AC
    GAC             1 0     7
 9    -TC
    CTC             1 0     9
    GTC             0 0     0
    NNG             1 1     10
10    -AG
11    -AG
    GAG             1 0     10
12    -CG
13    -CG
14    -TG
    NTG             1 0     14
15    -NT
    NNT             1 0     15
16    -CT
    ACT             1 0     16
17    -GT
18    -GT
    AGT             1 0     17
------------------------------

h=3                 R Q     Pos
 0    NNN
    NNN             1 0     0
 1    NNA
 2    ACA
    TCA             1 1     2
    NGA             1 1     3
 3    CGA
    TGA             1 1     0
 4    NTA
 5    GTA
 6    NNC
 7    CAC
    GAC             1 1     7
 8    TAC
 9    CTC
    CTC             1 0     9
    GTC             0 1     0
    NNG             1 1     10
10    NAG
11    GAG
    GAG             1 0     11
12    NCG
13    TCG
    NTG             1 1     14
14    GTG
15    NNT
    NNT             1 0     15
16    ACT
    ACT             1 0     16
17    AGT
    AGT             1 0     17
18    TGT
------------------------------

*/

class SmallMergeSequences : public testing::Test {
protected:
  SmallMergeSequences()
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

  ~SmallMergeSequences() { fs::remove_all(bufferPath); }

  size_t k;
  fs::path bufferPath, refBuffer, qryBuffer;
  ColouredGraph ref, qry;
};

TEST_F(SmallMergeSequences, InitialRefStructure) {
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

TEST_F(SmallMergeSequences, InitialQryStructure) {
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

TEST_F(SmallMergeSequences, H1InterleavingStructure) {
  ELMMergeSmall Ix(&qry, &ref, 1);

  // expected values
  std::vector<int64_t> expectedZ = {1,  6,  6,  6,  10, 10, 10,
                                    15, 15, 15, 19, 19, 19},
                       expectedP = {0, 1, 0, 0, 6, 0, 0, 10, 0, 0, 15, 0, 0};
  sdsl::bit_vector expectedR = {1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0},
                   expectedQ = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  // check
  const sdsl::bit_vector &actualR = Ix.refOverlap(), &actualQ = Ix.qryOverlap();

  ASSERT_THAT(Ix.previous(), ContainerEq(expectedZ));
  ASSERT_THAT(Ix.refBlocks(), ContainerEq(expectedP));
  ASSERT_THAT(actualR, ContainerEq(expectedR));
  ASSERT_THAT(actualQ, ContainerEq(expectedQ));
}

TEST_F(SmallMergeSequences, H1InterleavingStructureGsize5) {
  ELMMergeSmall Ix(&qry, &ref, 5);

  // expected values
  std::vector<int64_t> expectedZ = {1,  6,  6,  6,  10, 10, 10,
                                    15, 15, 15, 19, 19, 19},
                       expectedP = {0, 1, 0, 0, 6, 0, 0, 10, 0, 0, 15, 0, 0};
  sdsl::bit_vector expectedR = {1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0},
                   expectedQ = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  // check
  const sdsl::bit_vector &actualR = Ix.refOverlap(), &actualQ = Ix.qryOverlap();

  ASSERT_THAT(Ix.previous(), ContainerEq(expectedZ));
  ASSERT_THAT(Ix.refBlocks(), ContainerEq(expectedP));
  ASSERT_THAT(actualR, ContainerEq(expectedR));
  ASSERT_THAT(actualQ, ContainerEq(expectedQ));
}

TEST_F(SmallMergeSequences, LastRefEdge) {
  ELMMergeSmall Ix(&qry, &ref, 1);
  auto e = Ix.refEdge(19);

  ASSERT_EQ(e, 23);
}

TEST_F(SmallMergeSequences, H2InterleavingStructure) {
  ELMMergeSmall Ix(&qry, &ref, 1);
  Ix.DoOneInterleave();

  // expected values
  std::vector<int64_t> expectedZ = {1,  3,  4,  4,  9,  10, 10,
                                    10, 12, 15, 16, 17, 19},
                       expectedP = {0,  2,  3,  0,  7,  9, 0,
                                    10, 10, 14, 15, 16, 17};

  sdsl::bit_vector expectedR = {1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1},
                   expectedQ = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0};

  // check
  const sdsl::bit_vector &actualR = Ix.refOverlap(), &actualQ = Ix.qryOverlap();

  ASSERT_THAT(Ix.previous(), ContainerEq(expectedZ));
  ASSERT_THAT(Ix.refBlocks(), ContainerEq(expectedP));
  ASSERT_THAT(actualR, ContainerEq(expectedR));
  ASSERT_THAT(actualQ, ContainerEq(expectedQ));
}

TEST_F(SmallMergeSequences, H2InterleavingStructureGsize5) {
  ELMMergeSmall Ix(&qry, &ref, 5);
  Ix.DoOneInterleave();

  // expected values
  std::vector<int64_t> expectedZ = {1,  3,  4,  4,  9,  10, 10,
                                    10, 12, 15, 16, 17, 19},
                       expectedP = {0,  2,  3,  0,  7,  9, 0,
                                    10, 10, 14, 15, 16, 17};

  sdsl::bit_vector expectedR = {1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1},
                   expectedQ = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0};

  // check
  const sdsl::bit_vector &actualR = Ix.refOverlap(), &actualQ = Ix.qryOverlap();

  ASSERT_THAT(Ix.previous(), ContainerEq(expectedZ));
  ASSERT_THAT(Ix.refBlocks(), ContainerEq(expectedP));
  ASSERT_THAT(actualR, ContainerEq(expectedR));
  ASSERT_THAT(actualQ, ContainerEq(expectedQ));
}

TEST_F(SmallMergeSequences, H3InterleavingStructure) {
  ELMMergeSmall Ix(&qry, &ref, 1);
  Ix.Interleave();

  // expected values

  std::vector<int64_t> expectedZ = {1,  3,  3,  4,  8,  10, 10,
                                    10, 12, 14, 16, 17, 18};
  sdsl::bit_vector expectedQ = {0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0};

  // check
  const sdsl::bit_vector &actualQ = Ix.qryOverlap();

  ASSERT_THAT(Ix.previous(), ContainerEq(expectedZ));
  ASSERT_THAT(actualQ, ContainerEq(expectedQ));
}

TEST_F(SmallMergeSequences, FullInterleavingBiggerGrainsize) {
  ELMMergeSmall Ix(&qry, &ref, 5);
  Ix.Interleave();

  // expected values

  std::vector<int64_t> expectedZ = {1,  3,  3,  4,  8,  10, 10,
                                    10, 12, 14, 16, 17, 18};
  sdsl::bit_vector expectedQ = {0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0};

  // check
  sdsl::bit_vector const &actualQ = Ix.qryOverlap();

  ASSERT_THAT(Ix.previous(), ContainerEq(expectedZ));
  ASSERT_THAT(actualQ, ContainerEq(expectedQ));
}

/*

Ref 00 NNN A
Ref 01 NNN C
Ref 02 NNN T
Ref 03 NNA G
Ref 04 ACA C
Ref 05 CGA G
Ref 06 NTA C
Ref 07 GTA N
Ref 08 NNC G
Ref 09 CAC T
Ref 10 TAC A
Ref 11 TAC T
Ref 12 CTC G
Ref 13 NAG T
Ref 14 GAG T
Ref 15 NCG A
Ref 16 TCG N
Ref 17 GTG T
Ref 18 NNT A
Ref 19 ACT C
Ref 20 AGT A
Ref 21 AGT G
Ref 22 TGT A

Query
NNN G -> -1
NNN T -> 2
TCA N -> -1
NGA C -> -1
TGA G -> -1
GAC T -> -1
CTC A -> -1
GTC N -> -1
NNG A -> -1
GAG T -> 14
NTG A -> -1
NNT G -> -1
ACT C -> 19
AGT C -> -1

*/

TEST_F(SmallMergeSequences, Classify) {
  std::vector<int64_t> results = ClassifySmall(&qry, &ref, 5),
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

    ELMMergeSmall Ix(&g, &g, 1);
    Ix.Interleave();

    // expected values

    std::vector<int64_t> expectedZ(numQueries);
    std::iota(expectedZ.begin(), expectedZ.end(), 1);

    sdsl::bit_vector expectedQ(numQueries, 0);

    // check
    sdsl::bit_vector const &actualQ = Ix.qryOverlap();

    ASSERT_THAT(Ix.previous(), ContainerEq(expectedZ));
    ASSERT_THAT(actualQ, ContainerEq(expectedQ));
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

    size_t numQueries = qry.nodes();

    ELMMergeSmall Ix(&qry, &ref, 1);
    Ix.Interleave();

    // expected values
    sdsl::bit_vector expectedQ(numQueries, 1);
    expectedQ[0] = 0;

    sdsl::bit_vector const &actualQ = Ix.qryOverlap();
    ASSERT_THAT(actualQ, ContainerEq(expectedQ));
  }

  fs::path bufferPath, refBuffer, qryBuffer;
  std::string seqCG, seqAT;
  ColouredGraph ref, qry;
};

TEST_F(DisjointSequences, k3) { doTest(3); }

TEST_F(DisjointSequences, k5) { doTest(5); }

TEST_F(DisjointSequences, k11) { doTest(11); }
