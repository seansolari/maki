
#include "maki/core/utils/tempfile.hpp"
#include "test_common.hpp"

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
  
}