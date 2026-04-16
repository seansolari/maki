#include "maki/build/graph/build_colours.hpp"
#include "maki/build/graph/interleave_buffers.hpp"
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using ::testing::ContainerEq;
using ::testing::ElementsAreArray;

class GraphInsertTests : public testing::Test {
protected:
  GraphInsertTests() : k(7), num_colours(11) {}

  uint8_t k, num_colours;
};

TEST_F(GraphInsertTests, WriteColouredEdges) {
  // Prepare data
  packet pkt{};
  ++pkt.block;

  pkt.data.emplace_back(BufferValue((size_t)1ull, 0ull));
  pkt.data.emplace_back(BufferValue((size_t)2ull, 0ull));
  pkt.data.emplace_back(BufferValue((size_t)3ull, 1ull));
  pkt.data.emplace_back(BufferValue((size_t)4ull, 1ull));
  pkt.data.emplace_back(BufferValue((size_t)5ull, 2ull));
  pkt.data.emplace_back(BufferValue((size_t)6ull, 3ull));
  pkt.data.emplace_back(BufferValue((size_t)7ull, 3ull));
  pkt.data.emplace_back(BufferValue((size_t)8ull, 3ull));
  pkt.data.emplace_back(
      BufferValue((size_t)0ull,
                  KmerBuffer::terminalEdge)); // should be ignored
  pkt.data.emplace_back(
      BufferValue((size_t)0ull,
                  KmerBuffer::terminalEdge)); // should be ignored

  // Push to graph
  sdsl::int_vector<4> edges;
  sdsl::bit_vector succ;
  ArchivePayload carch;
  MetaColours cmap(8);

  pushNode(pkt, edges, succ, carch, cmap, 1u /* A */);

  // Check graph structure
  EXPECT_EQ(pkt.block, 0);

  sdsl::int_vector<4> Xedges = {0b1001, 0b1010, 0b1011, 0b1100};
  sdsl::bit_vector Xwplus = {0, 0, 0, 1};

  EXPECT_THAT(edges, ContainerEq(Xedges));
  EXPECT_THAT(succ, ContainerEq(Xwplus));

  // Check colours
  EXPECT_THAT(carch.raw.view(), ElementsAreArray({9, 10, 5, 11}));
  EXPECT_EQ(cmap.id({1, 2}), 9);
  EXPECT_EQ(cmap.id({3, 4}), 10);
  EXPECT_EQ(cmap.id({6, 7, 8}), 11);
}
