#include "maki/build/graph/build_colours.hpp"
#include "maki/build/graph/interleave_buffers.hpp"

#include "test_common.hpp"

class GraphInsertTests : public testing::Test {
protected:
  GraphInsertTests() : k(7), num_colours(11) {}

  uint8_t k, num_colours;
};

TEST_F(GraphInsertTests, WriteColouredEdges) {
  // Prepare data
  packet pkt{};
  ++pkt.block;

  pkt.emplace(BufferValue((size_t)8ull, 3ull));
  pkt.emplace(BufferValue((size_t)1ull, 0ull));
  pkt.emplace(BufferValue((size_t)7ull, 3ull));
  pkt.emplace(BufferValue((size_t)3ull, 1ull));
  pkt.emplace(BufferValue((size_t)2ull, 0ull));
  pkt.emplace(BufferValue((size_t)0ull,
                          KmerBuffer::terminalEdge)); // should be ignored
  pkt.emplace(BufferValue((size_t)5ull, 2ull));
  pkt.emplace(BufferValue((size_t)4ull, 1ull));
  pkt.emplace(BufferValue((size_t)6ull, 3ull));
  pkt.emplace(BufferValue((size_t)0ull,
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
  std::vector<uint64_t> cols = {9, 10, 5, 11};
  EXPECT_THAT(carch.raw.view(), ElementsAreArray(cols));
  EXPECT_EQ(cmap.id({1, 2}), 9);
  EXPECT_EQ(cmap.id({3, 4}), 10);
  EXPECT_EQ(cmap.id({6, 7, 8}), 11);
}

TEST_F(GraphInsertTests, InsertColouredKmers) {
  // Create k-mers to insert
  KmerBuffer buffer(value_size(ceil_log2(num_colours)), k, k,
                    {
                        {0b01111011, 0b00110010, 0b00000001}, // TGTCGAT -> C
                        {0b01111011, 0b00110010, 0b00001001}, // TGTCGAT -> C
                        {0b10010000, 0b00111000, 0b00000011}, // AACGAGT -> T
                        {0b10010000, 0b00111000, 0b00011011}, // AACGAGT -> T
                        {0b10010001, 0b00111000, 0b00000000}, // CACGAGT -> A
                        {0b10010001, 0b00111000, 0b00111000}, // CACGAGT -> A
                        {0b10010011, 0b00111000, 0b00001000}, // TACGAGT -> A
                        {0b10010011, 0b00111000, 0b00101000}, // TACGAGT -> A
                        {0b01110000, 0b00111010, 0b00000001}, // AATCGGT -> C
                        {0b01110011, 0b00111010, 0b00000001}, // TATCGGT -> C
                    });
  sdsl::int_vector<2> B = {BW_0_K, IS_0, BW_0_K, IS_0,   IS_K,
                           IS_0,   IS_K, IS_0,   BW_0_K, IS_K};

  // Push to graph
  sdsl::int_vector<4> edges;
  sdsl::bit_vector succ;
  ArchivePayload carch;
  MetaColours cmap(7);

  packet pkt{};
  pushNodes(pkt, buffer.begin(), buffer.end(), B.begin(), edges, succ, carch,
            cmap, 3u /* T */);

  // Check graph structure
  sdsl::int_vector<4> Xedges = {0b1010, 0b1100, 0b1001, 0b0001, 0b1010, 0b0010};
  sdsl::bit_vector Xwplus = {1, 1, 1, 1, 1, 1};
  EXPECT_THAT(edges, ContainerEq(Xedges));
  EXPECT_THAT(succ, ContainerEq(Xwplus));

  // Check colours
  std::vector<uint64_t> cols = {8, 9, 10, 11, 0, 0};
  EXPECT_THAT(carch.raw.view(), ElementsAreArray(cols));
  EXPECT_EQ(cmap.id({0, 1}), 8);
  EXPECT_EQ(cmap.id({0, 3}), 9);
  EXPECT_EQ(cmap.id({0, 7}), 10);
  EXPECT_EQ(cmap.id({1, 5}), 11);

  // Check blocks
  EXPECT_THAT(pkt.str.C, ElementsAreArray({0, 0, 0, 0, 6}));
  EXPECT_THAT(pkt.str.F, ElementsAreArray({0, 0, 0, 0, 6}));
}

/**
 * $$$$ -$$$ -> A
 * $$$$ -$$$ -> C
 * $$$$ -AAA -> G
 * AAA$ -AAA -> G
 * AAC$ -AAA -> T
 * $$$$ -TAA -> A
 * $$$$ -TAA -> C
 *
 */
TEST_F(GraphInsertTests, InsertTerminals) {
  TerminalBuffer terminals(k,
                           {{/* size */ 0b00000000, /* kmer */ 0b00000000,
                             0b00000000, /* edge */ 0b00000000},
                            {/* size */ 0b00000000, /* kmer */ 0b00000000,
                             0b00000000, /* edge */ 0b00000001},
                            {/* size */ 0b00000011, /* kmer */ 0b00000000,
                             0b00000000, /* edge */ 0b00000010},
                            {/* size */ 0b00000110, /* kmer */ 0b00000000,
                             0b00000000, /* edge */ 0b00000010},
                            {/* size */ 0b00000110, /* kmer */ 0b00000100,
                             0b00000000, /* edge */ 0b00000011},
                            {/* size */ 0b00000011, /* kmer */ 0b00000000,
                             0b00110000, /* edge */ 0b00000000},
                            {/* size */ 0b00000011, /* kmer */ 0b00000000,
                             0b00110000, /* edge */ 0b00000001}},
                           TerminalBuffer::autofit_tag);
  sdsl::int_vector<2> B = {BW_0_K, IS_0, BW_0_K, BW_0_K, BW_0_K, BW_0_K, IS_0};

  // Push to graph
  sdsl::int_vector<4> edges;
  sdsl::bit_vector succ;
  ArchivePayload carch;
  MetaColours cmap(7);

  auto view = terminals.asRange();
  auto ps = pushRange(view, B.begin(), edges, succ, carch, cmap, 3u /* T */);

  // Check graph structure
  sdsl::int_vector<4> Xedges = {0b1001, 0b1010, 0b1011, 0b1011,
                                0b1100, 0b1001, 0b1010};
  sdsl::bit_vector Xwplus = {0, 1, 1, 1, 1, 0, 1};
  EXPECT_THAT(edges, ContainerEq(Xedges));
  EXPECT_THAT(succ, ContainerEq(Xwplus));

  // Check colours
  std::vector<uint64_t> cols(7, 0);
  EXPECT_THAT(carch.raw.view(), ElementsAreArray(cols));

  // Check blocks
  EXPECT_THAT(ps.C, ElementsAreArray({0, 0, 0, 0, 5}));
  EXPECT_THAT(ps.F, ElementsAreArray({0, 0, 0, 0, 7}));
}

class FlushCompareTests : public testing::Test {
protected:
  FlushCompareTests()
      : k(10), k_eff(8), num_colours(7),
        kmers(value_size(ceil_log2(num_colours)), k, k_eff,
              {{/* kmer */ 0b01111011, 0b00001010, /* edge */ 0b00000001},
               {/* kmer */ 0b01111011, 0b00001010, /* edge */ 0b00001001},
               {/* kmer */ 0b10010000, 0b00011000, /* edge */ 0b00011011},
               {/* kmer */ 0b10010001, 0b00011000, /* edge */ 0b00000000},
               {/* kmer */ 0b10010001, 0b00011000, /* edge */ 0b00111000},
               {/* kmer */ 0b10010011, 0b00011000, /* edge */ 0b00001000},
               {/* kmer */ 0b10010011, 0b00011000, /* edge */ 0b00101000},
               {/* kmer */ 0b01110011, 0b00111010, /* edge */ 0b00000001}}),
        terminals(k,
                  {{/* size */ 0b00000010, /* kmer */ 0b00000000, 0b00000000,
                    0b00001111, /* edge */ 0b00000000},
                   {/* size */ 0b00000010, /* kmer */ 0b00000000, 0b00000000,
                    0b00001111, /* edge */ 0b00000010},
                   {/* size */ 0b00000010, /* kmer */ 0b00000000, 0b00000000,
                    0b00001111, /* edge */ 0b00000011},
                   {/* size */ 0b00000011, /* kmer */ 0b00000000, 0b00000000,
                    0b00001111, /* edge */ 0b00000000},
                   {/* size */ 0b00000110, /* kmer */ 0b00000000, 0b00000010,
                    0b00001111, /* edge */ 0b00000000},
                   {/* size */ 0b00000101, /* kmer */ 0b00000000, 0b00011000,
                    0b00001111, /* edge */ 0b00000000},
                   {/* size */ 0b00001000, /* kmer */ 0b10010000, 0b00011000,
                    0b00001111, /* edge */ 0b00000011},
                   {/* size */ 0b00001001, /* kmer */ 0b01110000, 0b00111010,
                    0b00001111, /* edge */ 0b00000001},
                   {/* size */ 0b00000011, /* kmer */ 0b00000000, 0b01000000,
                    0b00001111, /* edge */ 0b00000011}},
                  TerminalBuffer::autofit_tag),
        LessThan(k_eff, terminals.lengthBytes()) {}

  uint8_t k, k_eff;
  size_t num_colours;
  KmerBuffer kmers;
  TerminalBuffer terminals;
  MaskedAlignedBytesLessThan<RHS> LessThan;
};

/**
 * MaskedBytesNotEqual used to compared terminals (LHS, masked) to k-mers (RHS,
 * unmasked).
 *
 */
TEST_F(FlushCompareTests, MaskedBytesDiff_4k_bw0k) {
  MaskedBytesDiff Diff{8};
  std::vector<uint8_t> tmn = {/* size */ 4, /* kmer */ 0b00000000, 0b00010111,
                              /* edge */ 0b00000000},
                       kmr = {/* kmer */ 0b11000000, 0b00010111,
                              /* edge */ 0b00001000};
  EXPECT_EQ(Diff(tmn.data() + 1, kmr.data()), BW_0_K);
}

TEST_F(FlushCompareTests, MaskedBytesDiff_4k_isk) {
  MaskedBytesDiff Diff{8};
  std::vector<uint8_t> tmn = {/* size */ 4, /* kmer */ 0b00000000, 0b00010111,
                              /* edge */ 0b00000000},
                       kmr = {/* kmer */ 0b00000001, 0b00010111,
                              /* edge */ 0b00001000};
  EXPECT_EQ(Diff(tmn.data() + 1, kmr.data()), IS_K);
}

TEST_F(FlushCompareTests, MaskedBytesDiff_4k_is0) {
  MaskedBytesDiff Diff{8};
  std::vector<uint8_t> tmn = {/* size */ 7, /* kmer */ 0b11100000, 0b01101100,
                              /* edge */ 0b00000000},
                       kmr = {/* kmer */ 0b11100000, 0b01101100,
                              /* edge */ 0b00001000};
  EXPECT_EQ(Diff(tmn.data() + 1, kmr.data()), IS_0);
}

TEST_F(FlushCompareTests, MaskedBytesDiff_Truncated_bw0k) {
  MaskedBytesDiff Diff{7};
  std::vector<uint8_t> tmn = {/* size */ 10, /* kmer */ 0b10010011, 0b11001001,
                              0b00001111, /* edge */ 0b00000000},
                       kmr = {/* kmer */ 0b00010011, 0b00001001,
                              /* edge */ 0b00001000};
  EXPECT_EQ(Diff(tmn.data() + 1, kmr.data()), BW_0_K);
}

TEST_F(FlushCompareTests, MaskedBytesDiff_Truncated_isk) {
  MaskedBytesDiff Diff{7};
  std::vector<uint8_t> tmn = {/* size */ 10, /* kmer */ 0b10010011, 0b11001001,
                              0b00001111, /* edge */ 0b00000000},
                       kmr = {/* kmer */ 0b10010010, 0b00001001,
                              /* edge */ 0b00001000};
  EXPECT_EQ(Diff(tmn.data() + 1, kmr.data()), IS_K);
}

TEST_F(FlushCompareTests, MaskedBytesDiff_Truncated_is0) {
  MaskedBytesDiff Diff{7};
  std::vector<uint8_t> tmn = {/* size */ 10, /* kmer */ 0b10010011, 0b11001001,
                              0b00001111, /* edge */ 0b00000000},
                       kmr = {/* kmer */ 0b10010011, 0b00001001,
                              /* edge */ 0b00001000};
  EXPECT_EQ(Diff(tmn.data() + 1, kmr.data()), IS_0);
}

TEST_F(FlushCompareTests, LowerBound_KmerQuery_Seen) {
  auto it = std::lower_bound(kmers.begin(), kmers.end(), *terminals.constAt(6),
                             LessThan);
  ASSERT_EQ(it, kmers.at(2));
}

TEST_F(FlushCompareTests, LowerBound_KmerQuery_Begin) {
  auto it = std::lower_bound(kmers.begin(), kmers.end(),
                             *terminals.constBegin(), LessThan);
  ASSERT_EQ(it, kmers.begin());
}

TEST_F(FlushCompareTests, LowerBound_KmerQuery_End) {
  auto it = std::lower_bound(kmers.begin(), kmers.end(), *terminals.constAt(8),
                             LessThan);
  ASSERT_EQ(it, kmers.end());
}

TEST_F(FlushCompareTests, UpperBound_TerminalQuery_First) {
  auto it = std::upper_bound(terminals.constBegin(), terminals.constEnd(),
                             *kmers.begin(), LessThan);
  ASSERT_EQ(it, terminals.constAt(5));
}

TEST_F(FlushCompareTests, UpperBound_TerminalQuery_Seen) {
  auto it = std::upper_bound(terminals.constBegin(), terminals.constEnd(),
                             *kmers.at(2), LessThan);
  ASSERT_EQ(it, terminals.constAt(7));
}

TEST_F(FlushCompareTests, UpperBound_TerminalQuery_Begin) {
  std::vector<uint8_t> qry = {0b00000000, 0b00000000, 0, 0,
                              0,          0,          0, 0b00000001};
  auto it = std::upper_bound(terminals.constBegin(), terminals.constEnd(),
                             ndim::Span<uint8_t *>{qry.data(), 8}, LessThan);
  ASSERT_EQ(it, terminals.constAt(4));
}

TEST_F(FlushCompareTests, UpperBound_TerminalQuery_End) {
  std::vector<uint8_t> qry = {0b00000000, 0b01000000, 0, 0,
                              0,          0,          0, 0b00000001};
  auto it = std::upper_bound(terminals.constBegin(), terminals.constEnd(),
                             ndim::Span<uint8_t *>{qry.data(), 8}, LessThan);
  ASSERT_EQ(it, terminals.constEnd());
}

class GraphFlushingTests : public FlushCompareTests {
protected:
  GraphFlushingTests()
      : FlushCompareTests(), Bk(adjacentDifference(kmers)),
        Bt(adjacentDifference(terminals)) {}

  sdsl::int_vector<2> Bk, Bt;
};

TEST_F(GraphFlushingTests, FlushRange) {
  // Push to graph
  sdsl::int_vector<4> edges;
  sdsl::bit_vector succ;
  ArchivePayload carch;
  MetaColours cmap(7);

  auto terms = terminals.asRange();
  interleave(kmers, Bk.begin(), terms, Bt.begin(), edges, succ, carch, cmap,
             /* T */ 3);

  // Check graph structure
  sdsl::int_vector<4> Xedges = {0b1001, 0b1011, 0b1100, 0b1001, 0b1001,
                                0b1010, 0b1001, 0b1100, 0b1100, 0b1001,
                                0b0001, 0b1010, 0b0010, 0b1100};
  sdsl::bit_vector Xwplus = {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  ASSERT_THAT(edges, ContainerEq(Xedges));
  ASSERT_THAT(succ, ContainerEq(Xwplus));

  // Check colours
  std::vector<uint64_t> cols = {0, 0, 0, 0, 0, 8, 0, 0, 3, 9, 10, 0, 0, 0};
  EXPECT_THAT(carch.raw.view(), ElementsAreArray(cols));
  EXPECT_EQ(cmap.id({0, 1}), 8);
  EXPECT_EQ(cmap.id({0, 7}), 9);
  EXPECT_EQ(cmap.id({1, 5}), 10);
}
