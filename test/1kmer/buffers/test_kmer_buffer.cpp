
#include "maki/build/kmers/buffers/kmers.hpp"
#include "maki/core/seq/io.hpp"
#include "seqan3/alphabet/nucleotide/dna4.hpp"
#include <cstdint>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using ::testing::ElementsAreArray;
using seqan3::operator""_dna4;

static std::string kmerString(auto const &it, uint8_t k) {
  return toString(it.readKmer(k));
}

TEST(KmerEncodingTests, WriteKmerBits) {
  Dna4Sequence original{"ACGTACGTAC"_dna4};

  std::array<uint8_t, 3> buffer{};
  writeMerTo(original.begin(), buffer.data(), 10);

  auto decoded = readMerFrom(buffer.data(), 10);

  EXPECT_EQ(original, decoded);
}

TEST(KmerBufferInsertFragment, NonTerminalFragment) {
  constexpr uint8_t k = 3;

  KmerBuffer buffer(10, sizeof(uint32_t), k, k);

  std::vector<seqan3::dna4> seq = "ATACGT"_dna4;

  SequenceFragment fragment(std::move(seq), 17, false);

  auto it = buffer.begin();

  it = buffer.insertFragment(it, fragment);

  auto r = buffer.begin();

  EXPECT_EQ(kmerString(r, k), "ATA");
  EXPECT_EQ(r.readEdge(), seqan3::dna4{}.assign_char('C').to_rank());
  EXPECT_EQ(r.readValue().colour(), 17u);

  ++r;

  EXPECT_EQ(kmerString(r, k), "TAC");
  EXPECT_EQ(r.readEdge(), seqan3::dna4{}.assign_char('G').to_rank());
  EXPECT_EQ(r.readValue().colour(), 17u);

  ++r;

  EXPECT_EQ(kmerString(r, k), "ACG");
  EXPECT_EQ(r.readEdge(), seqan3::dna4{}.assign_char('T').to_rank());
  EXPECT_EQ(r.readValue().colour(), 17u);
}

TEST(KmerBufferInsertFragment, TerminalFragmentAddsTerminalEdge) {
  constexpr uint8_t k = 3;

  KmerBuffer buffer(10, sizeof(uint32_t), k, k);

  std::vector<seqan3::dna4> seq = "ATACGT"_dna4;

  SequenceFragment fragment(std::move(seq), 5, true);

  buffer.insertFragment(buffer.begin(), fragment);

  auto r = buffer.begin();

  EXPECT_EQ(kmerString(r, k), "ATA");
  ++r;

  EXPECT_EQ(kmerString(r, k), "TAC");
  ++r;

  EXPECT_EQ(kmerString(r, k), "ACG");
  ++r;

  EXPECT_EQ(kmerString(r, k), "CGT");

  EXPECT_EQ(r.readEdge(), KmerBuffer::terminalEdge);
  EXPECT_EQ(r.readValue().colour(), 0u);
}

TEST(KmerBufferInsertFragment, PreservesInputOrder) {
  constexpr uint8_t k = 3;

  KmerBuffer buffer(16, sizeof(uint32_t), k, k);

  std::vector<seqan3::dna4> seq = "AAAACCCC"_dna4;

  SequenceFragment fragment(std::move(seq), 1, false);

  buffer.insertFragment(buffer.begin(), fragment);

  std::vector<std::string> expected{"AAA", "AAA", "AAC", "ACC", "CCC"};

  auto r = buffer.begin();

  for (auto const &kmer : expected) {
    EXPECT_EQ(kmerString(r, k), kmer);
    ++r;
  }
}

TEST(KmerBufferInsertFragment, SuffixFilterSingleMatch) {
  constexpr uint8_t k = 3;

  KmerBuffer buffer(10, sizeof(uint32_t), k, k - 1);

  std::vector<seqan3::dna4> seq = "ATACGT"_dna4;

  SequenceFragment fragment(std::move(seq), 11, false);

  ShortSuffix suffix("G"_dna4);

  buffer.insertFragment(buffer.begin(), fragment, suffix);

  auto r = buffer.begin();

  EXPECT_EQ(kmerString(r, k - 1), "AC");

  EXPECT_EQ(r.readEdge(), seqan3::dna4{}.assign_char('T').to_rank());

  EXPECT_EQ(r.readValue().colour(), 11u);
}

TEST(KmerBufferInsertFragment, SuffixFilterMultipleMatches) {
  constexpr uint8_t k = 3;

  KmerBuffer buffer(10, sizeof(uint32_t), k, 1);

  std::vector<seqan3::dna4> seq = "ACGACGAC"_dna4;

  SequenceFragment fragment(std::move(seq), 9, false);

  ShortSuffix suffix("CG"_dna4);

  buffer.insertFragment(buffer.begin(), fragment, suffix);

  auto r = buffer.begin();

  EXPECT_EQ(kmerString(r, 1), "A");
  EXPECT_EQ(r.readEdge(), seqan3::dna4{}.assign_char('A').to_rank());

  ++r;

  EXPECT_EQ(kmerString(r, 1), "A");
  EXPECT_EQ(r.readEdge(), seqan3::dna4{}.assign_char('A').to_rank());
}

TEST(KmerBufferInsertFragment, SuffixFilterPreservesTerminalEdge) {
  constexpr uint8_t k = 3;

  KmerBuffer buffer(10, sizeof(uint32_t), k, 1);

  std::vector<seqan3::dna4> seq = "ACGT"_dna4;

  SequenceFragment fragment(std::move(seq), 3, true);

  ShortSuffix suffix("GT"_dna4);

  buffer.insertFragment(buffer.begin(), fragment, suffix);

  auto r = buffer.begin();

  EXPECT_EQ(kmerString(r, 1), "C");
  EXPECT_EQ(r.readEdge(), TerminalEdge);
  EXPECT_EQ(r.readValue().colour(), 0u);
}

TEST(KmerBufferInsertFragment, SuffixFilterAcrossByteBoundary) {
  constexpr uint8_t k = 9;

  KmerBuffer buffer(10, sizeof(uint32_t), k, 6);

  std::vector<seqan3::dna4> seq = "ACGTACGTAC"_dna4;

  SequenceFragment fragment(std::move(seq), 42, false);

  ShortSuffix suffix("GTA"_dna4);

  buffer.insertFragment(buffer.begin(), fragment, suffix);

  auto r = buffer.begin();

  EXPECT_EQ(kmerString(r, 6), "ACGTAC");

  EXPECT_EQ(r.readEdge(), seqan3::dna4{}.assign_char('C').to_rank());

  EXPECT_EQ(r.readValue().colour(), 42u);
}

// compare k-mers
// ----------------------------------------------------------------------------

/**
 * Test buffer design
 * ------------------
 *
 * This structure captures an edge case at index 1.
 *
 *  0: CCC C 0
 *  1: CCC C 0 <-- DUPLICATE
 *  2: CCC A 0
 *  3: CCA C 0
 *  4: CAC A 0
 *  5: ACA $ 0
 *  6: ACA C 1
 *  7: CAC C 1
 *  8: ACC $ 1
 *  9: GCC $ 2
 * 10: TCC $ 3
 * 11: TTT T 4
 * 12: TTT T 4 <-- DUPLICATE
 * 13: TTT $ 4
 *
 */
class KmerComparisonTests : public testing::Test {
protected:
  KmerComparisonTests() : k(3), buffer(14, 1, k), SmallDiff(1), LongDiff(3) {
    // insert sequences into buffer
    size_t i = 0;
    auto p = buffer.begin();

    std::vector sequences = {Dna4Sequence("CCCCCACA"_dna4),
                             Dna4Sequence("ACACC"_dna4),
                             Dna4Sequence("GCC"_dna4), Dna4Sequence("TCC"_dna4),
                             Dna4Sequence("TTTTT"_dna4)};

    for (auto &seq : sequences) {
      p = buffer.insertFragment(p, SequenceFragment{std::move(seq), i++, true});
    }
  }

  uint8_t k;
  KmerBuffer buffer;
  KmerDiff SmallDiff, LongDiff;
};

TEST_F(KmerComparisonTests, Short_IsZero) {
  auto r1 = buffer[0], r2 = buffer[1];

  ASSERT_EQ(SmallDiff(r1.data(), r2.data()), IS_0);
}

TEST_F(KmerComparisonTests, Short_BetweenZeroAndK) {
  auto r1 = buffer[2], r2 = buffer[3];

  ASSERT_EQ(SmallDiff(r1.data(), r2.data()), BW_0_K);
}

TEST_F(KmerComparisonTests, Short_IsK) {
  auto r1 = buffer[8], r2 = buffer[9];

  ASSERT_EQ(SmallDiff(r1.data(), r2.data()), IS_K);
}

TEST_F(KmerComparisonTests, Long_IsZero) {
  std::vector<uint8_t> r1 = {0b10101010, 0b10101010, 0b00001010},
                       r2 = {0b10101010, 0b10101010, 0b00001010};

  ASSERT_EQ(LongDiff(r1.data(), r2.data()), IS_0);
}

TEST_F(KmerComparisonTests, Long_BetweenZeroAndK) {
  std::vector<uint8_t> r1 = {0b10101010, 0b10101010, 0b00001111},
                       r2 = {0b10101010, 0b10111010, 0b00001111};

  ASSERT_EQ(LongDiff(r1.data(), r2.data()), BW_0_K);
}

TEST_F(KmerComparisonTests, Long_BetweenZeroAndK_LSB) {
  std::vector<uint8_t> r1 = {0b10101010, 0b10101010, 0b00001111},
                       r2 = {0b10101110, 0b10101010, 0b00001111};

  ASSERT_EQ(LongDiff(r1.data(), r2.data()), BW_0_K);
}

TEST_F(KmerComparisonTests, Long_IsK) {
  std::vector<uint8_t> r1 = {0b10101010, 0b10101010, 0b00001111},
                       r2 = {0b10101000, 0b10101010, 0b00001111};

  ASSERT_EQ(LongDiff(r1.data(), r2.data()), IS_K);
}

TEST_F(KmerComparisonTests, AdjacentDifference) {
  auto b = adjacentDifference(buffer);
  ASSERT_THAT(
      b,
      ElementsAreArray({BW_0_K, // dummy value
                        IS_0, IS_0, BW_0_K, BW_0_K, BW_0_K, IS_0, BW_0_K,
                        BW_0_K, IS_K, IS_K, BW_0_K, IS_0, IS_0}));
}
