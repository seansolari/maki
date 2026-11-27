
#include "maki/build/kmers/buffers/kmers.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using ::testing::ElementsAreArray;
using seqan3::operator""_dna4;

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