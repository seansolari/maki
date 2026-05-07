
#include "maki/build/io/fasta.hpp"
#include "maki/build/kmers/buffers/kmers.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/core/seq/io.hpp"
#include "test_common.hpp"
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using ::testing::ElementsAreArray;
using namespace seqan3::literals;

class BasicInsert : public testing::Test {
protected:
  BasicInsert() : _seq1("ACGTA"_dna4), _seq2("TGCAT"_dna4) {}

  Dna4Sequence _seq1, _seq2;
};

TEST_F(BasicInsert, Init) {
  KmerBuffer buffer(6, 1, 3);
  for (auto row : buffer) {
    EXPECT_TRUE(row == 0);
  }
}

TEST_F(BasicInsert, BulkInsert) {
  KmerBuffer buffer(6, 1, 3);

  auto p = buffer.begin();
  p = buffer.insert(p, _seq1.cbegin(), _seq1.cend(), 0, true);
  p = buffer.insert(p, _seq2.cbegin(), _seq2.cend(), 1, true);

  EXPECT_TRUE(p == buffer.end());

  auto expected = ndim::Matrix({{0b00100100, 0b00000011},
                                {0b00111001, 0b00000000},
                                {0b00001110, 0b00000111},
                                {0b00011011, 0b00001000},
                                {0b00000110, 0b00001011},
                                {0b00110001, 0b00000111}});
  ASSERT_THAT(buffer.memoryview(), MatrixEq(std::cref(expected)));
}

TEST(BasicInsertSuffix, SingleSequence) {
  Dna4Sequence sequence("ACGTACGTACGTACGT"_dna4);
  std::size_t k = 5, s = 2;
  KmerBuffer buffer(3, 1, k, k - s);
  buffer.insert(buffer.begin(), sequence.cbegin(), sequence.cend(), 1, true,
                ShortSuffix{s, 0b0011});

  ASSERT_THAT(buffer[0], ElementsAreArray({0b00100100, 0b00001001}));
  ASSERT_THAT(buffer[1], ElementsAreArray({0b00100100, 0b00001001}));
  ASSERT_THAT(buffer[2], ElementsAreArray({0b00100100, 0b00001001}));
}

TEST(BasicInsertSuffix, MultipleSequences) {
  Colours c;
  std::size_t k = 9, s = 2;
  std::vector<ChunkedDna4Genome> genomes(2);
  genomes[0] = parseFilterFNA(STRING(SMALL_SEQ), c, 1, k);
  genomes[1] = parseFilterFNA(STRING(SMALL_SEQ), c, 1, k);
  auto view = toView(genomes);

  // insert specific suffix
  ShortSuffix sfx(2, 0b1010);
  auto blocks = createSuffixPlan(view, k, s, k - s);

  KmerBuffer buffer(blocks.back()[sfx], 1, k, k - s);
  buffer.fill(view, blocks, sfx);

  // check
  auto expected = ndim::Matrix({{0b01111011, 0b00001010, 0b00001001},
                                {0b00100111, 0b00110001, 0b00001000},
                                {0b10010000, 0b00011000, 0b00001011},
                                {0b01100010, 0b00111010, 0b00001001},
                                {0b01111011, 0b00110010, 0b00001000},
                                {0b01111011, 0b00001010, 0b00010001},
                                {0b00100111, 0b00110001, 0b00010000},
                                {0b10010000, 0b00011000, 0b00010011},
                                {0b01100010, 0b00111010, 0b00010001},
                                {0b01111011, 0b00110010, 0b00010000},
                                {0b01111011, 0b00001010, 0b00001001},
                                {0b00100111, 0b00110001, 0b00001000},
                                {0b10010000, 0b00011000, 0b00001011},
                                {0b01100010, 0b00111010, 0b00001001},
                                {0b01111011, 0b00110010, 0b00001000},
                                {0b01111011, 0b00001010, 0b00010001},
                                {0b00100111, 0b00110001, 0b00010000},
                                {0b10010000, 0b00011000, 0b00010011},
                                {0b01100010, 0b00111010, 0b00010001},
                                {0b01111011, 0b00110010, 0b00010000}});
  ASSERT_THAT(buffer.memoryview(), MatrixEq(std::cref(expected)));
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

    for (const Dna4Sequence &s :
         {Dna4Sequence("CCCCCACA"_dna4), Dna4Sequence("ACACC"_dna4),
          Dna4Sequence("GCC"_dna4), Dna4Sequence("TCC"_dna4),
          Dna4Sequence("TTTTT"_dna4)}) {
      p = buffer.insert(p, s.cbegin(), s.cend(), i++, true);
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
