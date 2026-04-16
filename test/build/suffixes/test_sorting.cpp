#include "maki/build/io/fasta.hpp"
#include "maki/build/kmers/buffers/kmers.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"

#include <bitset>
#include <cstddef>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <maki/maki.h>
#include <vector>

using ::testing::Eq;
using ::testing::Value;

MATCHER_P(MatrixEq, expected_wrapper, "ndim::matrix") {
  auto &expected = expected_wrapper.get();
  if (arg.size() != expected.size()) {
    *result_listener << "arg.size() != expected.size() ";
    *result_listener << arg.size() << " vs " << expected.size();
    return false;
  }
  if (arg.rows() != expected.rows()) {
    *result_listener << "arg.rows() != expected.rows() ";
    *result_listener << arg.rows() << " vs " << expected.rows();
    return false;
  }
  if (arg.cols() != expected.cols()) {
    *result_listener << "arg.cols() != expected.cols() ";
    *result_listener << arg.cols() << " vs " << expected.cols();
    return false;
  }
  for (size_t r = 0; r < expected.rows(); ++r) {
    auto arg_row = arg[r];
    auto expected_row = expected[r];
    for (size_t c = 0; c < expected.cols(); ++c) {
      if (!Value(arg_row[c], Eq(expected_row[c]))) {
        *result_listener << "element[" << r << ", " << c << "] mismatch ";
        *result_listener << "0b" << std::bitset<8>(arg_row[c]) << " vs "
                         << "0b" << std::bitset<8>(expected_row[c]);
        return false;
      }
    }
  }
  return true;
}

TEST(TwoGenomesTests, Buffer_LSDRadixSort) {
  // Pre-amble
  Colours c;
  std::size_t k = 9;
  std::vector<ChunkedDna4Genome> genomes(2);
  genomes[0] = parseFilterFNA(STRING(SMALL_SEQ), c, 1, k);
  genomes[1] = parseFilterFNA(STRING(SMALL_SEQ), c, 1, k);
  auto view = toView(genomes);

  size_t num_kmers = 20;
  KmerBuffer buffer_a(num_kmers, 1, k, k - 2), buffer_b(num_kmers, 1, k, k - 2);

  ShortSuffix soi(2, 0b1010);
  auto blocks = createSuffixPlan(view, k, 2);
  buffer_a.fill(view, blocks, soi);

  // Sort and Test
  uint8_t *a_view = buffer_a.data(), *b_view = buffer_b.data();
  lsdRadixSort(a_view, b_view, num_kmers, buffer_a.recordBytes(),
               buffer_a.keyBytes(), 1);
  auto *result = buffer_a.keyBytes() % 2 ? &buffer_b : &buffer_a;

  // check
  auto expected = ndim::Matrix({{0b01111011, 0b00001010, 0b00001001},
                                {0b01111011, 0b00001010, 0b00001001},
                                {0b01111011, 0b00001010, 0b00010001},
                                {0b01111011, 0b00001010, 0b00010001},
                                {0b10010000, 0b00011000, 0b00001011},
                                {0b10010000, 0b00011000, 0b00001011},
                                {0b10010000, 0b00011000, 0b00010011},
                                {0b10010000, 0b00011000, 0b00010011},
                                {0b00100111, 0b00110001, 0b00001000},
                                {0b00100111, 0b00110001, 0b00001000},
                                {0b00100111, 0b00110001, 0b00010000},
                                {0b00100111, 0b00110001, 0b00010000},
                                {0b01111011, 0b00110010, 0b00001000},
                                {0b01111011, 0b00110010, 0b00001000},
                                {0b01111011, 0b00110010, 0b00010000},
                                {0b01111011, 0b00110010, 0b00010000},
                                {0b01100010, 0b00111010, 0b00001001},
                                {0b01100010, 0b00111010, 0b00001001},
                                {0b01100010, 0b00111010, 0b00010001},
                                {0b01100010, 0b00111010, 0b00010001}});
  ASSERT_THAT(result->memoryview(), MatrixEq(std::cref(expected)));
}
