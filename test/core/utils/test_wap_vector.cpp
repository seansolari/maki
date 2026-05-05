
#include "maki/core/utils/algo.hpp"
#include "maki/core/utils/wap_vector.hpp"
#include "test_common.hpp"

TEST(VectorAppend, AppendUint8) {
  std::vector<uint8_t> A = {1, 2, 3}, B = {4, 5, 6}, AB = {1, 2, 3, 4, 5, 6},
                       res = A;
  vector_append(res, B);
  ASSERT_EQ(res, AB);
}

TEST(VectorAppend, AppendUint64Add) {
  std::vector<uint64_t> A = {1, 2, 3}, B = {4, 5, 6}, AB = {1, 2, 3, 9, 10, 11},
                        res = A;
  vector_append_add(res, B, 5);
  ASSERT_EQ(res, AB);
}

TEST(WidthAdaptivePacking, Construct100Even) {
  std::mt19937_64 rng(42);
  constexpr std::size_t total = 100;
  auto truth = random_vector_equal_bits_bulk(rng, total);
  TestAppend(truth);
}

TEST(WidthAdaptivePacking, Construct100Unbalanced) {
  std::mt19937_64 rng(42);
  constexpr std::size_t total = 100;
  auto truth = random_vector_weighted_bits_bulk(rng, total, 4.0);
  TestAppend(truth);
}

TEST(WidthAdaptivePacking, Construct1e6Even) {
  std::mt19937_64 rng(42);
  constexpr std::size_t total = 1'000'000;
  auto truth = random_vector_equal_bits_bulk(rng, total);
  TestAppend(truth);
}

TEST(WidthAdaptivePacking, Construct1e6Unbalanced) {
  std::mt19937_64 rng(42);
  constexpr std::size_t total = 1'000'000;
  auto truth = random_vector_weighted_bits_bulk(rng, total, 4.0);
  TestAppend(truth);
}

TEST(WidthAdaptivePacking, Clear) {
  std::mt19937_64 rng(42);
  constexpr std::size_t total = 1'000;
  auto truth = random_vector_equal_bits_bulk(rng, total);

  wo_wap_vector interm(truth);
  EXPECT_EQ(interm.size(), total);
  interm.clear();
  EXPECT_EQ(interm.size(), 0);
}

TEST(WidthAdaptivePacking, Append100Even) {
  std::mt19937_64 rng(42);
  constexpr std::size_t totalA = 100, totalB = 1'000;
  TestAppendRaw(random_vector_equal_bits_bulk(rng, totalA),
                random_vector_equal_bits_bulk(rng, totalB));
}

TEST(WidthAdaptivePacking, Append100Unbalanced) {
  std::mt19937_64 rng(42);
  constexpr std::size_t totalA = 100, totalB = 1'000;
  TestAppendRaw(random_vector_weighted_bits_bulk(rng, totalA, 4.0),
                random_vector_weighted_bits_bulk(rng, totalB, 4.0));
}

TEST(WidthAdaptivePacking, CommutativeAppend100Even) {
  std::mt19937_64 rng(42);
  constexpr std::size_t totalA = 1'000, totalB = 10'000;
  TestAppendCommutative(random_vector_equal_bits_bulk(rng, totalA),
                        random_vector_equal_bits_bulk(rng, totalB));
}

TEST(WidthAdaptivePacking, CommutativeAppend100Unbalanced) {
  std::mt19937_64 rng(42);
  constexpr std::size_t totalA = 1'000, totalB = 10'000;
  TestAppendCommutative(random_vector_weighted_bits_bulk(rng, totalA, 4.0),
                        random_vector_weighted_bits_bulk(rng, totalB, 4.0));
}
