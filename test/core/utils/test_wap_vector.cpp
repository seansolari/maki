
#include "maki/core/utils/algo.hpp"
#include "maki/core/utils/tempfile.hpp"
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

class SerializeTests : public testing::Test {
protected:
  SerializeTests() : tmp(tempio::create_temporary_directory()) {}
  ~SerializeTests() { fs::remove_all(tmp); }
  fs::path tmp;

  static inline void save(const std::string &file, const ro_wap_vector &data) {
    std::ofstream os(file, std::ios::binary);
    cereal::BinaryOutputArchive archive( os );
    archive(data);
  }

  static inline void load(const std::string &file, ro_wap_vector &data) {
    std::ifstream os(file, std::ios::binary);
    cereal::BinaryInputArchive archive( os );
    archive(data);
  }
};

TEST_F(SerializeTests, SaveAndLoadRandom) {
  std::mt19937_64 rng(42);
  constexpr std::size_t total = 1'000;
  auto truth = random_vector_equal_bits_bulk(rng, total);

  ro_wap_vector data = MakeRoWapVector(truth);
  auto bin = tmp / "SaveAndLoadRandom.dat";
  save(bin, data);

  ro_wap_vector result;
  load(bin, result);

  EXPECT_EQ(result.size(), data.size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    EXPECT_EQ(result[i], data[i]);
  }
}
