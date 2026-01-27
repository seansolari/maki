
#include "maki/build/graph/wm/wm_construct.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/utils/tempfile.hpp"
#include <filesystem>
#include <gtest/gtest.h>
#include <random>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/io.hpp>

inline std::string MakeTempPath(const char *ext) {
  return tempio::create_temporary_file(std::filesystem::temp_directory_path(),
                                       ext);
}

template <typename Int = uint64_t>
inline std::vector<Int> GenVector(size_t n, Int mask) {
  std::vector<Int> v(n);
  std::mt19937_64 rng(12345);
  for (size_t i = 0; i < n; ++i)
    v[i] = static_cast<Int>(rng()) & mask;
  return v;
}

class WaveletMatrix : public testing::Test {
protected:
  WaveletMatrix() : file(MakeTempPath(".char")), arr(MakeTempPath(".b4")) {}
  ~WaveletMatrix() {
    std::filesystem::remove(file);
    std::filesystem::remove(arr);
  }
  std::string file, arr;

  void MakeWaveletMatrix(wavelet_matrix &wm, std::vector<uint8_t> &vals) const {
    // create vector on disk
    {
      sdsl::int_vector_buffer<4> ivb(file, std::ios::out);
      for (auto &v : vals)
        ivb.push_back(v);
      ivb.close();
    }

    // construct wavelet matrix
    initW(file, arr, std::filesystem::temp_directory_path());

    // load wavelet matrix
    sdsl::load_from_file(wm, arr);
  }
};

TEST_F(WaveletMatrix, SmallArray) {
  std::vector<uint8_t> vals = {0b0001, 0b1001, 0b0001, 0b1001,
                               0b0001, 0b1001, 0b0001, 0b1001};

  wavelet_matrix wm;
  MakeWaveletMatrix(wm, vals);

  // check values
  ASSERT_EQ(wm.size(), 8);
  ASSERT_LE(wm.sigma, 16);
  for (std::size_t i = 0; i < 8; ++i)
    ASSERT_EQ(wm[i], vals[i]);

  // check rank select
  ASSERT_EQ(wm.rank(vals.size(), 0b1001), 4);
  ASSERT_EQ(wm.rank(3, 0b1001), 1);
  ASSERT_EQ(wm.rank(vals.size(), 0b0001), 4);
  ASSERT_EQ(wm.select(1, 0b1001), 1);
  ASSERT_EQ(wm.select(4, 0b0001), 6);
}

TEST_F(WaveletMatrix, RandomValues) {
  std::size_t N = 1000;
  std::vector<uint8_t> vals = GenVector<uint8_t>(N, (1ull << 4) - 1);

  wavelet_matrix wm;
  MakeWaveletMatrix(wm, vals);

  // check values
  ASSERT_EQ(wm.size(), N);
  ASSERT_LE(wm.sigma, 16);
  for (std::size_t i = 0; i < N; ++i)
    ASSERT_EQ(wm[i], vals[i]);
}
