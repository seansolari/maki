#pragma once
#include "maki/core/utils/logging.hpp"
#include "maki/core/utils/wap_vector.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <random>
#include <vector>

template <typename URNG>
void sample_n_bitmasked(URNG &rng, unsigned bits, std::size_t n,
                        std::vector<uint64_t> &out) {
  auto dist = (bits >= 64) ? std::uniform_int_distribution<uint64_t>()
                           : std::uniform_int_distribution<uint64_t>(
                                 0, (uint64_t{1} << bits) - 1);
  for (std::size_t i = 0; i < n; ++i) {
    out.push_back(dist(rng));
  }
}

template <typename URNG>
std::array<std::size_t, 4>
weights_to_counts(URNG &rng, const std::array<double, 4> &weights,
                  std::size_t total) {
  std::discrete_distribution<std::size_t> dist(weights.begin(), weights.end());
  std::array<std::size_t, 4> counts{0, 0, 0, 0};
  for (std::size_t i = 0; i < total; ++i) {
    ++counts[dist(rng)];
  }
  return counts;
}

template <typename URNG>
std::array<std::size_t, 4> sample_equal_category_counts(URNG &rng,
                                                        std::size_t total) {
  std::array<double, 4> probs{1.0, 1.0, 1.0, 1.0};
  return weights_to_counts(rng, probs, total);
}

template <typename URNG>
std::array<std::size_t, 4>
sample_weighted_category_counts(URNG &rng, std::size_t total, double X) {
  std::array<double, 4> weights{
      X * X * X, // 8-bit
      X * X,     // 16-bit
      X,         // 32-bit
      1.0        // 64-bit
  };
  return weights_to_counts(rng, weights, total);
}

template <typename URNG>
std::vector<uint64_t> random_vector_equal_bits_bulk(URNG &rng,
                                                    std::size_t total) {
  static constexpr std::array<unsigned, 4> bit_sizes{8, 16, 32, 64};

  auto counts = sample_equal_category_counts(rng, total);

  std::vector<uint64_t> out;
  out.reserve(total);

  for (std::size_t i = 0; i < bit_sizes.size(); ++i) {
    sample_n_bitmasked(rng, bit_sizes[i], counts[i], out);
  }

  // Optional but recommended
  std::shuffle(out.begin(), out.end(), rng);
  return out;
}

template <typename URNG>
std::vector<uint64_t>
random_vector_weighted_bits_bulk(URNG &rng, std::size_t total, double X) {
  static constexpr std::array<unsigned, 4> bit_sizes{8, 16, 32, 64};

  auto counts = sample_weighted_category_counts(rng, total, X);

  std::vector<uint64_t> out;
  out.reserve(total);

  for (std::size_t i = 0; i < bit_sizes.size(); ++i) {
    sample_n_bitmasked(rng, bit_sizes[i], counts[i], out);
  }

  // Prevent category ordering artifacts
  std::shuffle(out.begin(), out.end(), rng);
  return out;
}

inline static void TestAppend(const std::vector<uint64_t> &truth) {
  // construct vector
  wo_wap_vector interm;
  for (const auto &x : truth)
    interm.push_back(x);
  EXPECT_EQ(interm.size(), truth.size());

  // compress vector
  ro_wap_vector result(std::move(interm));
  EXPECT_EQ(result.size(), truth.size());
  for (std::size_t i = 0; i < truth.size(); ++i) {
    EXPECT_EQ(result[i], truth[i]);
  }

  LOG_INFO() << "original vector=" << truth.size() * sizeof(uint64_t)
             << " bytes, compressed vector=" << detail::size_in_bytes(result)
             << " bytes";
}

inline static void TestAppendRaw(const std::vector<uint64_t> &a, const std::vector<uint64_t> &b) {
  wo_wap_vector wapA(a), wapB(b);
  wapA.append(wapB);
  EXPECT_EQ(wapA.size(), a.size() + b.size());

  ro_wap_vector result(wapA);
  auto truth = a;
  truth.insert(truth.end(), b.begin(), b.end());
  EXPECT_EQ(result.size(), truth.size());
  for (std::size_t i = 0; i < truth.size(); ++i)
    EXPECT_EQ(result[i], truth[i]);
}

inline static void TestAppendCommutative(const std::vector<uint64_t> &a, const std::vector<uint64_t> &b) {
  auto truth = a;
  truth.insert(truth.end(), b.begin(), b.end());

  wo_wap_vector AthenB(a), wapB(b);
  AthenB.append(wapB);
  EXPECT_EQ(AthenB.size(), truth.size());

  wo_wap_vector AB(truth);
  EXPECT_EQ(AB.size(), truth.size());

  ro_wap_vector first(std::move(AthenB)), second(std::move(AB));
  EXPECT_EQ(first.size(), truth.size());
  EXPECT_EQ(second.size(), truth.size());

  for (std::size_t i = 0; i < truth.size(); ++i) {
    EXPECT_EQ(first[i], truth[i]);
    EXPECT_EQ(second[i], truth[i]);
  }

  LOG_INFO() << "original vector=" << truth.size() * sizeof(uint64_t)
             << " bytes, first compressed vector=" << detail::size_in_bytes(first)
             << " bytes, second compressed vector=" << detail::size_in_bytes(second)
             << " bytes";
}
