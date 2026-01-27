
#include "maki/build/graph/build_colours.hpp"
#include "maki/core/graph/colours.hpp"
#include "gmock/gmock.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <variant>

using ::testing::ElementsAreArray;

void MakeColours(Colours &c, std::size_t n_) {
  for (std::size_t i = 0; i < n_; ++i) {
    c.getOrAssign(std::to_string(i));
  }
}

MetaColours MakeMetaColours(std::size_t n_, std::size_t c_) {
  Colours c;
  MakeColours(c, n_);
  return MetaColours(std::move(c.ids), c_);
}

std::vector<ColourVector> MakeMetaIds(colour_t n_) {
  std::vector<ColourVector> data;
  data.reserve(n_ - 1);
  for (colour_t i = 1; i < n_; ++i) {
    data.emplace_back(ColourVector({i - 1, i}));
  }
  return data;
}

void InsertMetaIds(MetaColours &cmap, std::size_t R) {
  auto cmb = MakeMetaIds(cmap.numColours());

  // first assign IDs in order
  for (const auto &c : cmb)
    cmap.insert(c);

  // multi-threaded increment
  oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<std::size_t>(0, cmb.size(), 10),
      [&](oneapi::tbb::blocked_range<std::size_t> &r) {
        std::size_t rep = 1;
        while (rep++ < R) {
          for (std::size_t i_ = r.begin(); i_ < r.end(); ++i_) {
            cmap.insert(cmb[i_]);
          }
        }
      });
}

TEST(ColourRegistry, SimpleMetaIds) {
  // construct meta colours
  std::size_t N = 100, R = 3;
  auto cmap = MakeMetaColours(N, N / 2);

  ASSERT_EQ(cmap.numMetaColours(), N+1);
  InsertMetaIds(cmap, R);
  ASSERT_EQ(cmap.numMetaColours(), 2*N+1);
  
  // convert to registry
  auto reg = toRegistry(std::move(cmap));

  // check seed names
  for (colour_t i = 1; i <= N; ++i) {
    ASSERT_EQ(reg.seed(i), std::to_string(i - 1));
  }

  // check single IDs
  for (uint64_t i = 0; i <= N; ++i) {
    auto ids = reg.colours(i);
    ASSERT_TRUE(std::holds_alternative<colour_t>(ids));
    ASSERT_EQ(std::get<colour_t>(ids), (colour_t)i);
  }

  // check meta ids
  for (uint64_t i = 1; i <= N; ++i) {
    auto ids = reg.colours(i + N);
    ASSERT_TRUE(std::holds_alternative<ColourRegistry::vector_ref>(ids));
    auto &v = std::get<ColourRegistry::vector_ref>(ids).get();
    ASSERT_THAT(v, ElementsAreArray({(colour_t)i - 1, (colour_t)i}));
  }

  // check colour occurrences
  for (colour_t i = 0; i <= N; ++i) {
    if (i == 0 || i == N) {
      ASSERT_EQ(reg.count(i), R);
    } else {
      ASSERT_EQ(reg.count(i), 2 * R);
    }
  }
}
