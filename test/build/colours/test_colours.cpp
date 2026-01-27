
#include "maki/build/graph/build_colours.hpp"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"
#include <cstddef>
#include <gtest/gtest.h>

TEST(InsertColours, RecognisesPreviousColour) {
  Colours c;

  ASSERT_EQ(c.size(), 1);
  auto id1 = c.getOrAssign("test");
  ASSERT_EQ(id1, 1);
  ASSERT_EQ(c.size(), 2);
  auto id2 = c.getOrAssign("test");
  ASSERT_EQ(id1, id2);
  ASSERT_EQ(c.size(), 2);
}

TEST(InsertColours, MultithreadedInsert) {
  // construct seed set
  const std::size_t N = 1000, R = 3;
  std::vector<std::string> seeds;
  seeds.reserve(N);
  for (std::size_t i = 0; i < N; ++i) {
    seeds.emplace_back(std::to_string(i));
  }

  // multithreaded insert seed set
  Colours c;
  colour_t globalMaxId = 0;
  std::mutex mtx;

  oneapi::tbb::parallel_for(
    oneapi::tbb::blocked_range<std::size_t>(0, N, N/100),
    [&](oneapi::tbb::blocked_range<std::size_t> &r){
      std::size_t rep = 0;
      colour_t maxId = 0;
      while (rep++ < R) {
        for (std::size_t i_ = r.begin(); i_ != r.end(); ++i_) {
          auto newId = c.getOrAssign(std::string{seeds[i_]});
          maxId = std::max(maxId, newId);
        }
      }
      {
        std::lock_guard lck{mtx};
        globalMaxId = std::max(globalMaxId, maxId);
      }
    }
  );

  ASSERT_EQ(globalMaxId, N);
}
