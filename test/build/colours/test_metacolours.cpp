
#include "maki/build/graph/build_colours.hpp"
#include "maki/core/graph/colours.hpp"
#include <cereal/external/rapidjson/reader.h>
#include <gtest/gtest.h>

void MakeColours(Colours &c, std::size_t n_) {
  for (std::size_t i = 0; i < n_; ++i) {
    c.getOrAssign(std::to_string(i));
  }
}

MetaColours MakeMetaColours(std::size_t n_) {
  Colours c;
  MakeColours(c, n_);
  return MetaColours(std::move(c.ids));
}

TEST(MetaColoursInsert, RecognisesPreviousColour) {
  auto cmap = MakeMetaColours(100);

  ASSERT_EQ(cmap.numMetaColours(), 101);
  auto id1 = cmap.id({1, 2, 3});
  ASSERT_EQ(id1, 101);
  ASSERT_EQ(cmap.numMetaColours(), 102);
  auto id2 = cmap.id({1, 2, 3});
  ASSERT_EQ(id1, id2);
  ASSERT_EQ(cmap.numMetaColours(), 102);
}

std::vector<ColourVector> MakeMetaIds(std::size_t n_) {
  std::vector<ColourVector> data;
  data.reserve(n_-1);
  for (std::size_t i = 1; i < n_; ++i) {
    data.emplace_back(ColourVector{(colour_t)i-1, (colour_t)i});
  }
  return data;
}

TEST(MetaColoursInsert, TrackMultiInsertSingleValues) {
  std::size_t N = 100, R = 3;
  auto cmap = MakeMetaColours(N);

  // multithreaded insert seed set
  Colours c;
  uint64_t globalMaxId = 0;
  std::mutex mtx;

  oneapi::tbb::parallel_for(
    oneapi::tbb::blocked_range<std::size_t>(1, N+1, N/100),
    [&](oneapi::tbb::blocked_range<std::size_t> &r){
      std::size_t rep = 0;
      uint64_t maxId = 0;
      while (rep++ < R) {
        for (std::size_t i_ = r.begin(); i_ != r.end(); ++i_) {
          auto newId = cmap.id({(colour_t)i_});
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

TEST(MetaColoursInsert, TrackMultiInsertVectors) {
  std::size_t N = 100, M = N-1, R = 3;
  auto cmap = MakeMetaColours(N);
  auto cmb = MakeMetaIds(N);

  // multithreaded insert seed set
  Colours c;
  uint64_t globalMaxId = 0;
  std::mutex mtx;

  oneapi::tbb::parallel_for(
    oneapi::tbb::blocked_range<std::size_t>(0, M, M/100),
    [&](oneapi::tbb::blocked_range<std::size_t> &r){
      std::size_t rep = 0;
      uint64_t maxId = 0;
      while (rep++ < R) {
        for (std::size_t i_ = r.begin(); i_ != r.end(); ++i_) {
          auto newId = cmap.id(cmb[i_]);
          maxId = std::max(maxId, newId);
        }
      }
      {
        std::lock_guard lck{mtx};
        globalMaxId = std::max(globalMaxId, maxId);
      }
    }
  );

  ASSERT_EQ(globalMaxId, N+M);
}
