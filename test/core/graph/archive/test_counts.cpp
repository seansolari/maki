
#include "maki/core/graph/archive/counts.hpp"
#include "maki/core/utils/tempfile.hpp"
#include "test_common.hpp"
#include <gtest/gtest.h>

TEST(CountingVectorTests, SmallBalanced) {
  NodeTypeDistribution dist{
      .single = 0.25, .uniform = 0.25, .delta = 0.25, .expl = 0.25};

  run_counting_test(16, dist);
}

TEST(CountingVectorTests, TinySingleNode) {
  NodeTypeDistribution dist{
      .single = 1.0, .uniform = 0.0, .delta = 0.0, .expl = 0.0};

  run_counting_test(1, dist);
}

TEST(CountingVectorTests, MediumBalanced) {
  NodeTypeDistribution dist{
      .single = 0.40, .uniform = 0.30, .delta = 0.20, .expl = 0.10};

  run_counting_test(1'000, dist);
}

TEST(CountingVectorTests, AlmostAllSingle) {
  NodeTypeDistribution dist{
      .single = 0.98, .uniform = 0.01, .delta = 0.01, .expl = 0.00};

  run_counting_test(5'000, dist);
}

TEST(CountingVectorTests, AlmostAllUniform) {
  NodeTypeDistribution dist{
      .single = 0.01, .uniform = 0.97, .delta = 0.01, .expl = 0.01};

  run_counting_test(5'000, dist);
}

TEST(CountingVectorTests, AlmostAllDelta) {
  NodeTypeDistribution dist{
      .single = 0.01, .uniform = 0.01, .delta = 0.97, .expl = 0.01};

  run_counting_test(5'000, dist);
}

TEST(CountingVectorTests, AlmostAllExplicit) {
  NodeTypeDistribution dist{
      .single = 0.01, .uniform = 0.01, .delta = 0.01, .expl = 0.97};

  run_counting_test(2'000, dist);
}

TEST(CountingVectorTests, StressBalanced10000) {
  NodeTypeDistribution dist{
      .single = 0.45, .uniform = 0.30, .delta = 0.20, .expl = 0.05};

  run_counting_test(10'000, dist);
}

TEST(CountingVectorTests, StressHighlySkewed10000) {
  NodeTypeDistribution dist{
      .single = 0.90, .uniform = 0.05, .delta = 0.04, .expl = 0.01};

  run_counting_test(10'000, dist);
}

TEST(CountingVectorTests, MultipleSeeds) {
  NodeTypeDistribution dist{
      .single = 0.50, .uniform = 0.25, .delta = 0.20, .expl = 0.05};

  for (uint32_t seed = 1; seed <= 10; ++seed) {
    run_counting_test(2'000, dist, seed);
  }
}

TEST(CountingVectorTests, AppendBasic) {
  std::array<uint64_t, 4> edges = {0, 0, 0, 0};

  CountBuffer a, b;

  edges[0] = 1;
  a.insert_node({.outdegree = 1, .edge_counts = edges.data()});
  edges = {2, 2, 0, 0};
  a.insert_node({.outdegree = 2, .edge_counts = edges.data()});
  edges = {3, 4, 5, 0};
  a.insert_node({.outdegree = 3, .edge_counts = edges.data()});
  edges = {6, 1'000'000, 7, 0};
  a.insert_node({.outdegree = 3, .edge_counts = edges.data()});

  edges[0] = 8;
  b.insert_node({.outdegree = 1, .edge_counts = edges.data()});
  edges = {9, 9, 0, 0};
  b.insert_node({.outdegree = 2, .edge_counts = edges.data()});
  edges = {10, 11, 12, 0};
  b.insert_node({.outdegree = 3, .edge_counts = edges.data()});
  edges = {13, 1'000'001, 14, 0};
  b.insert_node({.outdegree = 3, .edge_counts = edges.data()});

  CompressedCountBuffer xa(a), xb(b);

  EXPECT_EQ(xa.node_count(), 4);
  EXPECT_EQ(xa.edge_count(0, 0), 1);
  EXPECT_EQ(xa.edge_count(1, 0), 2);
  EXPECT_EQ(xa.edge_count(1, 1), 2);
  EXPECT_EQ(xa.edge_count(2, 0), 3);
  EXPECT_EQ(xa.edge_count(2, 1), 4);
  EXPECT_EQ(xa.edge_count(2, 2), 5);
  EXPECT_EQ(xa.edge_count(3, 0), 6);
  EXPECT_EQ(xa.edge_count(3, 1), 1'000'000);
  EXPECT_EQ(xa.edge_count(3, 2), 7);

  EXPECT_EQ(xb.node_count(), 4);
  EXPECT_EQ(xb.edge_count(0, 0), 8);
  EXPECT_EQ(xb.edge_count(1, 0), 9);
  EXPECT_EQ(xb.edge_count(1, 1), 9);
  EXPECT_EQ(xb.edge_count(2, 0), 10);
  EXPECT_EQ(xb.edge_count(2, 1), 11);
  EXPECT_EQ(xb.edge_count(2, 2), 12);
  EXPECT_EQ(xb.edge_count(3, 0), 13);
  EXPECT_EQ(xb.edge_count(3, 1), 1'000'001);
  EXPECT_EQ(xb.edge_count(3, 2), 14);

  a.append(b);

  CompressedCountBuffer xab(a);

  EXPECT_EQ(xab.node_count(), 8);
  EXPECT_EQ(xab.edge_count(0, 0), 1);
  EXPECT_EQ(xab.edge_count(1, 0), 2);
  EXPECT_EQ(xab.edge_count(1, 1), 2);
  EXPECT_EQ(xab.edge_count(2, 0), 3);
  EXPECT_EQ(xab.edge_count(2, 1), 4);
  EXPECT_EQ(xab.edge_count(2, 2), 5);
  EXPECT_EQ(xab.edge_count(3, 0), 6);
  EXPECT_EQ(xab.edge_count(3, 1), 1'000'000);
  EXPECT_EQ(xab.edge_count(3, 2), 7);
  EXPECT_EQ(xab.edge_count(4, 0), 8);
  EXPECT_EQ(xab.edge_count(5, 0), 9);
  EXPECT_EQ(xab.edge_count(5, 1), 9);
  EXPECT_EQ(xab.edge_count(6, 0), 10);
  EXPECT_EQ(xab.edge_count(6, 1), 11);
  EXPECT_EQ(xab.edge_count(6, 2), 12);
  EXPECT_EQ(xab.edge_count(7, 0), 13);
  EXPECT_EQ(xab.edge_count(7, 1), 1'000'001);
  EXPECT_EQ(xab.edge_count(7, 2), 14);
}

TEST(CountingVectorTests, AppendMultipleSeeds) {
  NodeTypeDistribution dist{
      .single = 0.25, .uniform = 0.25, .delta = 0.25, .expl = 0.25};

  for (uint32_t seed = 1; seed <= 10; ++seed) {
    run_append_test(2'000, 1'000, dist, seed);
  }
}

class SerializeTests : public testing::Test {
protected:
  SerializeTests() : tmp(tempio::create_temporary_directory()) {}
  ~SerializeTests() { fs::remove_all(tmp); }
  fs::path tmp;

  static inline void save(const std::string &file, const CompressedCountBuffer &data) {
    std::ofstream os(file, std::ios::binary);
    cereal::BinaryOutputArchive archive( os );
    archive(data);
  }

  static inline void load(const std::string &file, CompressedCountBuffer &data) {
    std::ifstream os(file, std::ios::binary);
    cereal::BinaryInputArchive archive( os );
    archive(data);
  }
};

TEST_F(SerializeTests, SaveAndLoadRandom) {
  NodeTypeDistribution dist{
      .single = 0.25, .uniform = 0.25, .delta = 0.25, .expl = 0.25};
  constexpr std::size_t size = 1'000;

  CountBuffer builder;
  builder.reserve(size);
  auto types = generate_node_types(size, dist, 1);
  DummyOracle oracle = build_oracle_and_builder(builder, types, 2);
  CompressedCountBuffer compressed(builder);
  builder.clear();
  
  fs::path bin = tmp / "SaveAndLoadRandom.bin";
  save(bin, compressed);

  CompressedCountBuffer result;
  load(bin, result);

  check_counts_equal(oracle, result);
}
