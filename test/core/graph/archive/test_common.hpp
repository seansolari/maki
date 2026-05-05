
#pragma once
#include "maki/core/graph/archive/counts.hpp"
#include "maki/core/utils/logging.hpp"
#include <gtest/gtest.h>
#include <random>

struct DummyNode {
  NodeCountClass type;
  std::vector<uint64_t> edge_counts;

  std::size_t outdegree() const { return edge_counts.size(); }

  uint64_t edge_count(uint32_t edge_rank) const {
    return edge_counts.at(edge_rank);
  }
};

struct NodeTypeDistribution {
  double single = 0.5;
  double uniform = 0.25;
  double delta = 0.20;
  double expl = 0.05;

  void normalize() {
    double tot = single + uniform + delta + expl;
    single /= tot;
    uniform /= tot;
    delta /= tot;
    expl /= tot;
  }
};

static inline std::vector<NodeCountClass>
generate_node_types(size_t n, NodeTypeDistribution &dist, uint32_t seed = 42) {
  dist.normalize();

  std::vector<NodeCountClass> types;
  types.reserve(n);

  std::mt19937 rng(seed);
  std::discrete_distribution<int> dd(
      {dist.single, dist.uniform, dist.delta, dist.expl});

  for (size_t i = 0; i < n; ++i) {
    switch (dd(rng)) {
    case 0:
      types.push_back(NodeCountClass::SINGLE);
      break;
    case 1:
      types.push_back(NodeCountClass::UNIFORM);
      break;
    case 2:
      types.push_back(NodeCountClass::DELTA);
      break;
    case 3:
      types.push_back(NodeCountClass::EXPLICIT);
      break;
    }
  }

  return types;
}

template <typename URNG>
static inline DummyNode make_dummy_node(NodeCountClass type, URNG &rng) {
  std::uniform_int_distribution<uint32_t> base_dist(1, 100);
  std::uniform_int_distribution<uint32_t> degree_dist(1, 4);
  std::uniform_int_distribution<int> delta_dist(-3, 3);

  DummyNode node;
  node.type = type;

  switch (type) {
  case NodeCountClass::NONE:
    break;

  case NodeCountClass::SINGLE: {
    node.edge_counts = {base_dist(rng)};
    break;
  }

  case NodeCountClass::UNIFORM: {
    uint32_t d = degree_dist(rng) + 1;
    uint32_t v = base_dist(rng);
    node.edge_counts.assign(d, v);
    break;
  }

  case NodeCountClass::DELTA: {
    uint32_t d = degree_dist(rng) + 1;
    uint32_t base = base_dist(rng);
    node.edge_counts.push_back(base);
    for (uint32_t i = 1; i < d; ++i)
      node.edge_counts.push_back(base + delta_dist(rng));
    break;
  }

  case NodeCountClass::EXPLICIT: {
    uint32_t d = degree_dist(rng) + 1;
    for (uint32_t i = 0; i < d; ++i)
      node.edge_counts.push_back(base_dist(rng));
    break;
  }
  }

  return node;
}

struct DummyOracle {
  std::vector<DummyNode> nodes;

  uint64_t edge_count(uint64_t node, uint32_t edge_rank) const {
    return nodes.at(node).edge_count(edge_rank);
  }

  std::size_t node_count() const { return nodes.size(); }
};

namespace detail {

static inline std::size_t size_in_bytes(const DummyOracle &oracle) {
  std::size_t bytes = 0;
  for (const auto &node : oracle.nodes) {
    bytes += sizeof(uint64_t) * node.outdegree();
  }
  return bytes;
}

} // namespace detail

template <typename URNG>
void fill_oracle_and_builder(DummyOracle &oracle, CountBuffer &builder,
                             const std::vector<NodeCountClass> &types,
                             URNG &rng) {
  for (auto t : types) {
    DummyNode node = make_dummy_node(t, rng);
    oracle.nodes.push_back(node);

    NodeInsertInfo info;
    info.outdegree = node.outdegree();
    info.edge_counts = node.edge_counts.data();

    builder.insert_node(info);
  }
}

static inline DummyOracle
build_oracle_and_builder(CountBuffer &builder,
                         const std::vector<NodeCountClass> &types,
                         uint32_t seed = 123) {
  DummyOracle oracle;
  oracle.nodes.reserve(types.size());

  std::mt19937 rng(seed);
  fill_oracle_and_builder(oracle, builder, types, rng);

  return oracle;
}

static inline void check_counts_equal(const DummyOracle &oracle,
                                      const CompressedCountBuffer &view) {
  ASSERT_EQ(oracle.node_count(), view.node_count());

  for (size_t n = 0; n < oracle.node_count(); ++n) {
    const auto &node = oracle.nodes[n];
    for (size_t e = 0; e < node.outdegree(); ++e) {
      ASSERT_EQ(oracle.edge_count(n, e), view.edge_count(n, e))
          << "Mismatch at node " << n << " edge " << e;
    }
  }
}

static inline void run_counting_test(size_t num_nodes,
                                     NodeTypeDistribution &dist,
                                     uint32_t seed = 12345) {
  auto types = generate_node_types(num_nodes, dist, seed);

  CountBuffer builder;
  builder.reserve(num_nodes);
  EXPECT_TRUE(builder.empty());

  DummyOracle oracle = build_oracle_and_builder(builder, types, seed + 1);
  EXPECT_FALSE(builder.empty());
  EXPECT_EQ(builder.node_count(), num_nodes);

  CompressedCountBuffer compressed(builder);
  builder.clear();
  EXPECT_TRUE(builder.empty());

  check_counts_equal(oracle, compressed);

  LOG_INFO() << "original vector=" << detail::size_in_bytes(oracle)
             << " bytes, compressed vector="
             << detail::size_in_bytes(compressed) << " bytes";
}

static inline void run_append_test(size_t num_nodes_a, size_t num_nodes_b,
                                   NodeTypeDistribution &dist,
                                   uint32_t seed = 12345) {
  auto types_a = generate_node_types(num_nodes_a, dist, seed),
       types_b = generate_node_types(num_nodes_b, dist, seed + 1);

  CountBuffer builder_a;
  builder_a.reserve(num_nodes_a);

  CountBuffer builder_b;
  builder_b.reserve(num_nodes_b);

  DummyOracle oracle;
  oracle.nodes.reserve(num_nodes_a + num_nodes_b);

  std::mt19937 rng(seed + 2);
  fill_oracle_and_builder(oracle, builder_a, types_a, rng);
  fill_oracle_and_builder(oracle, builder_b, types_b, rng);

  builder_a.append(builder_b);

  CompressedCountBuffer compressed(builder_a);
  check_counts_equal(oracle, compressed);

  LOG_INFO() << "original vector=" << detail::size_in_bytes(oracle)
             << " bytes, compressed vector="
             << detail::size_in_bytes(compressed) << " bytes";
}
