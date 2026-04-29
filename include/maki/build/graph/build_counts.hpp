
#pragma once
#include <cstddef>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/enc_vector.hpp>

struct RawCountConstruction {
  // Node-wise classification (raw)
  sdsl::bit_vector is_single;
  sdsl::bit_vector is_uniform;
  sdsl::bit_vector is_delta;
  sdsl::bit_vector is_explicit;

  // Offsets into value streams (monotone)
  std::vector<uint64_t> offsets;

  // Value streams (still uncompressed here)
  std::vector<uint32_t> base_counts;
  std::vector<int16_t> delta_counts;
  std::vector<uint32_t> explicit_counts;

  size_t node_count = 0;

  void begin_nodes(size_t expected_nodes) {
    is_single.resize(expected_nodes);
    is_uniform.resize(expected_nodes);
    is_delta.resize(expected_nodes);
    is_explicit.resize(expected_nodes);
  }
};

struct CompressedCountData {
  // Node class masks
  sdsl::rrr_vector<> is_single;
  sdsl::rrr_vector<> is_uniform;
  sdsl::rrr_vector<> is_delta;
  sdsl::rrr_vector<> is_explicit;

  // Rank supports
  sdsl::rrr_vector<>::rank_1_type rs_single;
  sdsl::rrr_vector<>::rank_1_type rs_uniform;
  sdsl::rrr_vector<>::rank_1_type rs_delta;
  sdsl::rrr_vector<>::rank_1_type rs_explicit;

  // Offsets compressed with Elias–Fano
  sdsl::enc_vector<> offsets;

  // Value streams (optionally packed later)
  std::vector<uint32_t> base_counts;     // width-adaptive encoding
  std::vector<int16_t> delta_counts;     // RLE(0) + zig-zag
  std::vector<uint32_t> explicit_counts; // wavelet tree
};

/*
 * Node classification following Karasikov et al.
 * Each node chooses exactly one class.
 */
enum class NodeCountClass : uint8_t {
  NONE = 0,    // outdegree == 0
  SINGLE = 1,  // outdegree == 1
  UNIFORM = 2, // outdegree >1, all counts equal
  DELTA = 3,   // small fan-out, store base + small deltas
  EXPLICIT = 4 // fallback
};

/*
 * Input for inserting a node.
 * The builder assumes:
 *  - nodes arrive in BOSS order
 *  - outgoing edges are already in BOSS edge order
 */
struct NodeInsertInfo {
  uint32_t outdegree;
  const uint32_t *edge_counts; // length = outdegree
};

/*
 * Policy object deciding how a node's counts are encoded.
 * Kept separate so heuristics can evolve.
 */
class NodeClassifier {
public:
  NodeCountClass classify(const NodeInsertInfo &n) const {
    if (n.outdegree == 0)
      return NodeCountClass::NONE;

    if (n.outdegree == 1)
      return NodeCountClass::SINGLE;

    bool uniform = true;
    for (uint32_t i = 1; i < n.outdegree; ++i) {
      if (n.edge_counts[i] != n.edge_counts[0]) {
        uniform = false;
        break;
      }
    }
    if (uniform)
      return NodeCountClass::UNIFORM;

    if (n.outdegree <= max_delta_degree && deltas_fit(n))
      return NodeCountClass::DELTA;

    return NodeCountClass::EXPLICIT;
  }

private:
  static constexpr uint32_t max_delta_degree = 4;

  bool deltas_fit(const NodeInsertInfo &n) const {
    int32_t base = static_cast<int32_t>(n.edge_counts[0]);
    for (uint32_t i = 1; i < n.outdegree; ++i) {
      int32_t d = static_cast<int32_t>(n.edge_counts[i]) - base;
      if (d < -32768 || d > 32767)
        return false;
    }
    return true;
  }
};

class CountingBOSSCountsBuilder {
public:
  explicit CountingBOSSCountsBuilder(size_t expected_nodes) {
    raw_.begin_nodes(expected_nodes);
  }

  void insert_node(const NodeInsertInfo &node) {
    const size_t i = raw_.node_count++;

    NodeCountClass cls = classifier_.classify(node);

    // clear bits
    raw_.is_single[i] = 0;
    raw_.is_uniform[i] = 0;
    raw_.is_delta[i] = 0;
    raw_.is_explicit[i] = 0;

    // set one class bit
    switch (cls) {
    case NodeCountClass::SINGLE:
      raw_.is_single[i] = 1;
      break;
    case NodeCountClass::UNIFORM:
      raw_.is_uniform[i] = 1;
      break;
    case NodeCountClass::DELTA:
      raw_.is_delta[i] = 1;
      break;
    case NodeCountClass::EXPLICIT:
      raw_.is_explicit[i] = 1;
      break;
    default:
      break;
    }

    raw_.offsets.push_back(current_offset(cls));
    encode_counts(cls, node);
  }

  const RawCountConstruction &raw_data() const { return raw_; }

private:
  RawCountConstruction raw_;
  NodeClassifier classifier_;

  uint64_t current_offset(NodeCountClass cls) const {
    if (cls == NodeCountClass::EXPLICIT)
      return raw_.explicit_counts.size();
    if (cls == NodeCountClass::DELTA || cls == NodeCountClass::SINGLE ||
        cls == NodeCountClass::UNIFORM)
      return raw_.base_counts.size();
    return 0;
  }

  void encode_counts(NodeCountClass cls, const NodeInsertInfo &n) {
    switch (cls) {
    case NodeCountClass::SINGLE:
    case NodeCountClass::UNIFORM:
      raw_.base_counts.push_back(n.edge_counts[0]);
      break;

    case NodeCountClass::DELTA: {
      uint32_t base = n.edge_counts[0];
      raw_.base_counts.push_back(base);
      for (uint32_t i = 1; i < n.outdegree; ++i)
        raw_.delta_counts.push_back(int16_t(n.edge_counts[i] - base));
      break;
    }

    case NodeCountClass::EXPLICIT:
      for (uint32_t i = 0; i < n.outdegree; ++i)
        raw_.explicit_counts.push_back(n.edge_counts[i]);
      break;

    default:
      break;
    }
  }
};

CompressedCountData compress_counts(const RawCountConstruction &raw) {
  CompressedCountData c;

  c.is_single = sdsl::rrr_vector<>(raw.is_single);
  c.is_uniform = sdsl::rrr_vector<>(raw.is_uniform);
  c.is_delta = sdsl::rrr_vector<>(raw.is_delta);
  c.is_explicit = sdsl::rrr_vector<>(raw.is_explicit);

  sdsl::util::init_support(c.rs_single, &c.is_single);
  sdsl::util::init_support(c.rs_uniform, &c.is_uniform);
  sdsl::util::init_support(c.rs_delta, &c.is_delta);
  sdsl::util::init_support(c.rs_explicit, &c.is_explicit);

  c.offsets = sdsl::enc_vector<>(raw.offsets.begin(), raw.offsets.end());

  // (Optional) further compression of value streams happens here
  c.base_counts = raw.base_counts;
  c.delta_counts = raw.delta_counts;
  c.explicit_counts = raw.explicit_counts;

  return c;
}

class CountingBOSSCountsView {
public:
  explicit CountingBOSSCountsView(const CompressedCountData &d) : data_(d) {}

  uint32_t edge_count(uint64_t node, uint32_t edge_rank) const {
    if (data_.is_single[node]) {
      return data_.base_counts[data_.offsets[node]];
    }

    if (data_.is_uniform[node]) {
      return data_.base_counts[data_.offsets[node]];
    }

    if (data_.is_delta[node]) {
      uint64_t off = data_.offsets[node];
      if (edge_rank == 0)
        return data_.base_counts[off];
      return uint32_t(int32_t(data_.base_counts[off]) +
                      data_.delta_counts[off + edge_rank - 1]);
    }

    if (data_.is_explicit[node]) {
      return data_.explicit_counts[data_.offsets[node] + edge_rank];
    }

    return 0;
  }

private:
  const CompressedCountData &data_;
};
