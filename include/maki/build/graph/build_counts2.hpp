
#pragma once
#include <cstdint>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/enc_vector.hpp>
#include <vector>

struct RawCountBuffer {
  sdsl::bit_vector is_single;
  sdsl::bit_vector is_uniform;
  sdsl::bit_vector is_delta;
  sdsl::bit_vector is_explicit;

  std::vector<uint32_t> base_counts; // VLE

  std::vector<uint64_t> delta_offsets; // monotone
  std::vector<int16_t> delta_counts; // RLE(0) + zig-zag

  std::vector<uint64_t> explicit_offsets; // monotone
  std::vector<uint32_t> explicit_counts; // VLE
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

struct CountBufferConstructor {
  void insert_node(const NodeInsertInfo &node) {
    NodeCountClass cls = classifier_.classify(node);

    raw_.is_single.push_back(0);
    raw_.is_uniform.push_back(0);
    raw_.is_delta.push_back(0);
    raw_.is_explicit.push_back(0);

    switch (cls) {
    case NodeCountClass::SINGLE:
      raw_.is_single.back() = 1;
      break;
    case NodeCountClass::UNIFORM:
      raw_.is_uniform.back() = 1;
      break;
    case NodeCountClass::DELTA:
      raw_.is_delta.back() = 1;
      break;
    case NodeCountClass::EXPLICIT:
      raw_.is_explicit.back() = 1;
      break;
    default:
      break;
    }

    encode_counts(cls, node);
  }

protected:
  RawCountBuffer raw_;
  NodeClassifier classifier_;

  void encode_counts(NodeCountClass cls, const NodeInsertInfo &node) {
    switch (cls) {
    case NodeCountClass::SINGLE:
    case NodeCountClass::UNIFORM:
      raw_.base_counts.push_back(node.edge_counts[0]);
      break;

    case NodeCountClass::DELTA: {
      uint32_t base = node.edge_counts[0];
      raw_.base_counts.push_back(base);
      raw_.delta_offsets.push_back(raw_.delta_counts.size());
      for (uint32_t i = 1; i < node.outdegree; ++i)
        raw_.delta_counts.push_back(int16_t(node.edge_counts[i] - base));
      break;
    }

    case NodeCountClass::EXPLICIT: {
      raw_.explicit_offsets.push_back(raw_.explicit_counts.size());
      for (uint32_t i = 0; i < node.outdegree; ++i)
        raw_.explicit_counts.push_back(node.edge_counts[i]);
      break;
    }
    
    default:
      break;
    }
  }
};

// Width-adaptive packing vector
class ro_wap_vector {
public:
  explicit ro_wap_vector(const std::vector<uint64_t> &v);

  uint64_t operator[](std::size_t i) const {
    if (is_u8[i])
      return u8[rs_u8(i)];
    if (is_u16[i])
      return u16[rs_u16(i)];
    if (is_u32[i])
      return u32[rs_u32(i)];
    if (is_u64[i])
      return u64[rs_u64(i)];
    return 0;
  }

protected:
  sdsl::rrr_vector<> is_u8;
  sdsl::rrr_vector<> is_u16;
  sdsl::rrr_vector<> is_u32;
  sdsl::rrr_vector<> is_u64;

  sdsl::rrr_vector<>::rank_1_type rs_u8;
  sdsl::rrr_vector<>::rank_1_type rs_u16;
  sdsl::rrr_vector<>::rank_1_type rs_u32;
  sdsl::rrr_vector<>::rank_1_type rs_u64;

  std::vector<uint8_t> u8;
  std::vector<uint16_t> u16;
  std::vector<uint32_t> u32;
  std::vector<uint64_t> u64;
};

struct CompressedCountBuffer {
  // Node class masks
  sdsl::rrr_vector<> is_single;
  sdsl::rrr_vector<> is_uniform;
  sdsl::rrr_vector<> is_delta;
  sdsl::rrr_vector<> is_explicit;

  // Rank supports
  sdsl::rrr_vector<>::rank_1_type rs_delta;
  sdsl::rrr_vector<>::rank_1_type rs_explicit;

  ro_wap_vector base_counts; // WAP

  sdsl::enc_vector<> delta_offsets; // monotone
  std::vector<int16_t> delta_counts; // RLE(0) + zig-zag

  sdsl::enc_vector<> explicit_offsets; // monotone
  ro_wap_vector explicit_counts; // WAP
};

CompressedCountBuffer compress_counts(const RawCountBuffer &raw);

class CompressedCountsView {
public:
  explicit CompressedCountsView(const CompressedCountBuffer &d) : data_(d) {}

  uint64_t edge_count(uint64_t node, uint32_t edge_rank) const {
    std::size_t exr = data_.rs_explicit(node), bci = node - exr;

    if (data_.is_single[node]) {
      return data_.base_counts[bci];
    }

    if (data_.is_uniform[node]) {
      return data_.base_counts[bci];
    }

    if (data_.is_delta[node]) {
      uint64_t base = data_.base_counts[bci];
      if (edge_rank == 0)
        return base;
      return uint64_t(int64_t(base) + data_.delta_counts[data_.delta_offsets[data_.rs_delta(node)] + edge_rank - 1]);
    }

    if (data_.is_explicit[node]) {
      return data_.explicit_counts[data_.explicit_offsets[exr] + edge_rank];
    }

    return 0;
  }

private:
  const CompressedCountBuffer &data_;
};
