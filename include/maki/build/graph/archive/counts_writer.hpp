
#pragma once
#include <cstdint>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/util.hpp>
#include <vector>

// Width-adaptive packing vector
class wo_wap_vector {
  friend class ro_wap_vector;

public:
  wo_wap_vector() = default;

  explicit wo_wap_vector(std::size_t reserve_) {
    is_u8.reserve(reserve_);
    is_u16.reserve(reserve_);
    is_u32.reserve(reserve_);
    is_u64.reserve(reserve_);
  }

  std::size_t size() const {
    assert(is_u16.size() == is_u8.size());
    assert(is_u32.size() == is_u8.size());
    assert(is_u64.size() == is_u8.size());
    return is_u8.size();
  }

  void push_back(uint64_t value_) {
    is_u8.push_back(0);
    is_u16.push_back(0);
    is_u32.push_back(0);
    is_u64.push_back(0);

    if (value_ <= std::numeric_limits<uint8_t>::max()) {
      is_u8.back() = 1;
      u8.push_back(value_);
    } else if (value_ <= std::numeric_limits<uint16_t>::max()) {
      is_u16.back() = 1;
      u16.push_back(value_);
    } else if (value_ <= std::numeric_limits<uint32_t>::max()) {
      is_u32.back() = 1;
      u32.push_back(value_);
    } else {
      is_u64.back() = 1;
      u64.push_back(value_);
    }
  }

private:
  sdsl::bit_vector is_u8, is_u16, is_u32, is_u64;
  std::vector<uint8_t> u8;
  std::vector<uint16_t> u16;
  std::vector<uint32_t> u32;
  std::vector<uint64_t> u64;
};

struct RawCountBuffer {
  sdsl::bit_vector is_single, is_uniform, is_delta, is_explicit;

  wo_wap_vector base_counts; // VLE

  std::vector<uint64_t> delta_offsets; // monotone
  std::vector<int16_t> delta_counts;   // RLE(0) + zig-zag

  std::vector<uint64_t> explicit_offsets; // monotone
  wo_wap_vector explicit_counts;          // VLE
};

/*
 * Node classification following Karasikov et al.
 * Each node belongs to one class.
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
  uint8_t outdegree;
  const uint64_t *edge_counts; // length = outdegree
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
    int64_t base = static_cast<int64_t>(n.edge_counts[0]);
    for (uint8_t i = 1; i < n.outdegree; ++i) {
      int64_t d = static_cast<int64_t>(n.edge_counts[i]) - base;
      if (d < -32768 || d > 32767)
        return false;
    }
    return true;
  }
};

struct CountsPayload {
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
      uint64_t base = node.edge_counts[0];
      raw_.base_counts.push_back(base);
      raw_.delta_offsets.push_back(raw_.delta_counts.size());
      for (uint8_t i = 1; i < node.outdegree; ++i)
        raw_.delta_counts.push_back(int16_t(node.edge_counts[i] - base));
      break;
    }

    case NodeCountClass::EXPLICIT: {
      raw_.explicit_offsets.push_back(raw_.explicit_counts.size());
      for (uint8_t i = 0; i < node.outdegree; ++i)
        raw_.explicit_counts.push_back(node.edge_counts[i]);
      break;
    }

    default:
      break;
    }
  }
};

class ro_wap_vector {
public:
  ro_wap_vector(const wo_wap_vector &v)
      : is_u8(v.is_u8), is_u16(v.is_u16), is_u32(v.is_u32), is_u64(v.is_u64),
        rs_u8(), rs_u16(), rs_u32(), rs_u64(), u8(v.u8), u16(v.u16), u32(v.u32),
        u64(v.u64) {
    init_support();
  }

  ro_wap_vector(wo_wap_vector &&v)
      : is_u8(v.is_u8), is_u16(v.is_u16), is_u32(v.is_u32), is_u64(v.is_u64),
        rs_u8(), rs_u16(), rs_u32(), rs_u64(), u8(std::move(v.u8)),
        u16(std::move(v.u16)), u32(std::move(v.u32)), u64(std::move(v.u64)) {
    init_support();
  }

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

  std::size_t size() const {
    assert(is_u16.size() == is_u8.size());
    assert(is_u32.size() == is_u8.size());
    assert(is_u64.size() == is_u8.size());
    return is_u8.size();
  }

protected:
  sdsl::rrr_vector<> is_u8, is_u16, is_u32, is_u64;
  sdsl::rrr_vector<>::rank_1_type rs_u8, rs_u16, rs_u32, rs_u64;
  std::vector<uint8_t> u8;
  std::vector<uint16_t> u16;
  std::vector<uint32_t> u32;
  std::vector<uint64_t> u64;

  void init_support() {
    sdsl::util::init_support(rs_u8, &is_u8);
    sdsl::util::init_support(rs_u16, &is_u16);
    sdsl::util::init_support(rs_u32, &is_u32);
    sdsl::util::init_support(rs_u64, &is_u64);
  }
};

class CompressedCountBuffer {
  CompressedCountBuffer(const RawCountBuffer &raw_)
      : is_single(raw_.is_single), is_uniform(raw_.is_uniform),
        is_delta(raw_.is_delta), is_explicit(raw_.is_explicit), rs_delta(),
        rs_explicit(), base_counts(raw_.base_counts),
        delta_offsets(raw_.delta_offsets), delta_counts(raw_.delta_counts),
        explicit_offsets(raw_.explicit_offsets),
        explicit_counts(raw_.explicit_counts) {
    init_support();
  }

  CompressedCountBuffer(RawCountBuffer &&raw_)
      : is_single(raw_.is_single), is_uniform(raw_.is_uniform),
        is_delta(raw_.is_delta), is_explicit(raw_.is_explicit), rs_delta(),
        rs_explicit(), base_counts(std::move(raw_.base_counts)),
        delta_offsets(raw_.delta_offsets),
        delta_counts(std::move(raw_.delta_counts)),
        explicit_offsets(raw_.explicit_offsets),
        explicit_counts(std::move(raw_.explicit_counts)) {
    init_support();
  }

  std::size_t size() const {
    assert(is_uniform.size() == is_single.size());
    assert(is_delta.size() == is_single.size());
    assert(is_explicit.size() == is_single.size());
    return is_single.size();
  }

  uint64_t edge_count(uint64_t node, uint32_t edge_rank) const {
    std::size_t exr = rs_explicit(node), bci = node - exr;

    if (is_single[node]) {
      return base_counts[bci];
    }

    if (is_uniform[node]) {
      return base_counts[bci];
    }

    if (is_delta[node]) {
      uint64_t base = base_counts[bci];
      if (edge_rank == 0)
        return base;
      return uint64_t(
          int64_t(base) +
          delta_counts[delta_offsets[rs_delta(node)] + edge_rank - 1]);
    }

    if (is_explicit[node]) {
      return explicit_counts[explicit_offsets[exr] + edge_rank];
    }

    return 0;
  }

private:
  // Node class masks
  sdsl::rrr_vector<> is_single, is_uniform, is_delta, is_explicit;

  // Rank supports
  sdsl::rrr_vector<>::rank_1_type rs_delta, rs_explicit;

  ro_wap_vector base_counts; // WAP

  sdsl::enc_vector<> delta_offsets;  // monotone
  std::vector<int16_t> delta_counts; // RLE(0) + zig-zag

  sdsl::enc_vector<> explicit_offsets; // monotone
  ro_wap_vector explicit_counts;       // WAP

  void init_support() {
    sdsl::util::init_support(rs_delta, &is_delta);
    sdsl::util::init_support(rs_explicit, &is_explicit);
  }
};
