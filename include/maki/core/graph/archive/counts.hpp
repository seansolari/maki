
#pragma once
#include "maki/core/utils/wap_vector.hpp"
#include <cstdint>
#include <sdsl/enc_vector.hpp>
#include <sdsl/rrr_vector.hpp>

/*
 * Node classification following Karasikov et al. - each node belongs to one
 * class.
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
 */
class NodeClassifier {
public:
  NodeCountClass classify(const NodeInsertInfo &n) const;

private:
  static constexpr uint32_t max_delta_degree = 4;
  bool deltas_fit(const NodeInsertInfo &n) const;
};

class CountBuffer {
  friend class CompressedCountBuffer;

public:
  void reserve(std::size_t size_);
  void clear();
  bool empty() const noexcept;
  std::size_t size() const;
  void insert_node(const NodeInsertInfo &node);
  void append(const CountBuffer &rhs);

protected:
  sdsl::bit_vector is_single, is_uniform, is_delta, is_explicit;
  wo_wap_vector base_counts;              // VLE
  std::vector<uint64_t> delta_offsets;    // monotone
  std::vector<int16_t> delta_counts;      // RLE(0) + zig-zag
  std::vector<uint64_t> explicit_offsets; // monotone
  wo_wap_vector explicit_counts;          // VLE
  NodeClassifier classifier_;

  void encode_counts(NodeCountClass cls, const NodeInsertInfo &node);
};

class CompressedCountBuffer {
  CompressedCountBuffer(const CountBuffer &raw_);
  CompressedCountBuffer(CountBuffer &&raw_);
  std::size_t size() const;
  uint64_t edge_count(uint64_t node, uint32_t edge_rank) const;

private:
  sdsl::rrr_vector<> is_single, is_uniform, is_delta, is_explicit;
  sdsl::rrr_vector<>::rank_1_type rs_delta, rs_explicit;
  ro_wap_vector base_counts;           // WAP
  sdsl::enc_vector<> delta_offsets;    // monotone
  std::vector<int16_t> delta_counts;   // RLE(0) + zig-zag
  sdsl::enc_vector<> explicit_offsets; // monotone
  ro_wap_vector explicit_counts;       // WAP

  void init_support();
};
