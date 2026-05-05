
#include "maki/core/graph/archive/counts.hpp"
#include "maki/core/utils/algo.hpp"
#include "maki/core/utils/wap_vector.hpp"
#include <cstdint>
#include <sdsl/io.hpp>

NodeCountClass NodeClassifier::classify(const NodeInsertInfo &n) const {
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

bool NodeClassifier::deltas_fit(const NodeInsertInfo &n) const {
  int64_t base = static_cast<int64_t>(n.edge_counts[0]);
  for (uint8_t i = 1; i < n.outdegree; ++i) {
    int64_t d = static_cast<int64_t>(n.edge_counts[i]) - base;
    if (d < -32768 || d > 32767)
      return false;
  }
  return true;
}

void CountBuffer::reserve(std::size_t size_) {
  is_single.reserve(size_);
  is_uniform.reserve(size_);
  is_delta.reserve(size_);
  is_explicit.reserve(size_);
}

void CountBuffer::clear() {
  is_single.clear();
  is_uniform.clear();
  is_delta.clear();
  is_explicit.clear();
  base_counts.clear();
  delta_offsets.clear();
  delta_counts.clear();
  explicit_offsets.clear();
  explicit_counts.clear();
}

bool CountBuffer::empty() const noexcept {
  assert(is_uniform.empty() == is_single.empty());
  assert(is_delta.empty() == is_single.empty());
  assert(is_explicit.empty() == is_single.empty());
  return is_single.empty();
}

std::size_t CountBuffer::node_count() const {
  assert(is_uniform.size() == is_single.size());
  assert(is_delta.size() == is_single.size());
  assert(is_explicit.size() == is_single.size());
  return is_single.size();
}

void CountBuffer::insert_node(const NodeInsertInfo &node) {
  NodeCountClass cls = classifier_.classify(node);

  is_single.push_back(0);
  is_uniform.push_back(0);
  is_delta.push_back(0);
  is_explicit.push_back(0);

  switch (cls) {
  case NodeCountClass::SINGLE:
    is_single.back() = 1;
    break;
  case NodeCountClass::UNIFORM:
    is_uniform.back() = 1;
    break;
  case NodeCountClass::DELTA:
    is_delta.back() = 1;
    break;
  case NodeCountClass::EXPLICIT:
    is_explicit.back() = 1;
    break;
  default:
    break;
  }

  encode_counts(cls, node);
}

void CountBuffer::encode_counts(NodeCountClass cls,
                                const NodeInsertInfo &node) {
  switch (cls) {
  case NodeCountClass::SINGLE:
  case NodeCountClass::UNIFORM:
    base_counts.push_back(node.edge_counts[0]);
    break;

  case NodeCountClass::DELTA: {
    uint64_t base = node.edge_counts[0];
    base_counts.push_back(base);
    delta_offsets.push_back(delta_counts.size());
    for (uint8_t i = 1; i < node.outdegree; ++i)
      delta_counts.push_back(int16_t(node.edge_counts[i] - base));
    break;
  }

  case NodeCountClass::EXPLICIT: {
    explicit_offsets.push_back(explicit_counts.size());
    for (uint8_t i = 0; i < node.outdegree; ++i)
      explicit_counts.push_back(node.edge_counts[i]);
    break;
  }

  default:
    break;
  }
}

void CountBuffer::append(const CountBuffer &rhs) {
  vector_append(is_single, rhs.is_single);
  vector_append(is_uniform, rhs.is_uniform);
  vector_append(is_delta, rhs.is_delta);
  vector_append(is_explicit, rhs.is_explicit);
  base_counts.append(rhs.base_counts);
  vector_append_add(delta_offsets, rhs.delta_offsets, delta_counts.size());
  vector_append(delta_counts, rhs.delta_counts);
  vector_append_add(explicit_offsets, rhs.explicit_offsets,
                    explicit_counts.size());
  explicit_counts.append(rhs.explicit_counts);
}

CompressedCountBuffer::CompressedCountBuffer(const CountBuffer &raw_)
    : is_single(raw_.is_single), is_uniform(raw_.is_uniform),
      is_delta(raw_.is_delta), is_explicit(raw_.is_explicit), rs_delta(),
      rs_explicit(), base_counts(raw_.base_counts),
      delta_offsets(raw_.delta_offsets), delta_counts(raw_.delta_counts),
      explicit_offsets(raw_.explicit_offsets),
      explicit_counts(raw_.explicit_counts) {
  init_support();
}

CompressedCountBuffer::CompressedCountBuffer(CountBuffer &&raw_)
    : is_single(raw_.is_single), is_uniform(raw_.is_uniform),
      is_delta(raw_.is_delta), is_explicit(raw_.is_explicit), rs_delta(),
      rs_explicit(), base_counts(std::move(raw_.base_counts)),
      delta_offsets(raw_.delta_offsets),
      delta_counts(std::move(raw_.delta_counts)),
      explicit_offsets(raw_.explicit_offsets),
      explicit_counts(std::move(raw_.explicit_counts)) {
  init_support();
}

std::size_t CompressedCountBuffer::node_count() const {
  assert(is_uniform.size() == is_single.size());
  assert(is_delta.size() == is_single.size());
  assert(is_explicit.size() == is_single.size());
  return is_single.size();
}

uint64_t CompressedCountBuffer::edge_count(uint64_t node,
                                           uint32_t edge_rank) const {
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

void CompressedCountBuffer::init_support() {
  sdsl::util::init_support(rs_delta, &is_delta);
  sdsl::util::init_support(rs_explicit, &is_explicit);
}

std::size_t detail::size_in_bytes(const CountBuffer &vec) {
  return sdsl::size_in_bytes(vec.is_single) +
         sdsl::size_in_bytes(vec.is_uniform) +
         sdsl::size_in_bytes(vec.is_delta) +
         sdsl::size_in_bytes(vec.is_explicit) +
         detail::size_in_bytes(vec.base_counts) +
         (sizeof(uint64_t) * vec.delta_offsets.capacity()) +
         (sizeof(int16_t) * vec.delta_counts.capacity()) +
         (sizeof(uint64_t) * vec.explicit_offsets.capacity()) +
         detail::size_in_bytes(vec.explicit_counts);
}

std::size_t detail::size_in_bytes(const CompressedCountBuffer &vec) {
  return sdsl::size_in_bytes(vec.is_single) +
         sdsl::size_in_bytes(vec.is_uniform) +
         sdsl::size_in_bytes(vec.is_delta) +
         sdsl::size_in_bytes(vec.is_explicit) +
         sdsl::size_in_bytes(vec.rs_delta) +
         sdsl::size_in_bytes(vec.rs_explicit) +
         detail::size_in_bytes(vec.base_counts) +
         sdsl::size_in_bytes(vec.delta_offsets) +
         (sizeof(int16_t) * vec.delta_counts.capacity()) +
         sdsl::size_in_bytes(vec.explicit_offsets) +
         detail::size_in_bytes(vec.explicit_counts);
}
