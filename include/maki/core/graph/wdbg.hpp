
#pragma once
#include "base.hpp"
#include "maki/core/graph/archive/counts.hpp"
#include <cereal/types/array.hpp>

struct WeightedGraphFiles : public DeBruijnGraphFiles {
  WeightedGraphFiles(fs::path base_) : DeBruijnGraphFiles(base_) {}
};

struct WeightedGraph : public DeBruijnGraph {
  // k-mer counts array
  CompressedCountBuffer occ;

  // load from disk
  static void FromDisk(WeightedGraph &g, const std::string &path);

  // count of an edge
  uint64_t edge_count(std::size_t i) const;

private:
  // serialise to disk
  friend class cereal::access;
  template <class Archive> void serialize(Archive &ar) { ar(k, F, C, occ); }
};
