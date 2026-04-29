
#pragma once
#include "base.hpp"
#include <cereal/types/array.hpp>
#include <filesystem>

namespace fs = std::filesystem;

struct WeightedGraphFiles {
  fs::path l, lR, lS, W, meta;
};

struct WeightedGraph : public DeBruijnGraph {
  static WeightedGraphFiles GraphFiles(fs::path base);
  static void FromDisk(WeightedGraph &g, const std::string &path);

private:
  // serialise to disk
  friend class cereal::access;
  template <class Archive> void serialize(Archive &ar) { ar(k, F, C); }
};
