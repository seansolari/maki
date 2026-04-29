
#pragma once
#include "base.hpp"
#include "colours.hpp"
#include "maki/core/graph/archive/archive_reader.hpp"
#include <cereal/types/array.hpp>
#include <filesystem>

namespace fs = std::filesystem;

struct ColouredGraphFiles {
  fs::path l, lR, lS, W, archive, meta;
};

struct ColouredGraph : public DeBruijnGraph {
  // Colour data
  std::unique_ptr<ArchiveReader> carch;
  ColourRegistry cmap;

  // load from disk
  static ColouredGraphFiles GraphFiles(fs::path base);
  static void FromDisk(ColouredGraph &g, const std::string &path);

private:
  // serialise to disk
  friend class cereal::access;
  template <class Archive> void serialize(Archive &ar) { ar(k, cmap, F, C); }
};
