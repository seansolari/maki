
#pragma once
#include "base.hpp"
#include "colours.hpp"
#include "maki/core/graph/archive/archive_reader.hpp"
#include <cereal/types/array.hpp>
#include <cereal/types/base_class.hpp>
#include <filesystem>

namespace fs = std::filesystem;

struct ColouredGraphFiles : public DeBruijnGraphFiles {
  ColouredGraphFiles(fs::path base_);
  fs::path archive;
};

struct ColouredGraph : public DeBruijnGraph {
  // Colour data
  std::unique_ptr<ArchiveReader> carch;
  ColourRegistry cmap;

  // load from disk
  static void FromDisk(ColouredGraph &g, const std::string &path);

private:
  // serialise to disk
  friend class cereal::access;
  template <class Archive> void serialize(Archive &ar) {
    ar(cereal::base_class<DeBruijnGraph>(this), cmap);
  }
};
