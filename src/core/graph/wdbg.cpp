
#include "maki/core/graph/wdbg.hpp"


void WeightedGraph::FromDisk(WeightedGraph &g, const std::string &path) {
  WeightedGraphFiles files(path);

  // Load metadata
  {
    std::ifstream is(files.meta, std::ios::binary);
    cereal::BinaryInputArchive iarchive(is);
    iarchive(g);
  }

  LoadBaseBuffers(g, files);
}

