
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

uint64_t WeightedGraph::edge_count(std::size_t i) const {
  std::size_t n = lR(i), r = i - (lS(n) + 1);
  return occ.edge_count(n, r);
}
