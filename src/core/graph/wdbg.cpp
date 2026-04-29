
#include "maki/core/graph/wdbg.hpp"

WeightedGraphFiles WeightedGraph::GraphFiles(fs::path base) {
  return WeightedGraphFiles{.l = base / "1.dat",
                            .lR = base / "2.dat",
                            .lS = base / "3.dat",
                            .W = base / "4.dat",
                            .meta = base / "6.dat"};
}

void WeightedGraph::FromDisk(WeightedGraph &g, const std::string &path) {
  auto files = GraphFiles(path);

  // Load metadata
  {
    std::ifstream is(files.meta, std::ios::binary);
    cereal::BinaryInputArchive iarchive(is);
    iarchive(g);
  }

  // Load `l` arrays
  sdsl::load_from_file(g.l, files.l);

  sdsl::load_from_file(g.lR, files.lR);
  g.lR.set_vector(&g.l);

  sdsl::load_from_file(g.lS, files.lS);
  g.lS.set_vector(&g.l);

  // Load `W` matrix
  sdsl::load_from_file(g.W, files.W);
}

