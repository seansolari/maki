
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/graph/base.hpp"
#include <cereal/archives/binary.hpp>

ColouredGraphFiles::ColouredGraphFiles(fs::path base_)
    : DeBruijnGraphFiles(base_), archive(base / "6.dat") {}

void ColouredGraph::FromDisk(ColouredGraph &g, const std::string &path) {
  ColouredGraphFiles files(path);

  // Load metadata
  {
    std::ifstream is(files.meta, std::ios::binary);
    cereal::BinaryInputArchive iarchive(is);
    iarchive(g);
  }

  LoadBaseBuffers(g, files);

  // Initialise colour archive reader
  g.carch = std::make_unique<ArchiveReader>(files.archive);
}
