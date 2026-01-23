
#include <cstdint>
#include <string>

#include "maki/build/graph/build_colours.hpp"
#include "maki/build/io/fasta.hpp"
#include "maki/core/seq/seq_io.hpp"

struct BuildParameters {
  std::string queryFile;
  std::string filterFile;
  std::string outputFolder;
  uint8_t k = 31;
  uint8_t s = 8;
  uint64_t threads = 16;
};

int main([[maybe_unused]] int argc, [[maybe_unused]] char const *argv[]) {
  BuildParameters params{};
  auto genomeFiles = readFilePaths(params.queryFile.data(), Gff3FileType);
  auto filterFiles = readFilePaths(params.filterFile.data(), FastaFileType);

  // parse input sequences
  Colours colours;
  auto genomes = parse(genomeFiles, parseGFF, colours, (std::size_t)params.k);
  auto filters = parse(filterFiles, parseFilterFNA, colours,
                       (std::size_t)params.threads, (std::size_t)params.k);

  // suffix-wise construct graph buffers
  auto view = combineViews(genomes, filters);

  return 0;
}
