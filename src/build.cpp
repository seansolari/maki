
#include <cstdint>
#include <string>

#include "maki/build/io/fasta.hpp"
#include "maki/build/graph/build_colours.hpp"
#include "maki/build/graph/construct_cdbg.hpp"
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
  uint64_t numFeatures = colours.size(); // mark filter colour codes
  auto filters = parse(filterFiles, parseFilterFNA, colours,
                       (std::size_t)params.threads, (std::size_t)params.k);

  // prepare input data
  auto view = combineViews(genomes, filters);
  MetaColours cmap(std::move(colours.ids), numFeatures);

  // suffix-wise buffer construction
  cdbg::BuildOptions ops{
    .kmer_size = params.k,
    .suffix_size = params.s,
    .out = params.outputFolder,
    .pool_size = 2 * params.threads
  };
  auto tmp = cdbg::construct(view, cmap, ops);

  return 0;
}
