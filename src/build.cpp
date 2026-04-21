
#include <cstdint>
#include <string>

#include "maki/build/graph/build_colours.hpp"
#include "maki/build/graph/construct_cdbg.hpp"
#include "maki/build/io/fasta.hpp"
#include "maki/core/seq/io.hpp"

struct BuildParameters {
  std::string queryFile;
  std::string outputFolder;
  uint8_t k = 31;
  uint8_t s = 8;
  uint64_t threads = 16;
  bool isFilter = false;
};

int gff_main(const BuildParameters &params) {
  auto files = readFilePaths(params.queryFile.data(), Gff3FileType);

  // parse input sequences
  Colours colours;
  auto genomes = parse(files, parseGFF, colours, (std::size_t)params.k);

  // suffix-wise buffer construction
  auto g = cdbg::construct(toView(genomes), MetaColours(std::move(colours.ids)),
                           {.kmer_size = params.k,
                            .suffix_size = params.s,
                            .out = params.outputFolder,
                            .pool_size = 2 * params.threads});

  return 0;
}

int filter_main(const BuildParameters &params) {
  auto files = readFilePaths(params.queryFile.data(), FastaFileType);

  // parse input sequences
  Colours colours;
  auto filters = parse(files, parseFilterFNA, colours,
                       (std::size_t)params.threads, (std::size_t)params.k);

  // suffix-wise buffer construction
  cdbg::BuildOptions ops;
  auto g = cdbg::construct(toView(filters), MetaColours(std::move(colours.ids)),
                           {.kmer_size = params.k,
                            .suffix_size = params.s,
                            .out = params.outputFolder,
                            .pool_size = 2 * params.threads});

  return 0;
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char const *argv[]) {
  BuildParameters params{};
  if (params.isFilter) {
    return filter_main(params);
  } else {
    return gff_main(params);
  }
}
