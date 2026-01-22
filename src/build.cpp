
#include <cstdint>
#include <string>

#include "maki/core/seq/seq_io.hpp"
#include "maki/build/io/fasta.hpp"
#include "maki/colour/encoder.hpp"

struct BuildParameters
{
  std::string queryFile;
  std::string filterFile;
  std::string outputFolder;
  uint8_t k = 31;
  uint8_t s = 8;
  uint64_t threads = 16;
};

int main(int argc, char const *argv[])
{
  BuildParameters params{};
  auto genomeFiles = readFilePaths(params.queryFile.data(), Gff3FileType);
  auto filterFiles = readFilePaths(params.filterFile.data(), FastaFileType)
  
  // parse input sequences
  Colours colours;
  auto genomes = parseWithColour(genomeFiles, colours, parseGFF, params.k);
  auto filters = parseWithColour(filterFiles, colours, parseFilterFNA, params.threads, params.k);

  // suffix-wise construct graph buffers
  auto view = combineViews(genomes, filters);

  return 0;
}
