
#include <cstdint>
#include <string>
#include "maki/io/fasta.hpp"
#include "maki/colour/encoder.hpp"

struct BuildParameters
{
  std::string queryFile;
  std::string filterFile;
  std::string outputFolder;
  uint8_t k = 31;
  uint8_t s = 8;
};

int main(int argc, char const *argv[])
{
  BuildParameters params{};
  auto genomeFiles = readFilePaths(params.queryFile.data());
  
  // parse input sequences
  Colours colours;
  auto genomes = loadGenomes(genomeFiles, colours, params.k);
  auto filters = loadFilter(params.filterFile, colours, params.k);

  // suffix-wise construct graph buffers

  return 0;
}
