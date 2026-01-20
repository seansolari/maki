
#include <cstdint>
#include <string>
#include "maki/io/fasta.hpp"

struct BuildParameters
{
  std::string queryFile;
  std::string filterFile;
  std::string outputFolder;
  uint8_t k = 0;
  uint8_t s = 0;
};

int main(int argc, char const *argv[])
{
  BuildParameters args{};
  auto genomeFiles = readFilePaths(params.queryFile.data());

  // parse reference genomes


  return 0;
}
