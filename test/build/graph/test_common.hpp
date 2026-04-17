#include "maki/build/io/fasta.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using ::testing::ContainerEq;
using ::testing::ElementsAreArray;

inline void ParseFastaToGenome(Dna4Genome &genome, const std::string &file, Colours &c, std::size_t k) {
  zstr::ifstream zis(file);
  std::string fastaData(std::istreambuf_iterator<char>(zis), {});

  // process stream data
  std::istringstream fs(fastaData);
  parseFastaStream(genome, fs, c, k);
}