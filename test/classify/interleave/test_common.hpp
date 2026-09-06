#pragma once
#include "maki/build/graph/construct_cdbg.hpp"
#include "maki/build/io/fasta.hpp"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <string>

using ::testing::ContainerEq;
using ::testing::ElementsAreArray;

inline std::vector<Dna4Genome> ReadGenome(const std::string &data, Colours &c,
                                          size_t k) {
  std::vector<Dna4Genome> result(1);
  std::istringstream datastream(data);
  parseFastaStream(result[0], datastream, c, k);
  return result;
}

inline void MakeGraph(const std::vector<Dna4Genome> &fna, Colours &c, size_t k,
                      size_t s, fs::path &bufferPath) {
  cdbg::construct(toView(fna), std::move(c.ids),
                  {k, s, bufferPath, 1});
}

inline void MakeGraph(std::string data, size_t k, size_t s,
                      fs::path &bufferPath) {
  Colours c;
  std::vector<Dna4Genome> fna = ReadGenome(data, c, k);
  MakeGraph(fna, c, k, s, bufferPath);
}

inline void MakeGraph(const std::vector<std::string> &kmers, size_t k, size_t s,
                      fs::path &bufferPath) {
  for (const auto &kmer : kmers) {
    ASSERT_EQ(kmer.size(), k + 1);
  }
  MakeGraph(">header\n" +
                std::accumulate(kmers.begin(), kmers.end(), std::string("N")) +
                "\n",
            k, s, bufferPath);
}

inline static std::string random_dna_sequence(size_t length, uint32_t seed) {
  static const char alphabet[] = {'A', 'C', 'G', 'T'};
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> dist(0, 3);

  std::string s;
  s.reserve(length);
  for (size_t i = 0; i < length; ++i) {
    s.push_back(alphabet[dist(rng)]);
  }
  return s;
}

inline static std::string apply_snp(const std::string &seq, std::size_t position) {
  std::string mutated = seq;
  assert(position < mutated.size());
  for (char c : std::array<char,4>{'A', 'C', 'G', 'T'}) {
    if (mutated[position] != c) {
      mutated[position] = c;
      break;
    }
  }
  return mutated;
}
