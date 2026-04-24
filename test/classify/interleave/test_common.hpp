#pragma once
#include "maki/build/io/fasta.hpp"
#include "maki/build/graph/construct_cdbg.hpp"
#include <random>
#include <string>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

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
                  {.kmer_size = k, .suffix_size = s, .out = bufferPath});
}

inline void MakeGraph(std::string data, size_t k, size_t s,
                      fs::path &bufferPath) {
  Colours c;
  std::vector<Dna4Genome> fna = ReadGenome(data, c, k);
  MakeGraph(fna, c, k, s, bufferPath);
}

inline std::vector<std::string>
GenerateRandomKmers(size_t k, size_t count, uint32_t seed,
                    const std::string &alphabet = "ACGT") {
  std::mt19937 rng(seed);
  std::uniform_int_distribution<> dist(0, alphabet.size() - 1);

  std::unordered_set<std::string> unique;
  while (unique.size() < count) {
    std::string km(k, 'A');
    for (size_t i = 0; i < k; ++i)
      km[i] = alphabet[dist(rng)];
    unique.insert(km);
  }

  return {unique.begin(), unique.end()};
}

// Guarantees controlled overlap between two sets
inline void GenerateOverlappingKmers(size_t k, size_t total,
                                     double overlap_ratio, uint32_t seed,
                                     std::vector<std::string> &out1,
                                     std::vector<std::string> &out2) {
  auto base = GenerateRandomKmers(k, total, seed);

  size_t shared = static_cast<size_t>(total * overlap_ratio);
  out1.assign(base.begin(), base.end());
  out2.assign(base.begin(), base.begin() + shared);

  auto extra = GenerateRandomKmers(k, total - shared, seed + 1);

  out2.insert(out2.end(), extra.begin(), extra.end());
}