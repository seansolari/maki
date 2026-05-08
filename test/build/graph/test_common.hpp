#include "maki/build/graph/construct_cdbg.hpp"
#include "maki/build/graph/construct_wdbg.hpp"
#include "maki/build/io/fasta.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/graph/wdbg.hpp"

#include "gmock/gmock.h"
#include <filesystem>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <random>
#include <sdsl/construct.hpp>
#include <type_traits>
#include <unordered_map>

using ::testing::ContainerEq;
using ::testing::ElementsAreArray;
using ::testing::Pointwise;
using ::testing::Eq;

static char complement(char c) {
  switch (c) {
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  default:
    return 'N';
  }
}

static inline std::string reverse_complement(const std::string &s) {
  std::string rc(s.rbegin(), s.rend());
  for (auto &c : rc)
    c = complement(c);
  return rc;
}

inline std::string random_dna(size_t length, double repeat_prob,
                              std::mt19937 &rng) {
  static const char bases[] = {'A', 'C', 'G', 'T'};
  std::uniform_int_distribution<int> base_dist(0, 3);
  std::uniform_real_distribution<double> prob(0.0, 1.0);

  std::string seq;
  seq.reserve(length);

  for (size_t i = 0; i < length; ++i) {
    // Introduce repetition: copy previous base
    if (i > 0 && prob(rng) < repeat_prob) {
      seq.push_back(seq[i - 1]);
    } else {
      seq.push_back(bases[base_dist(rng)]);
    }
  }
  return seq;
}

inline std::string motif_repeat(const std::string &motif, size_t length) {
  std::string seq;
  while (seq.size() < length) {
    seq += motif;
  }
  seq.resize(length);
  return seq;
}

inline std::vector<std::string> generate_sequences(size_t n_seqs, size_t length,
                                                   double repeatness,
                                                   std::mt19937 &rng) {
  std::vector<std::string> seqs;
  for (size_t i = 0; i < n_seqs; ++i) {
    seqs.push_back(random_dna(length, repeatness, rng));
  }
  return seqs;
}

inline std::unordered_map<std::string, int>
OracleKmerCount(const std::vector<std::string> &seqs, size_t k) {
  std::unordered_map<std::string, int> counts;

  for (const auto &seq : seqs) {
    if (seq.size() < k)
      continue;

    for (size_t i = 0; i <= seq.size() - k; ++i) {
      std::string kmer = seq.substr(i, k);
      std::string rc = reverse_complement(kmer);

      counts[kmer]++;
      counts[rc]++;
    }
  }

  return counts;
}

inline void BuildGraph(const std::vector<std::string> &seqs, size_t k, size_t s,
                       const fs::path &bufferPath, bool coloured = false) {
  std::string fasta = ">seqA\n";
  for (const auto &seq : seqs)
    fasta += seq + 'N';

  std::vector<Dna4Genome> genomes(1);
  std::istringstream datastream(fasta);
  Colours c;
  parseFastaStream(genomes[0], datastream, c, k);

  if (coloured) {
    cdbg::construct(toView(genomes), MetaColours(std::move(c.ids)),
                    {.kmer_size = k, .suffix_size = s, .out = bufferPath});
  } else {
    wdbg::construct(toView(genomes),
                    {.kmer_size = k, .suffix_size = s, .out = bufferPath});
  }
}

inline void ParseFastaToGenome(Dna4Genome &genome, const std::string &file,
                               Colours &c, std::size_t k) {
  zstr::ifstream zis(file);
  std::string fastaData(std::istreambuf_iterator<char>(zis), {});

  // process stream data
  std::istringstream fs(fastaData);
  parseFastaStream(genome, fs, c, k);
}

inline wavelet_matrix EdgesToWaveletMatrix(const sdsl::int_vector<4> &edges) {
  wavelet_matrix wm;
  sdsl::construct_im(wm, edges);
  return wm;
}

inline std::pair<std::vector<uint64_t>, std::vector<uint64_t>>
UnpackColours(const ColouredGraph &g) {
  std::pair<std::vector<uint64_t>, std::vector<uint64_t>> ca{};

  std::size_t C = g.carch->chunk_count();
  for (std::size_t c = 0; c < C; ++c) {
    auto view = g.carch->view(c);
    std::size_t L = view.size();
    for (std::size_t i = 0; i < L; ++i) {
      auto cols = g.cmap.colours(view.get(i));
      std::visit(
          [&](auto &&arg) {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, colour_t>) {
              // Push a single colour
              ca.first.push_back(arg);
              ca.second.push_back(1);
            } else if constexpr (std::is_same_v<T,
                                                ColourRegistry::vector_ref>) {
              // Push multiple colours
              const auto &colvec = arg.get();
              ca.first.insert(ca.first.end(), colvec.cbegin(), colvec.cend());
              ca.second.push_back(1);
              ca.second.insert(ca.second.end(), colvec.size() - 1, 0u);
            }
          },
          cols);
    }
  }

  return ca;
}

inline bool CheckKmers(const ColouredGraph &g) {
  for (std::size_t e = 0; e < g.edges(); ++e) {
    g.kmer(e);
  }
  return true;
}

inline std::unordered_map<std::string, int>
ExtractKmers(const WeightedGraph &g) {
  std::unordered_map<std::string, int> counts;
  for (std::size_t e = 0; e < g.edges(); ++e) {
    auto kmer = g.kmer(e);
    if (kmer.size() == static_cast<std::size_t>(g.k + 1)) {
      counts[toString(kmer)] = static_cast<int>(g.edge_count(e));
    }
  }
  return counts;
}

inline void VerifyGraph(const std::vector<std::string> &seqs, size_t k,
                        const fs::path &bufferPath) {
  ASSERT_TRUE(k > 3);

  auto cpath = bufferPath / "coloured";
  auto rpath = bufferPath / "counted";

  fs::create_directories(cpath);
  fs::create_directories(rpath);
  
  // create and check coloured graph

  BuildGraph(seqs, k - 1, 1, cpath, true);

  ColouredGraph g1;
  ColouredGraph::FromDisk(g1, cpath);
  CheckKmers(g1);
  
  // create and check counting graph

  BuildGraph(seqs, k - 1, 1, rpath, false);

  WeightedGraph g2;
  WeightedGraph::FromDisk(g2, rpath);

  // compare graph structures

  EXPECT_EQ(g1.l.size(), g2.l.size());
  EXPECT_THAT(g1.l, Pointwise(Eq(), g2.l));

  EXPECT_EQ(g1.W.size(), g2.W.size());
  EXPECT_THAT(g1.W, Pointwise(Eq(), g2.W));

  // compare k-mers to oracle

  auto graph_kmers = ExtractKmers(g2);
  auto oracle_map = OracleKmerCount(seqs, k);

  EXPECT_EQ(graph_kmers.size(), oracle_map.size());

  for (const auto &[kmer, expected_count] : oracle_map) {
    ASSERT_TRUE(graph_kmers.count(kmer)) << "Missing k-mer: " << kmer;
    EXPECT_EQ(graph_kmers[kmer], expected_count)
        << "Mismatch for k-mer: " << kmer;
  }
}
