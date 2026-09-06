#include <algorithm>
#include <string>
#include <unordered_set>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "maki/build/graph/construct_cdbg.hpp"
#include "maki/build/io/fasta.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/utils/logging.hpp"

using ::testing::ContainerEq;
using ::testing::ElementsAreArray;

// ==== Utilities to support testing
// =================================================

/**
 * Reverse complement of sequence.
 */

inline std::string ReverseComplement(const std::string &seq) {
  std::string rcomp = seq;
  // Reverse the sequence
  std::reverse(rcomp.begin(), rcomp.end());
  // Complement the bases
  for (char &nt : rcomp) {
    switch (nt) {
    case 'A':
      nt = 'T';
      break;
    case 'T':
      nt = 'A';
      break;
    case 'C':
      nt = 'G';
      break;
    case 'G':
      nt = 'C';
      break;
    }
  }
  return rcomp;
}

/**
 * Parse FASTA text into a vector of sequences. Lines starting with '>' are
 * headers. Sequences can span multiple lines; boundaries DO NOT cross (each
 * sequence is independent).
 */
inline static std::vector<std::string>
ParseFastaSequences(const std::string &fasta) {
  std::vector<std::string> sequences;
  std::string current;
  std::string line;
  for (size_t i = 0, n = fasta.size(); i <= n; ++i) {
    char c = (i < n ? fasta[i] : '\n');
    if (c == '\r')
      continue; // normalize CRLF
    if (c != '\n') {
      line.push_back(c);
    } else {
      if (!line.empty()) {
        if (!line.empty() && line[0] == '>') {
          if (!current.empty()) {
            sequences.push_back(current);
            current.clear();
          }
          // header line; ignore
        } else {
          // append sequence line (strip spaces, to upper)
          for (char ch : line) {
            if (!std::isspace(static_cast<unsigned char>(ch))) {
              current.push_back(static_cast<char>(
                  std::toupper(static_cast<unsigned char>(ch))));
            }
          }
        }
        line.clear();
      } else {
        // blank line: ignore
      }
    }
  }
  if (!current.empty()) {
    sequences.push_back(current);
  }
  return sequences;
}

/**
 * Compute all UNIQUE k-mers for a single sequence (does not cross boundaries).
 */
inline static std::unordered_set<std::string>
UniqueKMersForSequence(const std::string &seq, size_t k) {
  std::unordered_set<std::string> s;
  if (k == 0 || seq.size() < k)
    return s;
  for (size_t i = 0; i + k <= seq.size(); ++i) {
    s.insert(seq.substr(i, k));
  }
  return s;
}

/**
 * Compute all UNIQUE k-mers from a FASTA (union across sequences; no
 * cross-boundary k-mers).
 */
inline static std::unordered_set<std::string>
UniqueKMersFromFasta(const std::string &fasta, size_t k) {
  auto seqs = ParseFastaSequences(fasta);
  std::unordered_set<std::string> all;
  for (const auto &s : seqs) {
    // Forward Sequences
    {
      auto ks = UniqueKMersForSequence(s, k);
      all.insert(ks.begin(), ks.end());
    }
    // Reverse Complement
    {
      auto rks = UniqueKMersForSequence(ReverseComplement(s), k);
      all.insert(rks.begin(), rks.end());
    }
  }
  return all;
}

/**
 * Helper to compute edge count.
 */
inline static size_t edge_count(const ColouredGraph &g) { return g.edges(); }

/**
 * Build a FASTA string from sequences with trivial generated headers.
 */
inline static std::string BuildFasta(const std::vector<std::string> &seqs) {
  std::string fasta;
  for (size_t i = 0; i < seqs.size(); ++i) {
    fasta += ">seq" + std::to_string(i) + "\n";
    fasta += seqs[i] + "\n";
  }
  return fasta;
}

/**
 * Random DNA generator (A/C/G/T) with reproducible seed.
 */
inline static std::string RandomDNA(size_t length, std::mt19937 &rng) {
  static const char bases[4] = {'A', 'C', 'G', 'T'};
  std::uniform_int_distribution<int> dist(0, 3);
  std::string s;
  s.reserve(length);
  for (size_t i = 0; i < length; ++i)
    s.push_back(bases[dist(rng)]);
  return s;
}

template <class Dist>
inline static auto RandomFasta(size_t numSeqs, std::mt19937 &rng,
                               Dist &lenDist) {
  std::vector<std::string> seqs;
  seqs.reserve(numSeqs);
  size_t minLen = SIZE_MAX;
  for (size_t s = 0; s < numSeqs; ++s) {
    size_t len = static_cast<size_t>(lenDist(rng));
    seqs.push_back(RandomDNA(len, rng));
    minLen = std::min(minLen, len);
  }
  return std::make_pair(BuildFasta(seqs), minLen);
}

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

inline void PrintKmers(const ColouredGraph &g) {
  std::cout << "k-mers in graph:\n";
  const size_t E = edge_count(g);
  for (size_t i = 0; i < E; ++i) {
    std::string km = toString(g.kmer(i));
    std::cout << km << '\n';
  }
  std::cout << "end" << std::endl;
}

// ==== Core checker that targets both correctness and bwd sanity
// ====================

/**
 * Build graph from FASTA & k, then:
 *  - Check num_edges matches expected unique k-mers.
 *  - For each edge i:
 *      - kmer(i) length == k
 *      - kmer(i) is one of expected k-mers
 *      - bwd chain never yields an ID >= num_edges
 *  - Check the set of returned kmer(i) across all edges equals the expected
 * set.
 */
inline static void BuildAndCheckGraph(const std::string &fasta, size_t k,
                                      size_t s, fs::path &bufferPath) {
  MakeGraph(fasta, k - 1u, s, bufferPath);
  ColouredGraph g;
  ColouredGraph::FromDisk(g, bufferPath);

  const size_t E = edge_count(g);
  const auto expected_set = UniqueKMersFromFasta(fasta, k);

  ASSERT_GE(E, expected_set.size()) << "Edge count < unique k-mers in FASTA";

  std::unordered_set<std::string> observed;
  observed.reserve(E);

  for (size_t i = 0; i < E; ++i) {
    std::string km = toString(g.kmer(i));
    if (km.size() < k)
      continue; // Skip terminals

    ASSERT_EQ(km.size(), k)
        << "kmer(" << i << ") has length " << km.size() << ", expected " << k;
    ASSERT_TRUE(expected_set.find(km) != expected_set.end())
        << "kmer(" << i << ")='" << km << "' not found among expected k-mers.";

    // Targeted bwd sanity check: for up to k-1 steps, ensure any returned ID is
    // < E.
    size_t curr = i;
    for (size_t step = 0; step + 1 < k; ++step) { // up to k-1 backwards steps
      std::optional<size_t> prev = g.bwd(curr);
      if (!prev.has_value())
        break; // reached a boundary; valid
      ASSERT_LT(prev.value(), E)
          << "bwd(" << curr << ") returned out-of-range edge id "
          << prev.value() << " (E=" << E << ").";
      curr = prev.value();
    }

    observed.insert(km);
  }

  // Ensure bijection: all expected k-mers appear as edge labels
  ASSERT_EQ(observed.size(), expected_set.size());
  for (const auto &km : expected_set) {
    ASSERT_TRUE(observed.find(km) != observed.end())
        << "Expected k-mer '" << km << "' not found among graph edges.";
  }
}

/**
 * Check colour buffer structure.
 */
inline static void BuildAndCheckBufferEdgePositions(const std::string &fasta,
                                                    size_t k, size_t s,
                                                    fs::path &bufferPath) {
  MakeGraph(fasta, k - 1u, s, bufferPath);
  ColouredGraph g;
  ColouredGraph::FromDisk(g, bufferPath);

  std::size_t numBuffers = g.carch->chunk_count();
  LOG_DEBUG() << "checking " << numBuffers << " colour buffers";

  std::size_t edgeCount = 0;
  for (std::size_t c = 0; c < numBuffers; ++c) {
    auto meta = g.carch->meta(c);
    ASSERT_EQ(meta.start_index, edgeCount);

    auto view = g.carch->view(c);
    ASSERT_EQ(meta.elem_count, view.size());
    edgeCount += view.size();
  }
}
