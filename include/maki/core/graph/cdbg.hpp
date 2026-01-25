
#pragma once
#include <filesystem>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/select_support_mcl.hpp>
#include "maki/core/seq/seq_io.hpp"

namespace fs = std::filesystem;

struct ColouredGraphFiles {
  fs::path l, lR, lS, W, archive, meta;
};

ColouredGraphFiles graphFiles(fs::path base) {
  return ColouredGraphFiles{
    .l = base / "succ.dat",
    .lR = base / "succ-rank.dat",
    .lS = base / "succ-select.dat",
    .W = base / "edges.dat",
    .archive = base / "archive.dat",
    .meta = base / "graph.dat"
  };
}

struct ColouredGraph {
  const uint8_t k;
  const sdsl::bit_vector l;
  const sdsl::rank_support_v5<1, 1> lR;
  const sdsl::select_support_mcl<1, 1> lS;
  const sdsl::wt_int<> W;
  const ColourBufferRegistry archive;
  const std::array<size_t, 5> F, C;

  // number of nodes in graph
  inline size_t nodes() const { return lR(l.size()); }

  // number of edges in graph
  inline size_t edges() const { return W.size(); }

  struct IndexRange {
    IndexRange(size_t begin_, size_t end_) : begin(begin_), end(end_) {}
    size_t begin, end;
    inline size_t size() const noexcept { return end - begin; }
    inline bool empty() const noexcept { return begin == end; }
  };

  // step forward from `i`
  std::optional<size_t> fwd(size_t i) const;

  // step backward from `i`
  std::optional<size_t> bwd(size_t i) const;

  // returns edge value at position
  uint8_t edge(size_t i) const;

  // return k-mer represented by edge `i`
  Dna4Sequence kmer(size_t i) const;

  inline IndexRange getRoot() const {
    return IndexRange(0u, lastSucc(0u) + 1u);
  }

  // get Node/k-mer containing the edge at `i` in the form `[begin, end)`
  inline IndexRange getNode(size_t i) const {
    assert(i > 0u);
    return IndexRange(lastPred(i - 1) + 1u, lastSucc(i) + 1u);
  }

  // returns index of node by taking edge `e` from node at `i` (if the edge
  // exists)
  std::optional<size_t> outgoing(size_t i, uint8_t e) const;

protected:
  uint8_t block(size_t i) const;

  // Return index of first occurrence of `1` at or before `i`.
  inline size_t lastPred(size_t i) const { return lS(lR(i + 1)); }

  // Return index of first occurrence of `1` at or after `i`.
  inline size_t lastSucc(size_t i) const { return lS(lR(i) + 1); }

  // Return index of first occurrence of `c` at or before `i`.
  inline size_t WPred(size_t i, uint8_t c) const {
    return W.select(W.rank(i + 1, c), c);
  }

  // Return index of first occurrence of `c` at or after `i`.
  inline size_t WSucc(size_t i, uint8_t c) const {
    return W.select(W.rank(i, c) + 1, c);
  }
};
