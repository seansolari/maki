
#pragma once
#include "maki/core/seq/io.hpp"
#include <filesystem>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/wm_int.hpp>
#include <cereal/types/polymorphic.hpp>

namespace fs = std::filesystem;
using wavelet_matrix = sdsl::wm_int<>;

struct DeBruijnGraphFiles {
  DeBruijnGraphFiles(fs::path base_);
  fs::path base, l, lR, lS, W, meta;
};

struct DeBruijnGraph {
protected:
  static void LoadBaseBuffers(DeBruijnGraph &g, DeBruijnGraphFiles &files);

public:
  uint8_t k;

  // Structural data
  sdsl::bit_vector l;
  sdsl::rank_support_v5<1, 1> lR;
  sdsl::select_support_mcl<1, 1> lS;
  wavelet_matrix W;
  std::array<std::size_t, 5> F, C;

  // number of nodes in graph
  inline std::size_t nodes() const { return lR(l.size()); }

  // number of edges in graph
  inline std::size_t edges() const { return W.size(); }

  struct IndexRange {
    IndexRange(std::size_t begin_, std::size_t end_)
        : begin(begin_), end(end_) {}
    std::size_t begin, end;
    inline std::size_t size() const noexcept { return end - begin; }
    inline bool empty() const noexcept { return begin == end; }
  };

  // step forward from `i`
  std::optional<std::size_t> fwd(std::size_t i) const;

  // step backward from `i`
  std::optional<std::size_t> bwd(std::size_t i) const;

  // returns edge value at position
  uint8_t edge(std::size_t i) const;

  // return k-mer represented by edge `i`
  Dna4Sequence kmer(std::size_t i) const;

  // get start of node containing edge `i`
  inline std::size_t enopen(std::size_t i) const {
    return lastPred(i - 1) + 1u;
  }

  // get one-past-end of node containing edge `i`
  inline std::size_t enclose(std::size_t i) const { return lastSucc(i) + 1u; }

  inline IndexRange getRoot() const {
    return IndexRange(0u, lastSucc(0u) + 1u);
  }

  // get Node/k-mer containing the edge at `i` in the form `[begin, end)`
  inline IndexRange getNode(std::size_t i) const {
    assert(i > 0u);
    return IndexRange(enopen(i), enclose(i));
  }

  // returns index of node by taking edge `e` from node at `i` (if the edge
  // exists)
  std::optional<std::size_t> outgoing(std::size_t i, uint8_t e) const;

  // find edge encoding character `c` between edges `[l, r)`
  std::optional<std::size_t> findBetween(int64_t l, int64_t r, uint8_t c) const;

protected:
  uint8_t block(std::size_t i) const;

  // Return index of first occurrence of `1` at or before `i`.
  inline std::size_t lastPred(std::size_t i) const { return lS(lR(i + 1)); }

  // Return index of first occurrence of `1` at or after `i`.
  inline std::size_t lastSucc(std::size_t i) const { return lS(lR(i) + 1); }

  // Return index of first occurrence of `c` at or before `i`.
  inline std::size_t WPred(std::size_t i, uint8_t c) const {
    return W.select(W.rank(i + 1, c), c);
  }

  // Return index of first occurrence of `c` at or after `i`.
  inline std::size_t WSucc(std::size_t i, uint8_t c) const {
    return W.select(W.rank(i, c) + 1, c);
  }

public:
  friend class cereal::access;
  template <class Archive> void serialize(Archive &ar) { ar(k, F, C); }
};
