#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support_mcl.hpp>

#include "maki/classify/utils/locks.hpp"
#include "maki/core/graph/cdbg.hpp"

struct InterleavingOpts {
  std::size_t grainSize, querySize;
  std::size_t lockSampleRate = 8 * 1024;
};

using Input = std::array<int64_t, 5>;
using Blocks = std::array<int64_t, 5>;
using Output = std::array<int64_t, 5>;

// Small-query interleaving

class ELMMergeSmall {
public:
  ELMMergeSmall(const DeBruijnGraph *qry_, const DeBruijnGraph *ref_,
                std::size_t grainsize_);

  void Interleave();
  std::vector<int64_t> Fetch();

private:
  const DeBruijnGraph *qry, *ref;
  const InterleavingOpts opts;

  std::array<std::vector<int64_t>, 2> data;
  struct Overlap {
    using bit_vector = sdsl::bit_vector;
    bit_vector B, Bp;
    bit_vector::rank_1_type BpRank;
    bit_vector::select_1_type BpSelect;
    Overlap(std::size_t);
    void update();
  } qo, ro;
  std::vector<int64_t> refBPos, refBpPos;
  std::array<int64_t, 6> rC;
  LockedRegionManager locks;

  uint8_t h;

private:
  int64_t getPrevQryBlock(int64_t zPos) const;
  int64_t getPrevRefBlock(int64_t zPos) const;
  int64_t refFindBetween(int64_t l, int64_t r, uint8_t c) const;

  // get index corresponding to first edge of node at index `nodeIndex` in query
  // graph
  inline int64_t qryEdge(int64_t nodeIndex) {
    return (nodeIndex == 0) ? 0 : qry->lS(nodeIndex) + 1;
  }

  // get index corresponding to first edge of node at index `nodeIndex` in
  // reference graph
  inline int64_t refEdge(int64_t nodeIndex) {
    return (nodeIndex == 0) ? 0 : ref->lS(nodeIndex) + 1;
  }

  // output buffer for current interleaving iteration
  inline auto &current() { return data[h % 2]; }

  // input buffer for current interleaving iteration
  inline auto &previous() { return data[(h + 1) % 2]; }
  inline auto &refOverlap() const { return ro.Bp; }
  inline auto &qryOverlap() const { return qo.Bp; }
  inline auto &refBlocks() const { return refBpPos; }

private:
  Input From(int64_t);
  void InitPrevBlocks(int64_t, int64_t &, int64_t &, Blocks &, Output &);
  void InterleaveRange(int64_t, int64_t);
  void DoOneInterleave();
  void NextIteration();
  void FetchRange(int64_t, int64_t, std::vector<int64_t> &);
};

// Large-query interleaving



// Classify API

std::vector<int64_t> ClassifySmall(const DeBruijnGraph *qry,
                                   const DeBruijnGraph *ref,
                                   std::size_t grainsize);
