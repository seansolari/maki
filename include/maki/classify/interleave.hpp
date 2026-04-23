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

namespace detail {

using Input = std::array<int64_t, 5>;
using Blocks = std::array<int64_t, 5>;
using Output = std::array<int64_t, 5>;

struct Overlap {
  using bit_vector = sdsl::bit_vector;
  bit_vector B, Bp;
  bit_vector::rank_1_type BpRank;
  bit_vector::select_1_type BpSelect;
  Overlap(std::size_t);
  void update();
};

}; // namespace detail

// ELM merge base class
// ========================================================

class ELMMergeBase {
public:
  ELMMergeBase(const DeBruijnGraph *qry_, const DeBruijnGraph *ref_,
               std::size_t grainsize_);

  void Interleave();
  std::vector<int64_t> Fetch();

  virtual ~ELMMergeBase() = default;

protected:
  const DeBruijnGraph *qry, *ref;
  const InterleavingOpts opts;

  std::array<std::vector<int64_t>, 2> data;
  LockedRegionManager locks;
  uint8_t h;

protected:
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

  // find edge encoding character `c` between edges `[l, r)`
  int64_t refFindBetween(int64_t l, int64_t r, uint8_t c) const;

  // call `InterleaveRange()` on each range and then call `NextIteration()`
  void DoOneInterleave();

protected:
  virtual void InterleaveRange(int64_t, int64_t) = 0;
  virtual void NextIteration() = 0;
  virtual void FetchRange(int64_t, int64_t, std::vector<int64_t> &) = 0;
};

// Small-query interleaving
// ========================================================

class ELMMergeSmall : public ELMMergeBase {
public:
  ELMMergeSmall(const DeBruijnGraph *qry_, const DeBruijnGraph *ref_,
                std::size_t grainsize_);

private:
  detail::Overlap qo, ro;
  std::vector<int64_t> refBPos, refBpPos;
  std::array<int64_t, 6> rC;

  int64_t getPrevQryBlock(int64_t zPos) const;
  int64_t getPrevRefBlock(int64_t zPos) const;

  inline auto &refOverlap() const { return ro.Bp; }
  inline auto &qryOverlap() const { return qo.Bp; }
  inline auto &refBlocks() const { return refBpPos; }

  detail::Input From(int64_t);

  void InitPrevBlocks(int64_t, int64_t &, int64_t &, detail::Blocks &,
                      detail::Output &);

protected:
  virtual void InterleaveRange(int64_t, int64_t) override final;
  virtual void NextIteration() override final;
  virtual void FetchRange(int64_t, int64_t,
                          std::vector<int64_t> &) override final;
};

std::vector<int64_t> ClassifySmall(const DeBruijnGraph *qry,
                                   const DeBruijnGraph *ref,
                                   std::size_t grainsize);

// Large-query interleaving
// ========================================================

class ELMMergeLarge : public ELMMergeBase {
public:
  ELMMergeLarge(const DeBruijnGraph *qry_, const DeBruijnGraph *ref_,
                std::size_t grainsize_);

private:
  detail::Overlap b;

protected:
  virtual void InterleaveRange(int64_t, int64_t) override final;
  virtual void NextIteration() override final;
  virtual void FetchRange(int64_t, int64_t,
                          std::vector<int64_t> &) override final;
};

std::vector<int64_t> ClassifyLarge(const DeBruijnGraph *qry,
                                   const DeBruijnGraph *ref,
                                   std::size_t grainsize);
