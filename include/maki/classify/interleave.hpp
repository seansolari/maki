#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_invoke.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support_mcl.hpp>

#include "maki/classify/utils/locks.hpp"
#include "maki/core/graph/cdbg.hpp"

struct InterleavingOpts {
  std::size_t grainSize, querySize, referenceSize;
  std::size_t lockSampleRate = 8 * 1024;
  std::size_t size() const;
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

void SetToOne(sdsl::bit_vector&, std::size_t begin_, std::size_t end_);

}; // namespace detail

// ELM merge base class
// ========================================================

template <typename T> class ELMMergeBase {
public:
  ELMMergeBase(const DeBruijnGraph *qry_, const DeBruijnGraph *ref_,
               std::size_t grainsize_)
      : qry(qry_), ref(ref_), opts({.grainSize = grainsize_,
                                    .querySize = qry->nodes(),
                                    .referenceSize = ref->nodes()}),
        data(), locks(), h(0) {
    assert(qry->k == ref->k);
  }

  void Interleave() {
    while (h < ref->k)
      this->DoOneInterleave();
  }

  std::vector<int64_t> Fetch() {
    std::vector<int64_t> result(qry->W.size(), -1);
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<size_t>(0, opts.querySize, opts.grainSize),
        [&](oneapi::tbb::blocked_range<size_t> const &r) -> void {
          this->FetchRange(r.begin(), r.end(), result);
        });
    return result;
  }

  virtual ~ELMMergeBase() = default;

protected:
  const DeBruijnGraph *qry, *ref;
  const InterleavingOpts opts;

  std::array<T, 2> data;
  LockedRegionManager locks;
  uint8_t h;

protected:
  // get index corresponding to first edge of node at index `nodeIndex` in query
  // graph
  inline int64_t qryEdge(int64_t nodeIndex) const {
    return (nodeIndex == 0) ? 0 : qry->lS(nodeIndex) + 1;
  }

  // get index corresponding to first edge of node at index `nodeIndex` in
  // reference graph
  inline int64_t refEdge(int64_t nodeIndex) const {
    return (nodeIndex == 0) ? 0 : ref->lS(nodeIndex) + 1;
  }

  // output buffer for current interleaving iteration
  inline auto &current() { return data[h % 2]; }

  // input buffer for current interleaving iteration
  inline auto &previous() { return data[(h + 1) % 2]; }

  // call `InterleaveRange()` on each range and then call `NextIteration()`
  void DoOneInterleave() {
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<int64_t>(0, opts.querySize, opts.grainSize),
        [&](const oneapi::tbb::blocked_range<int64_t> &r) -> void {
          this->InterleaveRange(r.begin(), r.end());
        });
    this->NextIteration();
  }

protected:
  virtual void InterleaveRange(int64_t, int64_t) = 0;
  virtual void NextIteration() = 0;
  virtual void FetchRange(int64_t, int64_t, std::vector<int64_t> &) = 0;
};

// Small-query interleaving
// ========================================================

class ELMMergeSmall : public ELMMergeBase<std::vector<int64_t>> {
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

class ELMMergeLarge : public ELMMergeBase<sdsl::bit_vector> {
public:
  ELMMergeLarge(const DeBruijnGraph *qry_, const DeBruijnGraph *ref_,
                std::size_t grainsize_);

private:
  sdsl::bit_vector::rank_1_type ZpRank;
  detail::Overlap b;

  // Get position in interleaving corresponding to start of current block
  // within the range `[0, z_]`.
  int64_t getPrevBlock(int64_t z_) const;

  // Get edge positions in each graph corresponding to an index of the
  // interleaving
  std::pair<int64_t, int64_t> edgesAt(int64_t) const;

  // Get pointers for current position in the interleaving
  std::array<std::pair<const DeBruijnGraph*, int64_t>, 2> UpTo(int64_t) const;

  // Count occurrences of edges up to a position in the interleaving
  detail::Input From(int64_t, int64_t) const;

  // Check pressence of edges within a range of an interleaving block
  int64_t InitPrevBlocks(int64_t, int64_t, int64_t, detail::Blocks &) const;

protected:
  virtual void InterleaveRange(int64_t, int64_t) override final;
  virtual void NextIteration() override final;
  virtual void FetchRange(int64_t, int64_t,
                          std::vector<int64_t> &) override final;
};

std::vector<int64_t> ClassifyLarge(const DeBruijnGraph *qry,
                                   const DeBruijnGraph *ref,
                                   std::size_t grainsize);
