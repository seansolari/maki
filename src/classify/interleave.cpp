#include "maki/classify/interleave.hpp"
#include "maki/classify/utils/locks.hpp"
#include "maki/core/graph/cdbg.hpp"

#include <cstdint>

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_invoke.h>

void detail::Overlap::update() {
  oneapi::tbb::parallel_invoke(
      [&] { sdsl::util::set_to_value(B, 0); },
      [&] { sdsl::util::init_support(BpRank, &Bp); },
      [&] { sdsl::util::init_support(BpSelect, &Bp); });
}

detail::Overlap::Overlap(std::size_t size_)
    : B(size_, 0), Bp(size_, 0), BpRank(), BpSelect() {}

// ELM merge base class
// ========================================================

ELMMergeBase::ELMMergeBase(const DeBruijnGraph *qry_, const DeBruijnGraph *ref_,
                           std::size_t grainsize_)
    : qry(qry_), ref(ref_),
      opts({.grainSize = grainsize_, .querySize = qry->nodes()}), data(),
      locks(), h(0) {}

void ELMMergeBase::Interleave() {
  while (h < ref->k)
    this->DoOneInterleave();
}

std::vector<int64_t> ELMMergeBase::Fetch() {
  std::vector<int64_t> result(qry->W.size(), -1);
  oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<size_t>(0, opts.querySize, opts.grainSize),
      [&](oneapi::tbb::blocked_range<size_t> const &r) -> void {
        this->FetchRange(r.begin(), r.end(), result);
      });
  return result;
}

void ELMMergeBase::DoOneInterleave() {
  oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<int64_t>(0, opts.querySize, opts.grainSize),
      [&](const oneapi::tbb::blocked_range<int64_t> &r) -> void {
        this->InterleaveRange(r.begin(), r.end());
      });
  this->NextIteration();
}

int64_t ELMMergeBase::refFindBetween(int64_t l, int64_t r, uint8_t c) const {
  auto lrank = ref->W.rank(l, c), rrank = ref->W.rank(r, c);
  if (lrank < rrank) {
    return ref->W.select(lrank + 1, c);
  } else {
    c ^= 0b1000;
    lrank = ref->W.rank(l, c);
    rrank = ref->W.rank(r, c);
    return (lrank < rrank) ? ref->W.select(lrank + 1, c) : -1;
  }
}

// Small-query interleaving
// ========================================================

ELMMergeSmall::ELMMergeSmall(const DeBruijnGraph *qry_,
                             const DeBruijnGraph *ref_, std::size_t grainsize_)
    : ELMMergeBase(qry_, ref_, grainsize_), qo(opts.querySize),
      ro(opts.querySize), refBPos(opts.querySize), refBpPos(opts.querySize),
      rC({/* <N */ static_cast<int64_t>(ref->C[0]),
          /* <A */ static_cast<int64_t>(ref->C[1]),
          /* <C */ static_cast<int64_t>(ref->C[2]),
          /* <G */ static_cast<int64_t>(ref->C[3]),
          /* <T */ static_cast<int64_t>(ref->C[4]),
          /* <. */ static_cast<int64_t>(ref->nodes())}) {
  assert(qry->k == ref->k);
  data[0].resize(opts.querySize);
  data[1].resize(opts.querySize);
  locks.realloc(opts.querySize, opts.lockSampleRate);

  // h=1 interleaving
  std::array<std::size_t, 6> qC = {/* <N */ qry->C[0],
                                   /* <A */ qry->C[1],
                                   /* <C */ qry->C[2],
                                   /* <G */ qry->C[3],
                                   /* <T */ qry->C[4],
                                   /* <. */ opts.querySize};

  auto &Z = current(), &Zp = previous();
  std::size_t z = 0, zend;
  for (std::size_t c = 0; c < 5; ++c) {
    std::size_t qryNodes = qC[c + 1] - qC[c], refNodes = rC[c + 1] - rC[c];
    if (refNodes) {
      ro.B[z] = 1;
      refBPos[z] = rC[c];
    } else if (qryNodes) {
      qo.B[z] = 1;
    }
    zend = z + qryNodes;
    while (z < zend) {
      if (c == 0)
        Zp[z] = rC[c + 1];
      Z[z++] = rC[c + 1];
    }
  }

  // Update rank-select structures
  NextIteration();
}

int64_t ELMMergeSmall::getPrevQryBlock(int64_t zPos) const {
  auto zRank = qo.BpRank(zPos);
  if (zRank)
    return qo.BpSelect(zRank);
  else
    return 0;
}

int64_t ELMMergeSmall::getPrevRefBlock(int64_t zPos) const {
  auto zRank = ro.BpRank(zPos);
  if (zRank)
    return ro.BpSelect(zRank);
  else
    return 0;
}

detail::Input ELMMergeSmall::From(int64_t eq_) {
  return detail::Input{
      0 /* dummy value */,
      static_cast<int64_t>(qry->C[1] + qry->W.rank(eq_, 0b1001)),
      static_cast<int64_t>(qry->C[2] + qry->W.rank(eq_, 0b1010)),
      static_cast<int64_t>(qry->C[3] + qry->W.rank(eq_, 0b1011)),
      static_cast<int64_t>(qry->C[4] + qry->W.rank(eq_, 0b1100))};
}

void ELMMergeSmall::InitPrevBlocks(int64_t zqry, int64_t &prevBlock,
                                   int64_t &prevRefBlockStart,
                                   detail::Blocks &blocks, detail::Output &O) {
  int64_t prevQryBlock = getPrevQryBlock(zqry + 1),
          prevRefBlock = getPrevRefBlock(zqry + 1);

  if ((prevQryBlock == zqry) || (prevRefBlock == zqry))
    return;
  else if (prevRefBlock > prevQryBlock) { // in block started by reference
    prevBlock = prevRefBlock;
    prevRefBlockStart = refEdge(refBpPos[prevBlock]);
  } else { // in block started by query
    prevBlock = prevQryBlock;
    prevRefBlockStart = -1;
  }

  auto eref = refEdge(previous()[zqry]);
  for (uint64_t c = 1; c < 5; ++c) {
    auto cRankStart = qry->W.rank(qryEdge(prevBlock), c | 0b1000),
         cRankEnd = qry->W.rank(qryEdge(zqry), c | 0b1000);
    if (cRankStart < cRankEnd) {
      blocks[c] = prevBlock;
      O[c] = rC[c] + ref->W.rank(eref, c | 0b1000);
    }
  }
}

void ELMMergeSmall::InterleaveRange(int64_t zq_, int64_t end_) {
  auto &Zp = previous(), &Z = current();
  int64_t eq = qryEdge(zq_) /* current edge in query graph */,
          block = -1 /* current block */, refBlockStart = -1;
  detail::Input I =
      From(eq); // output position for each character in the interleaving
  detail::Blocks BlockId = {
      -1, -1, -1, -1, -1}; // block from which each character was last written
  detail::Output
      O; // cache current position for each character in reference graph
  InitPrevBlocks(zq_, block, refBlockStart, BlockId, O);

  LockedRegionManager::Accessor a(locks);
  while (zq_ < end_) {
    if (qo.Bp[zq_]) {
      block = zq_;
      refBlockStart = -1;
    } else if (ro.Bp[zq_]) {
      block = zq_;
      refBlockStart = refEdge(refBpPos[zq_]);
    }
    auto er = refEdge(Zp[zq_++]);
    do {
      uint8_t edge = qry->W[eq], c = edge & 0b0111;
      if ((edge & 0b1000) && c) {
        int64_t outz = I[c]++;
        a.access(outz);
        if (BlockId[c] != block) {
          BlockId[c] = block;
          auto cRankAtEnd = ref->W.rank(er, c | 0b1000);
          if (auto cRankAtStart = (refBlockStart == -1)
                                      ? cRankAtEnd
                                      : ref->W.rank(refBlockStart, c | 0b1000);
              cRankAtStart < cRankAtEnd) {
            ro.B[outz] = 1;
            refBPos[outz] = rC[c] + cRankAtStart;
          } else {
            qo.B[outz] = 1;
          }
          O[c] = rC[c] + cRankAtEnd;
        }
        Z[outz] = O[c];
      }
    } while (qry->l[eq++] == 0);
  }
}

void ELMMergeSmall::NextIteration() {
  oneapi::tbb::parallel_invoke([&] { qo.Bp |= qo.B; }, [&] { ro.Bp |= ro.B; },
                               [&] {
                                 oneapi::tbb::parallel_for(
                                     (size_t)0, opts.querySize,
                                     [&](const size_t &i) -> void {
                                       int64_t newv = refBPos[i];
                                       if (newv) {
                                         refBpPos[i] = newv;
                                         refBPos[i] = 0;
                                       }
                                     });
                               });
  oneapi::tbb::parallel_invoke([&] { qo.update(); }, [&] { ro.update(); });
  ++h;
}

void ELMMergeSmall::FetchRange(int64_t zq_, int64_t end_,
                               std::vector<int64_t> &out) {
  auto &Z = previous();
  int64_t q = qryEdge(zq_), qend;
  while (zq_ < end_) {
    qend = qryEdge(zq_ + 1);
    if (qo.Bp[zq_] == 0) {
      int64_t r = refEdge(Z[zq_] - 1), rend = refEdge(Z[zq_]);
      while (q < qend) {
        uint8_t c = qry->W[q];
        out[q] = refFindBetween(r, rend, c);
        ++q;
      }
    } else {
      q = qend;
    }
    ++zq_;
  }
}

std::vector<int64_t> ClassifySmall(const DeBruijnGraph *qry,
                                   const DeBruijnGraph *ref,
                                   std::size_t grainsize) {
  ELMMergeSmall Buffers(qry, ref, grainsize);
  Buffers.Interleave();
  return Buffers.Fetch();
}

// Large-query interleaving
// ========================================================

void ELMMergeLarge::InterleaveRange(int64_t, int64_t) {}

void ELMMergeLarge::NextIteration() {}

void ELMMergeLarge::FetchRange(int64_t, int64_t, std::vector<int64_t> &) {}

std::vector<int64_t> ClassifyLarge(const DeBruijnGraph *qry,
                                   const DeBruijnGraph *ref,
                                   std::size_t grainsize) {
  ELMMergeLarge Buffers(qry, ref, grainsize);
  Buffers.Interleave();
  return Buffers.Fetch();
}
