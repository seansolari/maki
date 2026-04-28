#include "maki/classify/interleave.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/utils/locks.hpp"
#include "maki/core/utils/logging.hpp"
#include "oneapi/tbb/parallel_invoke.h"

#include <cstdint>
#include <utility>

detail::Overlap::Overlap(std::size_t size_)
    : B(size_, 0), Bp(size_, 0), BpRank(), BpSelect() {

  LOG_DEBUG() << "Initialised overlap structure with size " << size_;
}

void detail::Overlap::update() {
  LOG_DEBUG() << "Updating overlap rank/select structures";

  oneapi::tbb::parallel_invoke(
      [&] { sdsl::util::set_to_value(B, 0); },
      [&] { sdsl::util::init_support(BpRank, &Bp); },
      [&] { sdsl::util::init_support(BpSelect, &Bp); });
}

void detail::SetToOne(sdsl::bit_vector &bv, std::size_t begin_,
                      std::size_t end_) {
  while (begin_ < end_)
    bv[begin_++] = 1;
}

// Small-query interleaving
// ========================================================

ELMMergeSmall::ELMMergeSmall(const DeBruijnGraph *qry_,
                             const DeBruijnGraph *ref_, std::size_t grainsize_)
    : ELMMergeBase(qry_, ref_, grainsize_), qo(opts.querySize),
      ro(opts.querySize), refBPos(opts.querySize), refBpPos(opts.querySize),
      rC({static_cast<int64_t>(ref->C[0]), static_cast<int64_t>(ref->C[1]),
          static_cast<int64_t>(ref->C[2]), static_cast<int64_t>(ref->C[3]),
          static_cast<int64_t>(ref->C[4]),
          static_cast<int64_t>(ref->nodes())}) {
  LOG_INFO() << "Initialising small merge";
  LOG_INFO() << "Query nodes     = " << opts.querySize;
  LOG_INFO() << "Reference nodes = " << opts.referenceSize;
  LOG_DEBUG() << "Grain size      = " << grainsize_;

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

    LOG_DEBUG() << "Symbol " << c << ": query=" << qryNodes
                << ", ref=" << refNodes;
  }

  NextIteration();
}

int64_t ELMMergeSmall::getPrevQryBlock(int64_t zPos) const {
  auto r = qo.BpRank(zPos);
  return r ? qo.BpSelect(r) : 0;
}

int64_t ELMMergeSmall::getPrevRefBlock(int64_t zPos) const {
  auto r = ro.BpRank(zPos);
  return r ? ro.BpSelect(r) : 0;
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
  LOG_DEBUG() << "ELMMergeSmall interleaving range [" << zq_ << ", " << end_
              << ")";

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
  oneapi::tbb::parallel_invoke(
      [&] { qo.Bp |= qo.B; }, [&] { ro.Bp |= ro.B; },
      [&] {
        oneapi::tbb::parallel_for((size_t)0,
                                  detail::interleavingSize<ELMMergeSmall>(opts),
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
        if (auto _re = ref->findBetween(r, rend, c); _re.has_value())
          out[q] = static_cast<int64_t>(_re.value());
        else
          out[q] = -1;
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
  LOG_INFO() << "Running small-query classification";

  ELMMergeSmall buffers(qry, ref, grainsize);
  buffers.Interleave();

  LOG_INFO() << "Small-query merge complete";
  return buffers.Fetch();
}

// Large-query interleaving
// ========================================================

ELMMergeLarge::ELMMergeLarge(const DeBruijnGraph *qry_,
                             const DeBruijnGraph *ref_, std::size_t grainsize_)
    : ELMMergeBase(qry_, ref_, grainsize_),
      b(detail::interleavingSize<ELMMergeLarge>(opts)) {
  LOG_INFO() << "Initialising large merge";
  LOG_INFO() << "Total nodes = "
             << detail::interleavingSize<ELMMergeLarge>(opts);

  // h=1 interleaving
  auto &Z = current(), &Zp = previous();
  std::size_t z = 0, zmid, zend;
  for (std::size_t c = 0; c < 5; ++c) {
    if (c < 4) {
      zmid = z + ref->C[c + 1] - ref->C[c];
      zend = zmid + qry->C[c + 1] - qry->C[c];
    } else {
      zmid = z + opts.referenceSize - ref->C[c];
      zend = zmid + opts.querySize - qry->C[c];
    }

    if (c == 0) {
      detail::SetToOne(Zp, zmid, zend);
    }
    detail::SetToOne(Z, zmid, zend);
    b.B[z] = 1;
    z = zend;

    LOG_DEBUG() << "Symbol " << c << ": block range [" << zmid << ", " << zend
                << ")";
  }

  // Update rank-select structures
  NextIteration();
}

int64_t ELMMergeLarge::getPrevBlock(int64_t z) const {
  return b.BpSelect(b.BpRank(z + 1));
}

std::pair<int64_t, int64_t> ELMMergeLarge::edgesAt(int64_t z_) const {
  std::size_t qn = ZpRank(z_), rn = z_ - qn;
  return std::make_pair(refEdge(rn), qryEdge(qn));
}

std::array<detail::GraphCursor, 2> ELMMergeLarge::UpTo(int64_t z_) const {
  auto x = edgesAt(z_);
  return {detail::GraphCursor{ref, x.first},
          detail::GraphCursor{qry, x.second}};
}

detail::Input ELMMergeLarge::From(int64_t re, int64_t qe) const {
  return detail::Input{
      0, /* dummy value */
      static_cast<int64_t>(ref->C[1] + ref->W.rank(re, 0b1001) + qry->C[1] +
                           qry->W.rank(qe, 0b1001)),
      static_cast<int64_t>(ref->C[2] + ref->W.rank(re, 0b1010) + qry->C[2] +
                           qry->W.rank(qe, 0b1010)),
      static_cast<int64_t>(ref->C[3] + ref->W.rank(re, 0b1011) + qry->C[3] +
                           qry->W.rank(qe, 0b1011)),
      static_cast<int64_t>(ref->C[4] + ref->W.rank(re, 0b1100) + qry->C[4] +
                           qry->W.rank(qe, 0b1100))};
}

int64_t ELMMergeLarge::InitPrevBlocks(int64_t z, int64_t re, int64_t qe,
                                      detail::Blocks &blocks) const {
  int64_t block = getPrevBlock(z);
  if (block == z)
    return block;
  auto [pre, pqe] = edgesAt(block);
  for (uint64_t c = 1; c < 5; ++c) {
    if (ref->findBetween(pre, re, c) || qry->findBetween(pqe, qe, c))
      blocks[c] = block;
  }
  return block;
}

void ELMMergeLarge::InterleaveRange(int64_t z_, int64_t zend_) {
  LOG_DEBUG() << "ELMMergeLarge interleaving range [" << z_ << ", " << zend_
              << ")";

  auto &Zp = previous(), &Z = current();
  std::array<detail::GraphCursor, 2> O = UpTo(z_);
  detail::Input I = From(O[0].i, O[1].i);
  detail::Blocks BlockId = {-1, -1, -1, -1, -1};
  int64_t block = InitPrevBlocks(z_, O[0].i, O[1].i, BlockId);

  LOG_DEBUG() << "Interleaving at position " << z_ << ", starting at G0{node "
              << ref->lR(O[0].i) << "[" << opts.referenceSize << "], edge "
              << O[0].i << "[" << ref->edges() << "]} "
              << "G1{node " << qry->lR(O[1].i) << "[" << opts.querySize
              << "], edge " << O[1].i << "[" << qry->edges()
              << "]}, output configured to {" << I[1] << ", " << I[2] << ", "
              << I[3] << ", " << I[4] << "}"
              << " at block=" << block << " [" << BlockId[1] << ", "
              << BlockId[2] << ", " << BlockId[3] << ", " << BlockId[4] << "];";

  LockedRegionManager::Accessor a(locks);
  while (z_ < zend_) {
    if (b.Bp[z_]) {
      block = z_;
    }
    auto j = Zp[z_++];
    auto &g = O[j];
    do {
      uint8_t edge = g.edge(), c = edge & 0b0111;
      if ((edge & 0b1000) && c) {
        int64_t outz = I[c]++;
        a.access(outz);
        Z[outz] = j;
        if (BlockId[c] != block) {
          BlockId[c] = block;
          b.B[outz] = 1;
        }
      }
    } while (g++ == 0);
  }
}

void ELMMergeLarge::NextIteration() {
  oneapi::tbb::parallel_invoke(
      [&] { b.Bp |= b.B; },
      [&] { sdsl::util::init_support(ZpRank, &current()); });
  b.update();
  ++h;
}

void ELMMergeLarge::FetchRange(int64_t z_, int64_t end_,
                               std::vector<int64_t> &out) {
  while (z_ < end_) {
    if (b.Bp[z_] == 0) {
      assert(z_ > 0);
      assert(previous()[z_-1] == 0);
      assert(previous()[z_] == 1);

      auto [r, q] = edgesAt(z_-1);
      decltype(r) rend = ref->enclose(r);
      do {
        uint8_t c = qry->W[q];
        if (auto rf = ref->findBetween(r, rend, c); rf.has_value())
          out[q] = static_cast<int64_t>(rf.value());
      } while (qry->l[q++] == 0);
    }
    ++z_;
  }
}

std::vector<int64_t> ClassifyLarge(const DeBruijnGraph *qry,
                                   const DeBruijnGraph *ref,
                                   std::size_t grainsize) {

  LOG_INFO() << "Running large-query classification";

  ELMMergeLarge buffers(qry, ref, grainsize);
  buffers.Interleave();

  LOG_INFO() << "Large-query merge complete";
  return buffers.Fetch();
}
