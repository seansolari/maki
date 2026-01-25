
#include "maki/build/graph/construct_cdbg.hpp"
#include "maki/build/graph/interleave_buffers.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/build/kmers/construct_terminals.hpp"

namespace cdbg {

// -----------------------------------------------------------------------------
// Intermediate data
// -----------------------------------------------------------------------------

void Buffers::collectKmers(const std::vector<const SequenceContainer *> &seqs_,
                           const std::vector<SuffixTable> &blocks_,
                           ShortSuffix s_) {
  std::size_t size = blocks_.back()[s_];
  kmers.resize(size);
  temp.resize(size);
  kmers.fill(seqs_, blocks_, s_);
  if (size > 1) {
    kmers.sort(&temp);
  }
  adjacentDifference(kmers, b);
}

void Buffers::setTerminals(TerminalRange &&data_) {
  terminals = std::move(data_);
  adjacentDifference(terminals, t);
}

// -----------------------------------------------------------------------------
// Output data
// -----------------------------------------------------------------------------

Sinks prepareSinks(BufferPaths &pths) {
  return std::make_tuple(EdgeSink(pths.edges), SuccSink(pths.succ),
                         ColourSink(pths.colours));
}

// -----------------------------------------------------------------------------
// Pipeline
// -----------------------------------------------------------------------------

SuffixwiseKmers::SuffixwiseKmers(
    const std::vector<const SequenceContainer *> &seqs,
    const TerminalRange &terminals, std::size_t k, std::size_t s,
    MetaColours *cmap)
    : s_(s), cmap_(cmap), seqs_(seqs), terminals_(terminals),
      blocks_(std::make_shared<std::vector<SuffixTable>>(
          createSuffixPlan(seqs, k, s))),
      buffers_(std::make_shared<BufferMaker>(blocks_->back().maxValue(),
                                             value_size(cmap->maxColourWidth()),
                                             k, k - s)),
      pool_() {}

void SuffixwiseKmers::setPool(std::shared_ptr<BundlePool> &p) { pool_ = p; }

std::unique_ptr<Bundle> SuffixwiseKmers::_getbundle(uint64_t id) const {
  auto bnd = pool_->acquire();
  bnd->id = id;
  return bnd;
}

std::unique_ptr<Bundle> SuffixwiseKmers::operator()(uint64_t rnk) const {
  auto sfx = ShortSuffix::fromRank(rnk);
  if (sfx.size() == s_) {
    return _extractPartialKmers(rnk, sfx);
  } else {
    return _extractKmers(rnk, sfx);
  }
}

std::unique_ptr<Bundle>
SuffixwiseKmers::_extractPartialKmers(uint64_t rnk, ShortSuffix sfx) const {
  // request work buffer
  auto bfr = buffers_->obtain();
  // collect terminals with this partial suffix
  bfr->setTerminals(terminals_.retrieve(sfx));
  // flush to output buffers
  auto bnd = _getbundle(rnk);
  auto &[edges, succ, carch] = bnd->payloads;
  pushRange(bfr->terminals, bfr->t.begin(), edges, succ, carch, *cmap_,
            sfx.msb());
  // clear and return buffers
  buffers_->release(bfr);
  return bnd;
}

std::unique_ptr<Bundle> SuffixwiseKmers::_extractKmers(uint64_t rnk,
                                                       ShortSuffix sfx) const {
  // request work buffer
  auto bfr = buffers_->obtain();
  // collect sequence data for suffix
  bfr->collectKmers(seqs_, *blocks_, sfx);
  bfr->setTerminals(terminals_.endsWith(sfx));
  // flush to output buffers
  auto bnd = _getbundle(rnk);
  auto &[edges, succ, carch] = bnd->payloads;
  interleave(bfr->kmers, bfr->b.begin(), bfr->terminals, bfr->t.begin(), edges,
             succ, carch, *cmap_, sfx.msb());
  // clear and return buffers
  buffers_->release(bfr);
  return bnd;
}

// -----------------------------------------------------------------------------
// API
// -----------------------------------------------------------------------------

BufferPaths construct(const std::vector<const SequenceContainer *> &data,
                      MetaColours &cmap, BuildOptions params) {
  BufferPaths outp{params.out / "edges.txt", params.out / "succ.sdsl",
                   params.out / "colours.maki"};
  // extract terminals and create suffix plan
  auto terminals = extractTerminalsSparse(data, params.kmer_size);
  // suffix-wise processing and sink chunks to disk
  ProcessChunks(SuffixwiseKmers(data, terminals.asRange(), params.kmer_size,
                                params.suffix_size, &cmap),
                prepareSinks(outp), params.pool_size, params.reserve_per_chunk,
                params.num_chunks);
  return outp;
}

} // namespace cdbg
