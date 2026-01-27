
#include "maki/build/graph/construct_cdbg.hpp"
#include "maki/build/graph/wm/wm_construct.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/build/kmers/construct_terminals.hpp"
#include <filesystem>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

namespace cdbg {

std::size_t BuildOptions::chunks() const {
  return ShortSuffix::numSuffixes(suffix_size);
}

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
// Pipeline
// -----------------------------------------------------------------------------

SuffixwiseKmers::SuffixwiseKmers(
    const std::vector<const SequenceContainer *> &seqs,
    const TerminalRange &terminals, std::size_t k, std::size_t s,
    MetaColours *cmap, push_summary *str)
    : s_(s), cmap_(cmap), str_(str), seqs_(seqs), terminals_(terminals),
      blocks_(std::make_shared<std::vector<SuffixTable>>(
          createSuffixPlan(seqs, k, s))),
      buffers_(std::make_shared<BufferMaker>(blocks_->back().maxValue(),
                                             value_size(cmap->colourWidth()),
                                             k, k - s)),
      pool_() {}

void SuffixwiseKmers::setPool(std::shared_ptr<BundlePool> &p) { pool_ = p; }

void SuffixwiseKmers::_count(push_summary &tkn) const {
  for (std::size_t i = 0; i < 5; ++i)
    std::atomic_ref{str_->F[i]} += tkn.F[i];
  for (std::size_t i = 0; i < 5; ++i)
    std::atomic_ref{str_->C[i]} += tkn.C[i];
}

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
  auto counts = pushRange(bfr->terminals, bfr->t.begin(), edges, succ, carch,
                          *cmap_, sfx.msb());
  _count(counts);
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
  auto counts =
      interleave(bfr->kmers, bfr->b.begin(), bfr->terminals, bfr->t.begin(),
                 edges, succ, carch, *cmap_, sfx.msb());
  _count(counts);
  // clear and return buffers
  buffers_->release(bfr);
  return bnd;
}

// -----------------------------------------------------------------------------
// Finalisation
// -----------------------------------------------------------------------------

void initSuccSupport(ColouredGraphFiles &outp) {
  sdsl::bit_vector succ;
  sdsl::load_from_file(succ, outp.l);
  {
    sdsl::rank_support_v5<1, 1> lRnk;
    sdsl::util::init_support(lRnk, &succ);
    sdsl::store_to_file(std::move(lRnk), outp.lR);
  }
  {
    sdsl::select_support_mcl<1, 1> lSel;
    sdsl::util::init_support(lSel, &succ);
    sdsl::store_to_file(std::move(lSel), outp.lS);
  }
}

ColouredGraphFiles finalise(TempBuffers inp, const std::string &out) {
  ColouredGraphFiles outp = graphFiles(out);

  // edges wavelet matrix
  initW(inp.files.edges, outp.W, out);
  std::filesystem::remove(inp.files.edges);

  // succ array
  std::filesystem::rename(inp.files.succ, outp.l);
  initSuccSupport(outp);

  // colour archive
  std::filesystem::rename(inp.files.colours, outp.archive);

  return outp;
}

// -----------------------------------------------------------------------------
// API
// -----------------------------------------------------------------------------

ColouredGraphFiles construct(const std::vector<const SequenceContainer *> &data,
                             MetaColours &cmap, BuildOptions params) {
  TempBuffers outp{.files = {.edges = params.out / "temp-edges.sdsl",
                             .succ = params.out / "temp-succ.sdsl",
                             .colours = params.out / "temp-colours.maki"},
                   .str = {}};
  // extract terminals and create suffix plan
  auto terminals = extractTerminalsSparse(data, params.kmer_size);
  // suffix-wise processing and sink chunks to disk
  Multi sinks{EdgeSink(outp.files.edges), SuccSink(outp.files.succ),
              ColourSink(outp.files.colours)};
  ProcessChunks(SuffixwiseKmers(data, terminals.asRange(), params.kmer_size,
                                params.suffix_size, &cmap, &outp.str),
                sinks, params.pool_size, params.reserve_per_chunk,
                params.chunks());
  sinks.finalize();
  return finalise(outp, params.out);
}

} // namespace cdbg
