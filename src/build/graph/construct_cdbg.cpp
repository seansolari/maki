
#include "maki/build/graph/construct_cdbg.hpp"
#include "maki/build/graph/build_colours.hpp"
#include "maki/build/graph/wm/wm_construct.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/build/kmers/construct_terminals.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/utils/logging.hpp"
#include <cereal/archives/binary.hpp>
#include <filesystem>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

namespace cdbg {

// -----------------------------------------------------------------------------
// Intermediate data
// -----------------------------------------------------------------------------

Buffers::Buffers(std::size_t length, std::size_t width, std::size_t k,
                 std::size_t keff)
    : kmers(length, width, k, keff), temp(length, width, k, keff), terminals(),
      b(), t() {}

void Buffers::collectKmers(const std::vector<const SequenceContainer *> &seqs_,
                           const std::vector<SuffixTable> &blocks_,
                           ShortSuffix s_) {
  std::size_t size = blocks_.back()[s_];

  LOG_DEBUG() << "Collecting k-mers for suffix " << s_.toString()
              << " ; count = " << size;

  kmers.resize(size);
  temp.resize(size);

  kmers.fill(seqs_, blocks_, s_);

  if (size > 1) {
    kmers.sort(&temp);
  }

  adjacentDifference(kmers, b);
}

void Buffers::setTerminals(TerminalRange &&data_) {
  const std::size_t n = data_.size();

  if (n == 0) {
    LOG_DEBUG() << "No terminal nodes for this suffix";
  }

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
                                             value_size(cmap->colourWidth()), k,
                                             k - s)),
      pool_() {

  LOG_INFO() << "Initialised suffix-wise k-mer extractor"
             << " (k=" << k << ", suffix=" << s << ")";

  LOG_INFO() << "Suffix plan maximum value = " << blocks_->back().maxValue();

  LOG_INFO() << "Terminal ranges total = " << terminals.size();
}

void SuffixwiseKmers::setPool(std::shared_ptr<BundlePool> &p) { pool_ = p; }

void SuffixwiseKmers::_count(push_summary &tkn) const {
  for (std::size_t i = 0; i < 5; ++i)
    std::atomic_ref(str_->F[i]) += tkn.F[i];
  for (std::size_t i = 0; i < 5; ++i)
    std::atomic_ref(str_->C[i]) += tkn.C[i];
}

std::unique_ptr<Bundle> SuffixwiseKmers::_getbundle(uint64_t id) const {
  auto bnd = pool_->acquire();
  bnd->id = id;
  return bnd;
}

std::unique_ptr<Bundle> SuffixwiseKmers::operator()(uint64_t idx) const {
  auto sfx = ShortSuffix::fromIndex(idx, s_);

  LOG_DEBUG() << "Processing suffix index " << idx << " (size=" << sfx.size()
              << ")";

  if (sfx.size() < s_) {
    return _extractPartialKmers(idx, sfx);
  } else {
    return _extractKmers(idx, sfx);
  }
}

std::unique_ptr<Bundle>
SuffixwiseKmers::_extractPartialKmers(uint64_t idx, ShortSuffix sfx) const {

  auto bfr = buffers_->obtain();
  bfr->setTerminals(terminals_.retrieve(sfx));

  const std::size_t n = bfr->terminals.size();

  LOG_DEBUG() << "Extracted " << n << " terminal (partial) k-mers for suffix "
              << sfx.toString();

  auto bnd = _getbundle(idx);
  auto &[edges, succ, carch] = bnd->payloads;

  auto counts = pushRange(bfr->terminals, bfr->t.begin(), edges, succ,
                          sfx.msb(), carch, *cmap_);

  _count(counts);

  buffers_->release(bfr);
  return bnd;
}

std::unique_ptr<Bundle> SuffixwiseKmers::_extractKmers(uint64_t idx,
                                                       ShortSuffix sfx) const {

  auto bfr = buffers_->obtain();
  bfr->collectKmers(seqs_, *blocks_, sfx);
  bfr->setTerminals(terminals_.endsWith(sfx));

  LOG_DEBUG() << "Interleaving " << bfr->kmers.size() << " k-mers and "
              << bfr->terminals.size() << " terminals";

  auto bnd = _getbundle(idx);
  auto &[edges, succ, carch] = bnd->payloads;

  auto counts =
      interleave(bfr->kmers, bfr->b.begin(), bfr->terminals, bfr->t.begin(),
                 edges, succ, sfx.msb(), carch, *cmap_);

  _count(counts);

  buffers_->release(bfr);
  return bnd;
}

// -----------------------------------------------------------------------------
// Finalisation
// -----------------------------------------------------------------------------

void initSuccSupport(ColouredGraphFiles &outp) {
  LOG_INFO() << "Initialising rank/select support for successor bitvector";

  sdsl::bit_vector succ;
  sdsl::load_from_file(succ, outp.l);

  LOG_INFO() << "Successor bitvector length = " << succ.size();

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

ColouredGraphFiles finalise(TempBuffers inp, std::size_t k, MetaColours &&cols,
                            const std::string &out) {

  LOG_INFO() << "Finalising coloured de Bruijn graph (k=" << k << ")";

  ColouredGraphFiles outp = ColouredGraph::GraphFiles(out);

  LOG_INFO() << "Constructing edge wavelet matrix";
  initW(inp.files.edges, outp.W, out);

  LOG_INFO() << "Removing temporary edge file: " << inp.files.edges;
  std::filesystem::remove(inp.files.edges);

  LOG_INFO() << "Moving successor array to final location";
  std::filesystem::rename(inp.files.succ, outp.l);

  initSuccSupport(outp);

  LOG_INFO() << "Moving colour archive";
  std::filesystem::rename(inp.files.colours, outp.archive);

  ColouredGraph g;
  g.k = k;
  g.cmap = toRegistry(std::move(cols));

  g.F[0] = 0u;
  g.F[1] = inp.str.F[0];
  g.F[2] = inp.str.F[0] + inp.str.F[1];
  g.F[3] = inp.str.F[0] + inp.str.F[1] + inp.str.F[2];
  g.F[4] = inp.str.F[0] + inp.str.F[1] + inp.str.F[2] + inp.str.F[3];

  g.C[0] = 0u;
  g.C[1] = inp.str.C[0];
  g.C[2] = inp.str.C[0] + inp.str.C[1];
  g.C[3] = inp.str.C[0] + inp.str.C[1] + inp.str.C[2];
  g.C[4] = inp.str.C[0] + inp.str.C[1] + inp.str.C[2] + inp.str.C[3];

  LOG_INFO() << "Serialising graph metadata";
  std::ofstream os(outp.meta, std::ios::binary);
  if (!os) {
    LOG_ERROR() << "Failed to open metadata file: " << outp.meta;
  }

  cereal::BinaryOutputArchive oarchive(os);
  oarchive(g);

  LOG_INFO() << "Graph construction complete";
  return outp;
}

// -----------------------------------------------------------------------------
// API
// -----------------------------------------------------------------------------

ColouredGraphFiles construct(const std::vector<const SequenceContainer *> &data,
                             MetaColours &&cmap, dbg::BuildOptions params) {

  LOG_INFO() << "Starting CDBG construction";
  LOG_INFO() << "Input sequences: " << data.size();
  LOG_INFO() << "Output directory: " << params.out;

  TempBuffers outp{.files = {.edges = params.out / "temp-edges.sdsl",
                             .succ = params.out / "temp-succ.sdsl",
                             .colours = params.out / "temp-colours.maki"},
                   .str = {}};

  LOG_INFO() << "Extracting terminal k-mers";
  auto terminals = extractTerminalsSparse(data, params.kmer_size);
  LOG_INFO() << "Total terminal entries = " << terminals.size();

  Multi sinks{EdgeSink(outp.files.edges), SuccSink(outp.files.succ),
              ColourSink(outp.files.colours)};

  ProcessChunks(SuffixwiseKmers(data, terminals.asRange(), params.kmer_size,
                                params.suffix_size, &cmap, &outp.str),
                sinks, params.pool_size, params.reserve_per_chunk,
                params.chunks());

  LOG_INFO() << "Finalising sinks";
  sinks.finalize();

  return finalise(outp, params.kmer_size, std::move(cmap), params.out);
}

} // namespace cdbg
