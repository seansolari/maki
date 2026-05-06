
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
// Pipeline
// -----------------------------------------------------------------------------

SuffixwiseKmers::SuffixwiseKmers(
    const std::vector<const SequenceContainer *> &seqs,
    const TerminalRange &terms, std::size_t k, std::size_t s, MetaColours *cmap,
    push_summary *str)
    : dbg::Suffixwise<KmerBuffer, CDBG_SINK_SET>(
          seqs, terms,
          std::make_shared<std::vector<SuffixTable>>(
              createSuffixPlan(seqs, k, s)),
          value_size(cmap->colourWidth()), k, s, str),
      colourMap(cmap) {}

std::unique_ptr<Bundle> SuffixwiseKmers::operator()(uint64_t idx) const {
  auto sfx = ShortSuffix::fromIndex(idx, suffixSize);

  LOG_DEBUG() << "Processing suffix index " << idx << " (size=" << sfx.size()
              << ")";

  if (sfx.size() < suffixSize) {
    return extractPartialKmers(idx, sfx);
  } else {
    return extractKmers(idx, sfx);
  }
}

std::unique_ptr<Bundle>
SuffixwiseKmers::extractPartialKmers(uint64_t idx, ShortSuffix sfx) const {

  auto bfr = kmerBuffers->obtain();
  bfr->setTerminals(terminals.retrieve(sfx));

  const std::size_t n = bfr->terminals.size();

  LOG_DEBUG() << "Extracted " << n << " terminal (partial) k-mers for suffix "
              << sfx.toString();

  auto bnd = getBundle(idx);
  auto &[edges, succ, carch] = bnd->payloads;

  auto counts = pushRange(bfr->terminals, bfr->t.begin(), edges, succ,
                          sfx.msb(), carch, *colourMap);

  pushRegionStructure(counts);

  kmerBuffers->release(bfr);
  return bnd;
}

std::unique_ptr<Bundle> SuffixwiseKmers::extractKmers(uint64_t idx,
                                                      ShortSuffix sfx) const {

  auto bfr = kmerBuffers->obtain();
  bfr->collectKmers(sequences, *suffixCounts, sfx);
  bfr->setTerminals(terminals.endsWith(sfx));

  LOG_DEBUG() << "Interleaving " << bfr->kmers.size() << " k-mers and "
              << bfr->terminals.size() << " terminals";

  auto bnd = getBundle(idx);
  auto &[edges, succ, carch] = bnd->payloads;

  auto counts =
      interleave(bfr->kmers, bfr->b.begin(), bfr->terminals, bfr->t.begin(),
                 edges, succ, sfx.msb(), carch, *colourMap);

  pushRegionStructure(counts);

  kmerBuffers->release(bfr);
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
