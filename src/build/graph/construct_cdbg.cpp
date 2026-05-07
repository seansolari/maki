
#include "maki/build/graph/construct_cdbg.hpp"
#include "maki/build/graph/build_colours.hpp"
#include "maki/build/graph/construct_common.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/build/kmers/construct_terminals.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/utils/logging.hpp"
#include <filesystem>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

namespace cdbg {

// -----------------------------------------------------------------------------
// Pipeline
// -----------------------------------------------------------------------------

SuffixwiseKmers::SuffixwiseKmers(
    const std::vector<const SequenceContainer *> &seqs,
    const TerminalRange &terms,
    std::shared_ptr<std::vector<SuffixTable>> &&suffixPlan,
    std::shared_ptr<dbg::BufferMaker<KmerBuffer>> &&buffers, std::size_t s,
    MetaColours *cmap, push_summary *str)
    : dbg::Suffixwise<KmerBuffer, CDBG_SINK_SET>(
          seqs, terms, std::move(suffixPlan), std::move(buffers), s, str),
      colourMap(cmap) {}

SuffixwiseKmers SuffixwiseKmers::FromSequences(
    const std::vector<const SequenceContainer *> &seqs,
    const TerminalRange &terms, std::size_t k, std::size_t s, MetaColours *cmap,
    push_summary *str) {
  auto suffixPlan =
      std::make_shared<std::vector<SuffixTable>>(createSuffixPlan(seqs, k, s));
  auto bufferFactory = std::make_shared<dbg::BufferMaker<KmerBuffer>>(
      suffixPlan->back().maxValue(), value_size(cmap->colourWidth()), k, k - s);
  return SuffixwiseKmers(seqs, terms, std::move(suffixPlan),
                         std::move(bufferFactory), s, cmap, str);
}

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


ColouredGraphFiles finalise(TempBuffers inp, std::size_t k, MetaColours &&cols,
                            const std::string &out) {

  LOG_INFO() << "Finalising coloured de Bruijn graph (k=" << k << ")";

  ColouredGraphFiles outp(out);
  dbg::detail::finaliseGraphBuffers(inp.files.edges, inp.files.succ, outp);

  LOG_INFO() << "Moving colour archive";
  std::filesystem::rename(inp.files.colours, outp.archive);

  ColouredGraph g;
  g.k = k;
  dbg::detail::finaliseGraphStructure(g, inp.str);
  g.cmap = toRegistry(std::move(cols));
  
  dbg::detail::serialize(g, outp.meta);

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

  ProcessChunks(SuffixwiseKmers::FromSequences(
                    data, terminals.asRange(), params.kmer_size,
                    params.suffix_size, &cmap, &outp.str),
                sinks, params.pool_size, params.reserve_per_chunk,
                params.chunks());

  LOG_INFO() << "Finalising sinks";
  sinks.finalize();

  return finalise(outp, params.kmer_size, std::move(cmap), params.out);
}

} // namespace cdbg
