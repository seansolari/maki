
#include "maki/build/graph/construct_wdbg.hpp"
#include "maki/build/io/fastq.hpp"
#include "maki/build/kmers/construct_terminals.hpp"
#include "maki/core/graph/archive/counts.hpp"
#include "maki/core/utils/logging.hpp"
#include <oneapi/tbb/global_control.h>

namespace wdbg {

// -----------------------------------------------------------------------------
// Pipeline
// -----------------------------------------------------------------------------

SuffixwiseTerminals SuffixwiseTerminals::FromSequences(
    const std::vector<const SequenceContainer *> &seqs,
    const TerminalRange &terminals, std::size_t k, std::size_t s,
    push_summary *str) {
  auto suffixPlan =
      std::make_shared<std::vector<SuffixTable>>(createSuffixPlan(seqs, k, s));
  auto bufferFactory = std::make_shared<dbg::BufferMaker<TerminalBuffer>>(
      suffixPlan->back().maxValue(), k, k - s);
  return SuffixwiseTerminals(seqs, terminals, std::move(suffixPlan),
                             std::move(bufferFactory), s, str);
}

std::unique_ptr<Bundle> SuffixwiseTerminals::operator()(uint64_t idx) const {
  auto sfx = ShortSuffix::fromIndex(idx, suffixSize);

  LOG_DEBUG() << "Processing suffix index " << idx << " (size=" << sfx.size()
              << ")";

  if (sfx.size() < suffixSize) {
    return extractPartialSuffix(idx, sfx);
  } else {
    return extractSuffix(idx, sfx);
  }
}

std::unique_ptr<Bundle>
SuffixwiseTerminals::extractPartialSuffix(uint64_t idx, ShortSuffix sfx) const {
  auto bfr = kmerBuffers->obtain();
  bfr->setTerminals(terminals.retrieve(sfx));

  const std::size_t n = bfr->terminals.size();

  LOG_DEBUG() << "Extracted " << n << " terminal (partial) k-mers for suffix "
              << sfx.toString();

  auto bnd = getBundle(idx);
  auto &[edges, succ, edgeCounts] = bnd->payloads;

  auto counts = pushRange(bfr->terminals, bfr->t.begin(), edges, succ,
                          sfx.msb(), edgeCounts);

  pushRegionStructure(counts);

  kmerBuffers->release(bfr);
  return bnd;
}

std::unique_ptr<Bundle>
SuffixwiseTerminals::extractSuffix(uint64_t idx, ShortSuffix sfx) const {
  auto bfr = kmerBuffers->obtain();

  bfr->collectKmers(sequences, *suffixCounts, sfx);

  const std::size_t n = bfr->kmers.size();

  LOG_DEBUG() << "Extracted " << n << " k-mers for suffix " << sfx.toString();

  auto bnd = getBundle(idx);
  auto &[edges, succ, edgeCounts] = bnd->payloads;

  auto counts =
      pushRange(bfr->kmers, bfr->b.begin(), edges, succ, sfx.msb(), edgeCounts);

  pushRegionStructure(counts);

  kmerBuffers->release(bfr);
  return bnd;
}

// -----------------------------------------------------------------------------
// Finalisation
// -----------------------------------------------------------------------------

WeightedGraphFiles finalise(TempBuffers inp, std::size_t k,
                            CountBuffer &&counts, const std::string &out) {

  LOG_INFO() << "Finalising weighted de Bruijn graph (k=" << k << ")";

  WeightedGraphFiles outp(out);
  dbg::detail::finaliseGraphBuffers(inp.files.edges, inp.files.succ, outp);

  WeightedGraph g;
  g.k = k;
  dbg::detail::finaliseGraphStructure(g, inp.str);

  LOG_INFO() << "Compressing k-mer counts buffer";
  const std::size_t orig_counts_size = detail::size_in_bytes(counts);
  g.occ = CompressedCountBuffer(std::move(counts));
  const std::size_t new_counts_size = detail::size_in_bytes(g.occ);
  LOG_INFO() << "Compression complete, uncmp=" << orig_counts_size
             << ", cmp=" << new_counts_size;

  dbg::detail::serialize(g, outp.meta);

  LOG_INFO() << "Graph construction complete";
  return outp;
}

// -----------------------------------------------------------------------------
// API
// -----------------------------------------------------------------------------

WeightedGraphFiles construct(const std::vector<const SequenceContainer *> &data,
                             dbg::BuildOptions params) {

  LOG_INFO() << "Starting WDBG construction, node size=" << params.kmer_size;
  LOG_INFO() << "Input chunks: " << data.size();
  LOG_INFO() << "Output directory: " << params.out;

  TempBuffers outp{.files = {.edges = params.out / "temp-edges.sdsl",
                             .succ = params.out / "temp-succ.sdsl"},
                   .str = {}};

  // Graph edge counts
  CountBuffer rawCounts;

  LOG_INFO() << "Extracting terminal k-mers";
  auto terminals = extractTerminalsDense(data, params.suffix_size);
  LOG_INFO() << "Total terminal entries = " << terminals.size();

  Multi sinks{EdgeSink(outp.files.edges), SuccSink(outp.files.succ),
              CountSink(&rawCounts)};

  ProcessChunks(SuffixwiseTerminals::FromSequences(
                    data, terminals.asRange(), params.kmer_size,
                    params.suffix_size, &outp.str),
                sinks, params.pool_size, params.reserve_per_chunk,
                params.chunks());

  LOG_INFO() << "Finalising sinks";
  sinks.finalize();

  return finalise(outp, params.kmer_size, std::move(rawCounts), params.out);
}

WeightedGraphFiles construct(const DataFilePair &fp, dbg::BuildOptions params) {
  oneapi::tbb::global_control global_limit(
      oneapi::tbb::global_control::max_allowed_parallelism, params.threads);

  auto chunks =
      chunkReads(detail::parsePairedFastq(fp, params.kmer_size,
                                          params.kmer_size, params.threads),
                 10 * params.threads);
  auto view = toView(chunks);
  return construct(view, params);
}

} // namespace wdbg
