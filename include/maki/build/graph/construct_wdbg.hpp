
#pragma once

#include "construct_common.hpp"
#include "maki/build/graph/archive/counts_writer.hpp"
#include "maki/build/graph/archive/sdsl_writer.hpp"
#include "maki/build/graph/archive/sink_manager.hpp"
#include "maki/build/graph/interleave_buffers.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/core/graph/wdbg.hpp"
#include "maki/core/seq/concepts.hpp"

// -----------------------------------------------------------------------------
// Construct abundance-weighted succinct de Bruijn graph
// -----------------------------------------------------------------------------
// Data extracted from input sequences is a collection of coloured k-mers that
// are interleaved and flushed to produce:
//    1) Edge array - in-memory byte array used to construct wavelet matrix.
//    2) Last array - in-memory bit vector representing node boundaries for
//        graph traversal.
//    3) Weights array - in-memory array holding occurrence counts for each
//        k-mer in the sample graph.
// -----------------------------------------------------------------------------

namespace dbg {

template <>
struct BufferMaker<TerminalBuffer> : public Factory<Buffers<TerminalBuffer>> {
  BufferMaker(std::size_t length_, std::size_t k_, std::size_t keff_)
      : Factory<Buffers<TerminalBuffer>>(), length(length_), k(k_),
        keff(keff_) {}

  std::size_t length, k, keff;

  inline std::unique_ptr<Buffers<TerminalBuffer>> obtain() {
    return Factory<Buffers<TerminalBuffer>>::obtain(
        length, k, keff, TerminalBuffer::autofit_tag);
  }
};

} // namespace dbg

namespace wdbg {

using Buffers = dbg::Buffers<TerminalBuffer>;
using BufferMaker = dbg::BufferMaker<TerminalBuffer>;

// -----------------------------------------------------------------------------
// Output data
// -----------------------------------------------------------------------------

using EdgeSink = SdslIntVectorOnDiskSink<4>;
using SuccSink = SdslIntVectorOnDiskSink<1>;
using CountSink = CountBufferSink;

#define WDBG_SINK_SET EdgeSink, SuccSink, CountSink

using Sinks = std::tuple<WDBG_SINK_SET>;
using Bundle = ChunkBundleT<WDBG_SINK_SET>;
using BundlePool = ::BundlePool<WDBG_SINK_SET>;
using Multi = MultiSink<WDBG_SINK_SET>;

struct BufferPaths {
  std::filesystem::path edges;
  std::filesystem::path succ;
};

struct TempBuffers {
  BufferPaths files;
  push_summary str;
};

// -----------------------------------------------------------------------------
// Pipeline
// -----------------------------------------------------------------------------

template <container_span T>
struct SuffixwiseTerminals
    : public dbg::Suffixwise<T, TerminalBuffer, WDBG_SINK_SET> {

  using dbg::Suffixwise<T, TerminalBuffer, WDBG_SINK_SET>::Suffixwise;

  static SuffixwiseTerminals FromSequences(const T &seqs,
                                           const TerminalRange &terminals,
                                           std::size_t k, std::size_t s,
                                           push_summary *str) {
    auto suffixPlan = std::make_shared<std::vector<SuffixTable>>(
        createSuffixPlan(seqs, k, s));
    auto bufferFactory = std::make_shared<dbg::BufferMaker<TerminalBuffer>>(
        suffixPlan->back().maxValue(), k, k - s);
    return SuffixwiseTerminals(seqs, terminals, std::move(suffixPlan),
                               std::move(bufferFactory), s, str);
  }

  std::unique_ptr<Bundle> operator()(uint64_t idx) const {
    auto sfx = ShortSuffix::fromIndex(idx, this->suffixSize);

    LOG_DEBUG() << "Processing suffix index " << idx << " (size=" << sfx.size()
                << ")";

    if (sfx.size() < this->suffixSize) {
      return extractPartialSuffix(idx, sfx);
    } else {
      return extractSuffix(idx, sfx);
    }
  }

  std::unique_ptr<Bundle> extractSuffix(uint64_t idx, ShortSuffix sfx) const {
    auto bfr = this->kmerBuffers->obtain();

    bfr->collectKmers(this->sequences, *this->suffixCounts, sfx);

    const std::size_t n = bfr->kmers.size();

    LOG_DEBUG() << "Extracted " << n << " k-mers for suffix " << sfx.toString();

    auto bnd = this->getBundle(idx);
    auto &[edges, succ, edgeCounts] = bnd->payloads;

    auto counts = pushRange(bfr->kmers, bfr->b.begin(), edges, succ, sfx.msb(),
                            edgeCounts);

    this->pushRegionStructure(counts);

    this->kmerBuffers->release(bfr);
    return bnd;
  }

  std::unique_ptr<Bundle> extractPartialSuffix(uint64_t idx,
                                               ShortSuffix sfx) const {
    auto bfr = this->kmerBuffers->obtain();
    bfr->setTerminals(this->terminals.retrieve(sfx));

    const std::size_t n = bfr->terminals.size();

    LOG_DEBUG() << "Extracted " << n << " terminal (partial) k-mers for suffix "
                << sfx.toString();

    auto bnd = this->getBundle(idx);
    auto &[edges, succ, edgeCounts] = bnd->payloads;

    auto counts = pushRange(bfr->terminals, bfr->t.begin(), edges, succ,
                            sfx.msb(), edgeCounts);

    this->pushRegionStructure(counts);

    this->kmerBuffers->release(bfr);
    return bnd;
  }
};

// -----------------------------------------------------------------------------
// Finalisation
// -----------------------------------------------------------------------------

WeightedGraphFiles finalise(TempBuffers inp, std::size_t k,
                            CountBuffer &&counts, const std::string &out);

// -----------------------------------------------------------------------------
// API
// -----------------------------------------------------------------------------

template <container_span T>
WeightedGraphFiles construct(const T &data,
                             dbg::BuildOptions params = {}) {
  LOG_INFO() << "Starting WDBG construction, node size=" << params.kmer_size;
  LOG_INFO() << "Input chunks: " << data.size();
  
  std::filesystem::create_directories(params.out);
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

  ProcessChunks(SuffixwiseTerminals<T>::FromSequences(
                    data, terminals.asRange(), params.kmer_size,
                    params.suffix_size, &outp.str),
                sinks, params.pool_size, params.reserve_per_chunk,
                params.chunks());

  LOG_INFO() << "Finalising sinks";
  sinks.finalize();

  return finalise(outp, params.kmer_size, std::move(rawCounts), params.out);
}

} // namespace wdbg
