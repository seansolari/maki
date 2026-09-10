
#pragma once

#include "construct_common.hpp"
#include "maki/build/graph/archive/counts_writer.hpp"
#include "maki/build/graph/archive/sdsl_writer.hpp"
#include "maki/build/graph/archive/sink_manager.hpp"
#include "maki/build/graph/interleave_buffers.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/core/graph/wdbg.hpp"
#include "maki/core/seq/concepts.hpp"
#include "maki/core/seq/io.hpp"


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

struct SuffixwiseTerminals
    : public dbg::Suffixwise<TerminalBuffer, WDBG_SINK_SET> {

  using dbg::Suffixwise<TerminalBuffer, WDBG_SINK_SET>::Suffixwise;

  static SuffixwiseTerminals
  FromSequences(const std::vector<const SequenceContainer *> &seqs,
                const TerminalRange &terminals, std::size_t k, std::size_t s,
                push_summary *);

  std::unique_ptr<Bundle> operator()(uint64_t) const;
  std::unique_ptr<Bundle> extractPartialSuffix(uint64_t, ShortSuffix) const;
  std::unique_ptr<Bundle> extractSuffix(uint64_t, ShortSuffix) const;
};

// -----------------------------------------------------------------------------
// API
// -----------------------------------------------------------------------------

WeightedGraphFiles construct(const std::vector<const SequenceContainer *> &data,
                             dbg::BuildOptions params = {});

} // namespace wdbg
