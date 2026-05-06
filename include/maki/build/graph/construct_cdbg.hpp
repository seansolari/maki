
#pragma once

#include "construct_common.hpp"
#include "interleave_buffers.hpp"
#include "maki/build/graph/archive/archive_writer.hpp"
#include "maki/build/graph/archive/sdsl_writer.hpp"
#include "maki/build/graph/archive/sink_manager.hpp"
#include "maki/build/graph/build_colours.hpp"
#include "maki/build/kmers/buffers/kmers.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/seq/concepts.hpp"
#include <cstddef>
#include <filesystem>
#include <memory>
#include <sdsl/int_vector.hpp>
#include <vector>

namespace cdbg {

// -----------------------------------------------------------------------------
// Construct coloured succinct de Bruijn graph
// -----------------------------------------------------------------------------
// Data extracted from input sequences is a collection of coloured k-mers that
// are interleaved and flushed to produce:
//    1) Edge array - byte array on disk, where each byte represents the 4-bit
//        encoded edge value. Used for semi-external construction of wavelet
//        matrix.
//    2) Last array - Node boundaries for graph traversal. Bit vector
//        constructed sequentially on disk using SDSL int_vec_mapper<>.
//    3) Colour archive - Holds colour information in chunked format for
//        disk-based traversal and lookup.
// -----------------------------------------------------------------------------

using Buffers = dbg::Buffers<KmerBuffer>;
using BufferMaker = dbg::BufferMaker<KmerBuffer>;

// -----------------------------------------------------------------------------
// Output data
// -----------------------------------------------------------------------------

using EdgeSink = SdslIntVectorOnDiskSink<4>;
using SuccSink = SdslIntVectorOnDiskSink<1>;
using ColourSink = ArchiveWriter<>;

#define CDBG_SINK_SET EdgeSink, SuccSink, ColourSink

using Sinks = std::tuple<CDBG_SINK_SET>;
using Bundle = ChunkBundleT<CDBG_SINK_SET>;
using BundlePool = ::BundlePool<CDBG_SINK_SET>;
using Multi = MultiSink<CDBG_SINK_SET>;

struct BufferPaths {
  std::filesystem::path edges;
  std::filesystem::path succ;
  std::filesystem::path colours;
};

struct TempBuffers {
  BufferPaths files;
  push_summary str;
};

// -----------------------------------------------------------------------------
// Pipeline
// -----------------------------------------------------------------------------

struct SuffixwiseKmers : public dbg::Suffixwise<KmerBuffer, CDBG_SINK_SET> {
  SuffixwiseKmers(const std::vector<const SequenceContainer *> &seqs,
                  const TerminalRange &terminals, std::size_t k, std::size_t s,
                  MetaColours *cmap, push_summary *);

  MetaColours *colourMap;
  std::unique_ptr<Bundle> operator()(uint64_t) const;
  std::unique_ptr<Bundle> extractKmers(uint64_t, ShortSuffix) const;
  std::unique_ptr<Bundle> extractPartialKmers(uint64_t, ShortSuffix) const;
};

// -----------------------------------------------------------------------------
// API
// -----------------------------------------------------------------------------

ColouredGraphFiles construct(const std::vector<const SequenceContainer *> &data,
                             MetaColours &&cmap, dbg::BuildOptions params = {});

} // namespace cdbg
