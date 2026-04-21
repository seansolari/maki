
#pragma once

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

struct BuildOptions {
  // Algorithm parameters
  std::size_t kmer_size = 31;
  std::size_t suffix_size = 8;
  // Output parameters
  std::filesystem::path out;
  // Space parameters
  std::size_t pool_size = 16;
  std::size_t reserve_per_chunk = 0;
  std::size_t chunks() const;
};

// -----------------------------------------------------------------------------
// Intermediate data
// -----------------------------------------------------------------------------

struct Buffers {
  Buffers(std::size_t length, std::size_t width, std::size_t k,
          std::size_t keff)
      : kmers(length, width, k, keff), temp(length, width, k, keff),
        terminals(), b(), t() {}

  /**
   * Collect k-mers with given suffix `s_` from input sequences.
   */
  void collectKmers(const std::vector<const SequenceContainer *> &seqs_,
                    const std::vector<SuffixTable> &blocks_, ShortSuffix s_);

  /**
   * Load terminals data.
   */
  void setTerminals(TerminalRange &&t_);

  KmerBuffer kmers;
  KmerBuffer temp;
  TerminalRange terminals;
  sdsl::int_vector<2> b, t;
};

struct BufferMaker : public Factory<Buffers> {
  BufferMaker(std::size_t length_, std::size_t width_, std::size_t k_,
              std::size_t keff_)
      : Factory<Buffers>(), length(length_), width(width_), k(k_), keff(keff_) {
  }

  std::size_t length;
  std::size_t width;
  std::size_t k;
  std::size_t keff;

  inline std::unique_ptr<Buffers> obtain() {
    return Factory<Buffers>::obtain(length, width, k, keff);
  }
};

// -----------------------------------------------------------------------------
// Output data
// -----------------------------------------------------------------------------

using EdgeSink = SdslIntVectorOnDiskSink<4>;
using SuccSink = SdslIntVectorOnDiskSink<1>;
using ColourSink = ArchiveWriter<>;

using Sinks = std::tuple<EdgeSink, SuccSink, ColourSink>;
using Bundle = ChunkBundleT<EdgeSink, SuccSink, ColourSink>;
using BundlePool = ::BundlePool<EdgeSink, SuccSink, ColourSink>;
using Multi = MultiSink<EdgeSink, SuccSink, ColourSink>;

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

struct SuffixwiseKmers {
  SuffixwiseKmers(const std::vector<const SequenceContainer *> &seqs,
                  const TerminalRange &terminals, std::size_t k, std::size_t s,
                  MetaColours *cmap, push_summary *);

  void setPool(std::shared_ptr<BundlePool> &p);
  std::unique_ptr<Bundle> operator()(uint64_t) const;

protected:
  void _count(push_summary &) const;
  std::unique_ptr<Bundle> _getbundle(uint64_t id) const;
  std::unique_ptr<Bundle> _extractKmers(uint64_t, ShortSuffix) const;
  std::unique_ptr<Bundle> _extractPartialKmers(uint64_t, ShortSuffix) const;

protected:
  // input data
  std::size_t s_;
  MetaColours *cmap_;
  push_summary *str_;

  // input buffers
  const std::vector<const SequenceContainer *> &seqs_;
  const TerminalRange &terminals_;

  // shared auxilliary data
  std::shared_ptr<std::vector<SuffixTable>> blocks_;
  std::shared_ptr<BufferMaker> buffers_;

  // output buffers
  std::shared_ptr<BundlePool> pool_;
};

// -----------------------------------------------------------------------------
// API
// -----------------------------------------------------------------------------

ColouredGraphFiles construct(const std::vector<const SequenceContainer *> &data,
                             MetaColours &&cmap, BuildOptions params = {});

} // namespace cdbg
