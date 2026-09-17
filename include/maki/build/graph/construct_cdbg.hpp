
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
#include "maki/build/kmers/construct_terminals.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/seq/concepts.hpp"
#include <cstddef>
#include <filesystem>
#include <memory>
#include <sdsl/int_vector.hpp>
#include <vector>

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

namespace dbg {

template <>
struct BufferMaker<KmerBuffer> : public Factory<Buffers<KmerBuffer>> {
  BufferMaker(std::size_t length_, std::size_t width_, std::size_t k_,
              std::size_t keff_)
      : Factory<Buffers<KmerBuffer>>(), length(length_), width(width_), k(k_),
        keff(keff_) {}

  std::size_t length, width, k, keff;

  inline std::unique_ptr<Buffers<KmerBuffer>> obtain() {
    return Factory<Buffers<KmerBuffer>>::obtain(length, width, k, keff);
  }
};

} // namespace dbg

namespace cdbg_detail {

using Buffers = dbg::Buffers<KmerBuffer>;
using BufferMaker = dbg::BufferMaker<KmerBuffer>;

// -----------------------------------------------------------------------------
// Output data
// -----------------------------------------------------------------------------

using EdgeSink = SdslIntVectorOnDiskSink<4>;
using SuccSink = SdslIntVectorOnDiskSink<1>;
using ColourSink = ArchiveWriter<>;

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

template <container_span T, class... Sinks>
struct SuffixwiseKmers : public dbg::Suffixwise<T, KmerBuffer, Sinks...> {
  SuffixwiseKmers(const T &seqs, const TerminalRange &terms,
                  std::shared_ptr<std::vector<SuffixTable>> &&suffixPlan,
                  std::shared_ptr<dbg::BufferMaker<KmerBuffer>> &&buffers,
                  std::size_t s, MetaColours *cmap, push_summary *str)
      : dbg::Suffixwise<T, KmerBuffer, Sinks...>(
            seqs, terms, std::move(suffixPlan), std::move(buffers), s, str),
        colourMap(cmap) {}

  static SuffixwiseKmers FromSequences(const T &seqs,
                                       const TerminalRange &terms,
                                       std::size_t k, std::size_t s,
                                       MetaColours *cmap, push_summary *str) {
    assert(k > s);
    auto suffixPlan = std::make_shared<std::vector<SuffixTable>>(
        createSuffixPlan(seqs, k, s, k - s));
    auto bufferFactory = std::make_shared<dbg::BufferMaker<KmerBuffer>>(
        suffixPlan->back().maxValue(), value_size(cmap->colourWidth()), k,
        k - s);
    return SuffixwiseKmers(seqs, terms, std::move(suffixPlan),
                           std::move(bufferFactory), s, cmap, str);
  }

  std::unique_ptr<ChunkBundleT<Sinks...>> operator()(uint64_t idx) const {
    auto sfx = ShortSuffix::fromIndex(idx, this->suffixSize);

    LOG_DEBUG() << "Processing suffix index " << idx << " (size=" << sfx.size()
                << ")";

    if (sfx.size() < this->suffixSize) {
      return extractPartialKmers(idx, sfx);
    } else {
      return extractKmers(idx, sfx);
    }
  }

  std::unique_ptr<ChunkBundleT<Sinks...>> extractKmers(uint64_t idx,
                                                       ShortSuffix sfx) const {
    auto bfr = this->kmerBuffers->obtain();

    bfr->collectKmers(this->sequences, *this->suffixCounts, sfx);
    bfr->setTerminals(this->terminals.endsWith(sfx));

    LOG_DEBUG() << "Interleaving " << bfr->kmers.size() << " k-mers and "
                << bfr->terminals.size() << " terminals";

    auto bnd = this->getBundle(idx);
    auto &[edges, succ, carch] = bnd->payloads;

    auto counts =
        interleave(bfr->kmers, bfr->b.begin(), bfr->terminals, bfr->t.begin(),
                   edges, succ, sfx.msb(), carch, *colourMap);

    this->pushRegionStructure(counts);

    this->kmerBuffers->release(bfr);
    return bnd;
  }

  std::unique_ptr<ChunkBundleT<Sinks...>>
  extractPartialKmers(uint64_t idx, ShortSuffix sfx) const {
    auto bfr = this->kmerBuffers->obtain();
    bfr->setTerminals(this->terminals.retrieve(sfx));

    const std::size_t n = bfr->terminals.size();

    LOG_DEBUG() << "Extracted " << n << " terminal (partial) k-mers for suffix "
                << sfx.toString();

    auto bnd = this->getBundle(idx);
    auto &[edges, succ, carch] = bnd->payloads;

    auto counts = pushRange(bfr->terminals, bfr->t.begin(), edges, succ,
                            sfx.msb(), carch, *colourMap);

    this->pushRegionStructure(counts);

    this->kmerBuffers->release(bfr);
    return bnd;
  }

  MetaColours *colourMap;
};

// -----------------------------------------------------------------------------
// Finalisation
// -----------------------------------------------------------------------------

ColouredGraphFiles finalise(TempBuffers inp, std::size_t k, MetaColours &&cols,
                            const std::string &out);

} // namespace cdbg_detail

// -----------------------------------------------------------------------------
// API
// -----------------------------------------------------------------------------

namespace cdbg {

using Multi = MultiSink<cdbg_detail::EdgeSink, cdbg_detail::SuccSink,
                        cdbg_detail::ColourSink>;

template <container_span T>
using SuffixwiseKmers = cdbg_detail::SuffixwiseKmers<
    T, cdbg_detail::EdgeSink, cdbg_detail::SuccSink, cdbg_detail::ColourSink>;

template <container_span T>
ColouredGraphFiles construct(const T &data, MetaColours &&cmap,
                             dbg::BuildOptions params = {}) {
  LOG_INFO() << "Starting CDBG construction";
  LOG_INFO() << "Input sequences: " << data.size();

  std::filesystem::create_directories(params.out);
  LOG_INFO() << "Output directory: " << params.out;

  cdbg_detail::TempBuffers outp{
      .files = {.edges = params.out / "temp-edges.sdsl",
                .succ = params.out / "temp-succ.sdsl",
                .colours = params.out / "temp-colours.maki"},
      .str = {}};

  LOG_INFO() << "Extracting terminal k-mers";
  auto terminals = extractTerminalsSparse(data, params.kmer_size);
  LOG_INFO() << "Total terminal entries = " << terminals.size();

  LOG_INFO() << "Initialising sinks";
  Multi sinks{cdbg_detail::EdgeSink(outp.files.edges),
              cdbg_detail::SuccSink(outp.files.succ),
              cdbg_detail::ColourSink(outp.files.colours)};

  LOG_INFO() << "Processing chunks";
  ProcessChunks(SuffixwiseKmers<T>::FromSequences(
                    data, terminals.asRange(), params.kmer_size,
                    params.suffix_size, &cmap, &outp.str),
                sinks, params.pool_size, params.reserve_per_chunk,
                params.chunks());

  LOG_INFO() << "Finalising sinks";
  sinks.finalize();

  return finalise(outp, params.kmer_size, std::move(cmap), params.out);
}

} // namespace cdbg
