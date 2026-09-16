
#pragma once

#include "construct_common.hpp"
#include "maki/build/graph/archive/counts_writer.hpp"
#include "maki/build/graph/archive/sdsl_writer.hpp"
#include "maki/build/graph/archive/sink_manager.hpp"
#include "maki/build/graph/interleave_buffers.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/core/graph/wdbg.hpp"
#include "maki/core/seq/concepts.hpp"
#include "maki/core/utils/algo.hpp"

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

namespace dbg_detail {

using Buffers = dbg::Buffers<TerminalBuffer>;
using BufferMaker = dbg::BufferMaker<TerminalBuffer>;

// -----------------------------------------------------------------------------
// Output data
// -----------------------------------------------------------------------------

using EdgeSink = SdslIntVectorOnDiskSink<4>;
using SuccSink = SdslIntVectorOnDiskSink<1>;
using CountSink = CountBufferSink;

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

template <container_span T, class... Sinks>
struct SuffixwiseTerminals
    : public dbg::Suffixwise<T, TerminalBuffer, Sinks...> {

  using dbg::Suffixwise<T, TerminalBuffer, Sinks...>::Suffixwise;

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

  std::unique_ptr<ChunkBundleT<Sinks...>> operator()(uint64_t idx) const {
    auto sfx = ShortSuffix::fromIndex(idx, this->suffixSize);

    LOG_DEBUG() << "Processing suffix index " << idx << " (size=" << sfx.size()
                << ")";

    if (sfx.size() < this->suffixSize) {
      return extractPartialSuffix(idx, sfx);
    } else {
      return extractSuffix(idx, sfx);
    }
  }

  std::unique_ptr<ChunkBundleT<Sinks...>> extractSuffix(uint64_t idx,
                                                        ShortSuffix sfx) const {
    auto bfr = this->kmerBuffers->obtain();

    bfr->collectKmers(this->sequences, *this->suffixCounts, sfx);

    const std::size_t n = bfr->kmers.size();

    LOG_DEBUG() << "Extracted " << n << " k-mers for suffix " << sfx.toString();

    auto bnd = this->getBundle(idx);
    auto &edges = std::get<0>(bnd->paylods);
    auto &succ = std::get<1>(bnd->payloads);
    auto counts = call_tail(
        [&](auto &&...args) {
          return pushRange(bfr->kmers, bfr->b.begin(), edges, succ, sfx.msb(),
                           std::forward<decltype(args)>(args)...);
        },
        bnd->payloads);

    this->pushRegionStructure(counts);

    this->kmerBuffers->release(bfr);
    return bnd;
  }

  std::unique_ptr<ChunkBundleT<Sinks...>>
  extractPartialSuffix(uint64_t idx, ShortSuffix sfx) const {
    auto bfr = this->kmerBuffers->obtain();
    bfr->setTerminals(this->terminals.retrieve(sfx));

    const std::size_t n = bfr->terminals.size();

    LOG_DEBUG() << "Extracted " << n << " terminal (partial) k-mers for suffix "
                << sfx.toString();

    auto bnd = this->getBundle(idx);
    auto &edges = std::get<0>(bnd->paylods);
    auto &succ = std::get<1>(bnd->payloads);
    auto counts = call_tail(
        [&](auto &&...args) {
          return pushRange(bfr->terminals, bfr->t.begin(), edges, succ,
                           sfx.msb(), std::forward<decltype(args)>(args)...);
        },
        bnd->payloads);

    this->pushRegionStructure(counts);

    this->kmerBuffers->release(bfr);
    return bnd;
  }
};

// -----------------------------------------------------------------------------
// Finalisation
// -----------------------------------------------------------------------------

DeBruijnGraphFiles finalise(TempBuffers inp, std::size_t k,
                            const std::string &out);

WeightedGraphFiles finalise(TempBuffers inp, std::size_t k,
                            CountBuffer &&counts, const std::string &out);

} // namespace dbg_detail

// -----------------------------------------------------------------------------
// API
// -----------------------------------------------------------------------------

namespace dbg {

using Multi = MultiSink<dbg_detail::EdgeSink, dbg_detail::SuccSink>;

template <container_span T>
using SuffixwiseTerminals =
    dbg_detail::SuffixwiseTerminals<T, dbg_detail::EdgeSink,
                                    dbg_detail::SuccSink>;

template <container_span T>
DeBruijnGraphFiles construct(const T &data, dbg::BuildOptions params = {}) {
  LOG_INFO() << "Starting DBG construction, node size=" << params.kmer_size;
  LOG_INFO() << "Input chunks: " << data.size();

  std::filesystem::create_directories(params.out);
  LOG_INFO() << "Output directory: " << params.out;

  dbg_detail::TempBuffers outp{
      .files = {.edges = params.out / "temp-edges.sdsl",
                .succ = params.out / "temp-succ.sdsl"},
      .str = {}};

  LOG_INFO() << "Extracting terminal k-mers";
  auto terminals = extractTerminalsDense(data, params.suffix_size);
  LOG_INFO() << "Total terminal entries = " << terminals.size();

  Multi sinks{dbg_detail::EdgeSink(outp.files.edges),
              dbg_detail::SuccSink(outp.files.succ)};

  ProcessChunks(SuffixwiseTerminals<T>::FromSequences(
                    data, terminals.asRange(), params.kmer_size,
                    params.suffix_size, &outp.str),
                sinks, params.pool_size, params.reserve_per_chunk,
                params.chunks());

  LOG_INFO() << "Finalising sinks";
  sinks.finalize();

  return dbg_detail::finalise(outp, params.kmer_size, params.out);
}

} // namespace dbg

namespace wdbg {

using Multi = MultiSink<dbg_detail::EdgeSink, dbg_detail::SuccSink,
                        dbg_detail::CountSink>;

template <container_span T>
using SuffixwiseTerminals = dbg_detail::SuffixwiseTerminals<
    T, dbg_detail::EdgeSink, dbg_detail::SuccSink, dbg_detail::CountSink>;

template <container_span T>
WeightedGraphFiles construct(const T &data, dbg::BuildOptions params = {}) {
  LOG_INFO() << "Starting WDBG construction, node size=" << params.kmer_size;
  LOG_INFO() << "Input chunks: " << data.size();

  std::filesystem::create_directories(params.out);
  LOG_INFO() << "Output directory: " << params.out;

  dbg_detail::TempBuffers outp{
      .files = {.edges = params.out / "temp-edges.sdsl",
                .succ = params.out / "temp-succ.sdsl"},
      .str = {}};

  // Graph edge counts
  CountBuffer rawCounts;

  LOG_INFO() << "Extracting terminal k-mers";
  auto terminals = extractTerminalsDense(data, params.suffix_size);
  LOG_INFO() << "Total terminal entries = " << terminals.size();

  Multi sinks{dbg_detail::EdgeSink(outp.files.edges),
              dbg_detail::SuccSink(outp.files.succ),
              dbg_detail::CountSink(&rawCounts)};

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
