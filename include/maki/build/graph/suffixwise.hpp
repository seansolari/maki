
#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>
#include <memory>
#include <mutex>
#include <oneapi/tbb/parallel_pipeline.h>

#include "maki/build/graph/archive/sink_manager.hpp"

struct BuildOptions
{
  // Algorithm parameters
  size_t suffix_size = 8;
  //
  size_t pool_size = 8;
  size_t reserve_bytes_per_chunk = 8ull << 20; // 8 MiB default reserve
  uint64_t expected_chunks = 0;                //
};

// Pipeline Stage 1: Get Next Suffix
//  void -> Suffix

struct GetNextSuffix
{
  size_t *i;      // current suffix
  const size_t S; // total number of suffixes

  std::unique_ptr<Suffix> operator()(oneapi::tbb::flow_control &fc_) const
  {
    if (*i < S)
    {
      return std::make_unique(Suffix::as_suffix(*i++));
    }
    else
    {
      fc_.stop();
      return {static_cast<Suffix *>(nullptr)};
    }
  }
};

// Pipeline Stage 2: Process Suffix Data
// Suffix -> Bundle

// Thread-safe stack of buffers
template <typename Buffer>
class ThreadSafeBufferStack
{
private:
  std::vector<std::unique_ptr<Buffer>> data;
  std::mutex mtx;

public:
  void push(std::unique_ptr<Buffer> buffer)
  {
    std::lock_guard<std::mutex> lock(mtx);
    data.push_back(std::move(buffer));
  }

  std::unique_ptr<Buffer> pop()
  {
    std::lock_guard<std::mutex> lock(mtx);
    if (data.empty())
    {
      return {}; // Or throw exception
    }
    std::unique_ptr<Buffer> buffer = std::move(data.back());
    data.pop_back();
    return buffer;
  }
};

template <class... Sinks>
using SharedPipelinePointer = std::shared_ptr<PipelineFacade<Sinks...>>;

template <class Buffer, class... Sinks>
struct ProcessSuffixData
{
  std::shared_ptr<ThreadSafeBufferStack<Buffer>> _buffers;
  SharedPipelinePointer<Sinks...> _mgr;

  std::unique_ptr<ChunkBundleT<Sinks...>> operator()(std::unique_ptr<Suffix> sx_) const
  {
    // get free buffer for k-mer processing
    auto buf = _buffers->pop();
    buf->extractSuffix(sx_);
    // flush to output
    auto bdl = _mgr->pool().acquire();
    bdl->id = sx_.as_index();
    buf->flush(bdl->payloads);
    // clear k-mer buffer and put back on stack
    buf->clear();
    _buffers->push(std::move(buf));
    // pass payload for archiving
    return bdl;
  }
};

// Pipeline Stage 3: Flush Bundle
// Bundle -> void

template <class... Sinks>
struct FlushBundle
{
  SharedPipelinePointer<Sinks...> _mgr;

  void operator()(std::unique_ptr<ChunkBundleT<Sinks...>> bdl_) const
  {
    _mgr->queue().submit(std::move(bdl_));
  }
};

/**
 * Suffixwise Fill
 * ---------------
 *
 * Construct buffers in a suffix-wise fashion.
 *
 *  1) From the input data, extract k-mers and terminals with a given suffix
 *      into temporary buffers (multi-threaded).
 *  2) Sort k-mers and establish overlap (multi-threaded).
 *  3) Flush k-mers into graph buffers (single thread).
 *      a) Edges written to flat disk array.
 *      b) Node boundaries written to flat disk array.
 *      c) Colours encoded and written to archive.
 *      d) k-mer counts written to in-memory array.
 *
 */
template <typename Cursor_, typename Unwind_, typename Kmers_>
void suffixwiseFill(Cursor_ &crs_,
                    std::vector<range_t<Unwind_>> const &data,
                    Unwind_ Apply,
                    std::vector<suffix::SuffixTable> const &blocks,
                    Kmers_ &kmers,
                    Kmers_ &temp,
                    terminals::TerminalRange &tmls,
                    suffix::SmallRollingNuclSeq suffix_,
                    uint8_t s_,
                    uint32_t threads)
{
  if (suffix_.size() == s_)
  {
    LOG(INFO) << "inserting suffix " << suffix_.toString();

    extractAndSortSuffix(data, Apply, blocks, suffix_, kmers, temp, threads);
    terminals::TerminalRange tx = tmls.endsWith(suffix_);

    flush::flush(crs_, kmers, tx,
                 BaseGraph::LOCK_SAMPLE_RATE,
                 static_cast<uint8_t>(mers::dna4ToDna5(suffix_.msb())));
  }
  else
  {
    insertPartialSuffix(crs_, tmls, suffix_);

    /* A */ suffixwiseFill(crs_, data, Apply, blocks, kmers, temp, tmls,
                           suffix_ + uint64_t(0b00), s_, threads);
    /* C */ suffixwiseFill(crs_, data, Apply, blocks, kmers, temp, tmls,
                           suffix_ + uint64_t(0b01), s_, threads);
    /* G */ suffixwiseFill(crs_, data, Apply, blocks, kmers, temp, tmls,
                           suffix_ + uint64_t(0b10), s_, threads);
    /* T */ suffixwiseFill(crs_, data, Apply, blocks, kmers, temp, tmls,
                           suffix_ + uint64_t(0b11), s_, threads);
  }
}

template <typename Unwind_, typename Kmers_>
void extractAndSortSuffix(std::vector<range_t<Unwind_>> const &data,
                          Unwind_ Apply,
                          std::vector<suffix::SuffixTable> const &counts,
                          suffix::SmallRollingNuclSeq suffix_,
                          Kmers_ &kmers,
                          Kmers_ &temp,
                          uint32_t threads)
{
  // block start for each genome

  auto blocks = accumulatePartials(counts, [&](suffix::SuffixTable const &st) -> size_t
                                   { return st[suffix_]; });
  size_t nRecs = blocks.back();

  LOG(INFO) << "inserting " << nRecs << " k-mers with this suffix";

  // extract and sort k-mers with this suffix

  kmers.resize(nRecs);
  temp.resize(nRecs);
  kmers.insertKmers(data, Apply, blocks, suffix_);
  if (nRecs > 1)
    kmers.sort(&temp, threads);
}

template <typename Cursor_>
void insertPartialSuffix(Cursor_ &crs_,
                         terminals::TerminalRange &tmls,
                         suffix::SmallRollingNuclSeq suffix_)
{
  auto vals = tmls.retrieve(suffix_);
  if (!vals.empty())
  {
    flush::flush(crs_, vals, vals.size());
  }
}