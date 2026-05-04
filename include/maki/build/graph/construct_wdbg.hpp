
#pragma once

#include "construct_common.hpp"
#include "maki/build/graph/archive/sdsl_writer.hpp"
#include "maki/build/graph/archive/vector_writer.hpp"
#include "maki/build/graph/interleave_buffers.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/core/graph/wdbg.hpp"
#include "maki/core/seq/concepts.hpp"

namespace wdbg {

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

// -----------------------------------------------------------------------------
// Intermediate data
// -----------------------------------------------------------------------------

struct Buffers {
  Buffers(std::size_t length, std::size_t width, std::size_t k,
          std::size_t keff);

  /**
   * Collect k-mers with given suffix `s_` from input sequences.
   */
  void collectKmers(const std::vector<const SequenceContainer *> &seqs_,
                    const std::vector<SuffixTable> &blocks_, ShortSuffix s_);

  /**
   * Load terminals data.
   */
  void setTerminals(TerminalRange &&t_);

  TerminalBuffer kmers, temp;
  TerminalRange terminals;
  sdsl::int_vector<2> b, t;
};

// -----------------------------------------------------------------------------
// Output data
// -----------------------------------------------------------------------------

using EdgeSink = SdslIntVectorOnDiskSink<4>;
using SuccSink = SdslIntVectorOnDiskSink<1>;


using WDBGSinks =
    std::tuple<SdslIntVectorInMemorySink<4>, // edges: sdsl::int_vector<4> ->
                                             // sdsl::int_vector<4>
               SdslIntVectorInMemorySink<1>, // succ: sdsl::bit_vector ->
                                             // sdsl::bit_vector
               VectorInMemorySink<uint64_t>  // counts: std::vector<uint64_t> ->
                                             // std::vector<uint64_t>
               >;

// -----------------------------------------------------------------------------
// Pipeline
// -----------------------------------------------------------------------------

struct SuffixwiseTerminals {

  void setPool(std::shared_ptr<BundlePool> &p);
  std::unique_ptr<Bundle> operator()(uint64_t) const;

protected:
  // input data
  std::size_t s_;
  push_summary *str_;

  // input buffers
  const std::vector<const SequenceContainer *> &seqs_;
  const TerminalRange &terminals_;
};

// -----------------------------------------------------------------------------
// API
// -----------------------------------------------------------------------------

WeightedGraphFiles construct(const std::vector<const SequenceContainer *> &data,
                             dbg::BuildOptions params = {});

} // namespace wdbg
