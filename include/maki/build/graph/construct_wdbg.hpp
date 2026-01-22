
#pragma once

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

#include <utility>

#include "maki/build/graph/archive/sdsl_writer.hpp"
#include "maki/build/graph/archive/vector_writer.hpp"

using WDBGSinks = std::tuple<
    SdslIntVectorInMemorySink<4>, // edges: sdsl::int_vector<4> -> sdsl::int_vector<4>
    SdslIntVectorInMemorySink<1>, // succ: sdsl::bit_vector -> sdsl::bit_vector
    VectorInMemorySink<uint64_t>  // counts: std::vector<uint64_t> -> std::vector<uint64_t>
    >;

struct Buffers
{
  sdsl::int_vector<4> edges;
  sdsl::bit_vector succ;
  std::vector<uint64_t> counts;
};
