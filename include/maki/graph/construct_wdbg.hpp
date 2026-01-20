
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

#include <filesystem>

namespace fs = std::filesystem;

struct Buffers
{
  sdsl::int_vector<4> edges;
  sdsl::bit_vector succ;
  std::vector<uint64_t> counts;
};

using WDBGSinks = std::tuple<
    , // edges: sdsl::int_vector<4> -> sdsl::int_vector<4>
    , // succ: sdsl::bit_vector -> sdsl::bit_vector
    , // counts: std::vector<uint64_t> -> std::vector<uint64_t>
    >;
