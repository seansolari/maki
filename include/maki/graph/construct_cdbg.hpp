
#pragma once

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

#include <filesystem>

namespace fs = std::filesystem;

struct BufferPaths
{
  fs::path base;
  fs::path edges;
  fs::path succ;
  fs::path colours;
};

using CDBGSinks = std::tuple<
    , // edges: std::vector<uint8_t> -> FILE
    , // succ: sdsl::bit_vector -> sdsl::int_vector_handle
    , // colours: std::vector<uint8_t> -> Archive
    >;

struct ColouredKmerBuffers
{
  // payload global position
  suffix s;
  // k-mer data buffers
  KmerBuffer kmers;
  KmerBuffer temp;
  TerminalRange terminals;
  overlap_vector;
};

void constructCDBG()
{
  BufferPaths outp { main };

  CDBGSinks Sinks{
    { outp.edges }, // consume edges
    { outp.succ },  // consume succ
    { outp.colours }
  };

  for (each suffix)
  {

  }
}
