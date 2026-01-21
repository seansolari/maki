
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

#include "maki/io/fasta.hpp"
#include "maki/graph/archive/byte_writer.hpp"
#include "maki/graph/archive/sdsl_writer.hpp"
#include "maki/graph/archive/archive_writer.hpp"

namespace fs = std::filesystem;

using CDBGSinks = std::tuple<
    ByteArraySink<>,            // edges: std::vector<uint8_t> -> FILE
    SdslIntVectorOnDiskSink<1>, // succ: sdsl::bit_vector -> sdsl::int_vector_handle
    ArchiveWriter<>             // colours: std::vector<uint8_t> -> Archive
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

struct BufferPaths
{
  fs::path edges;
  fs::path succ;
  fs::path colours;
};

struct BuildParams
{
  size_t k;
  size_t s;
  fs::path out;
};

BufferPaths constructCDBG(const std::vector<SequenceContainer*> &data, BuildParams params)
{
  BufferPaths outp{
    params.out / "edges.txt",
    params.out / "succ.sdsl",
    params.out / "colours.maki"
  };

  // extract terminals
  

  // create suffix plan

  // prepare output buffers

  CDBGSinks Sinks{
      {outp.edges},    // consume edges
      {outp.succ},     // consume succ
      {outp.colours}}; // consume colours

  // suffix-wise processing

  for (each suffix)
  {
  }

  return outp;
}
