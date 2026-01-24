
#pragma once

#include <cstddef>
#include <filesystem>
#include <vector>

#include "maki/build/kmers/buffers/kmers.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/build/kmers/construct_terminals.hpp"
#include "maki/build/graph/archive/archive_writer.hpp"
#include "maki/build/graph/archive/byte_writer.hpp"
#include "maki/build/graph/archive/sdsl_writer.hpp"
#include "maki/core/seq/seq_concepts.hpp"

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
  std::size_t threads = 16;
  std::size_t pool_size = 16;
  std::size_t reserve_per_chunk = 0;
};

struct ColouredKmerBuffers {
  KmerBuffer data;
  KmerBuffer temp;
};

/*

using CDBGSinks =
    std::tuple<ByteArraySink<>, // edges: std::vector<uint8_t> -> FILE
               SdslIntVectorOnDiskSink<1>, // succ: sdsl::bit_vector ->
                                           // sdsl::int_vector_handle
               ArchiveWriter<> // colours: std::vector<uint8_t> -> Archive
               >;

*/

struct BufferPaths {
  std::filesystem::path edges;
  std::filesystem::path succ;
  std::filesystem::path colours;
};

BufferPaths constructCDBG(const std::vector<const SequenceContainer *> &data,
                          BuildOptions params = {}) {
  BufferPaths outp{params.out / "edges.txt", params.out / "succ.sdsl",
                   params.out / "colours.maki"};

  // extract terminals and create suffix plan
  auto terminals = extractTerminalsSparse(data, params.kmer_size);
  auto splan = createSuffixPlan(data, params.kmer_size, params.suffix_size);
  std::size_t requiredBufferSize = splan.back().maxValue();

  return outp;
}
