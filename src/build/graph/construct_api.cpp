
#include "maki/build/graph/construct_api.hpp"
#include "maki/build/graph/construct_cdbg.hpp"
#include "maki/build/graph/construct_wdbg.hpp"
#include "maki/build/io/filter.hpp"
#include "maki/build/io/gff.hpp"
#include <oneapi/tbb/global_control.h>

// -----------------------------------------------------------------------------
// GFF data
// -----------------------------------------------------------------------------

ColouredGraphFiles constructGffColouredDbg(const GenomeManifest &manifest,
                                           dbg::BuildOptions params) {

  oneapi::tbb::global_control global_limit(
      oneapi::tbb::global_control::max_allowed_parallelism, params.threads);

  // parse data
  Colours colours;
  auto genomes =
      parse(manifest.files, parseGFF, colours, (std::size_t)params.kmer_size);

  // construct graph
  return cdbg::construct(genomes, MetaColours(std::move(colours.ids)), params);
}

// -----------------------------------------------------------------------------
// Fasta data
// -----------------------------------------------------------------------------

ColouredGraphFiles constructFnaColouredDbg(const fs::path &fastaFile,
                                           dbg::BuildOptions params) {
  // parse data
  Colours colours;
  auto genome = parseGenome(fastaFile, colours, params.kmer_size);

  // chunk genome
  std::size_t chunksize =
      ChunkedDna4Genome::calculateChunksize(genome, params.threads);
  auto xgenome = chunkFNA(std::move(genome), chunksize, params.kmer_size);

  // construct graph
  return cdbg::construct(xgenome, MetaColours(std::move(colours.ids)), params);
}

// -----------------------------------------------------------------------------
// Fastq data
// -----------------------------------------------------------------------------

DeBruijnGraphFiles constructFqDbg(const reads::FastqDatasetChunkView &data,
                                  dbg::BuildOptions params) {
  return dbg::construct(data, params);
}

WeightedGraphFiles
constructFqWeightedDbg(const reads::FastqDatasetChunkView &data,
                       dbg::BuildOptions params) {
  return wdbg::construct(data, params);
}
