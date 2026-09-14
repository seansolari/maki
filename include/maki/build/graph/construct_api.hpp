#pragma once

#include "maki/build/graph/construct_common.hpp"
#include "maki/build/io/all.hpp"

// -----------------------------------------------------------------------------
// GFF data
// -----------------------------------------------------------------------------

ColouredGraphFiles constructGffColouredDbg(const GenomeManifest &manifest, dbg::BuildOptions params) {
  // parse data
  Colours colours;
  auto genomes = parse(manifest.files, parseGFF, colours, (std::size_t)params.kmer_size);

  // construct graph
  // ...
}

// -----------------------------------------------------------------------------
// Fasta data
// -----------------------------------------------------------------------------

ColouredGraphFiles constructFnaColouredDbg(const fs::path &fastaFile, dbg::BuildOptions params) {
  Colours colours;
  auto genome = parseFilterFNA(fastaFile, colours, params.threads, params.kmer_size);

  // construct graph
  // ...
}

// -----------------------------------------------------------------------------
// Fastq data
// -----------------------------------------------------------------------------

BaseGraphFiles constructFqDbg(const reads::FastqDatasetChunkView &data, dbg::BuildOptions params) {

}

WeightedGraphFiles constructFqWeightedDbg(const reads::FastqDatasetChunkView &data, dbg::BuildOptions params) {

}
