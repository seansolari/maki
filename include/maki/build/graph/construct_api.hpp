#pragma once

#include "maki/build/graph/build_colours.hpp"
#include "maki/build/graph/construct_cdbg.hpp"
#include "maki/build/graph/construct_common.hpp"
#include "maki/build/graph/construct_wdbg.hpp"
#include "maki/build/io/filter.hpp"
#include "maki/build/io/gff.hpp"
#include "maki/build/io/fastq.hpp"

// -----------------------------------------------------------------------------
// GFF data
// -----------------------------------------------------------------------------

ColouredGraphFiles constructGffColouredDbg(const GenomeManifest &manifest,
                                           dbg::BuildOptions params) {
  // parse data
  Colours colours;
  auto genomes =
      parse(manifest.files, parseGFF, colours, (std::size_t)params.kmer_size);

  // construct graph
  std::span<const Dna4Genome> view = genomes;
  return cdbg::construct(view, MetaColours(std::move(colours.ids)), params);
}

// -----------------------------------------------------------------------------
// Fasta data
// -----------------------------------------------------------------------------

ColouredGraphFiles constructFnaColouredDbg(const fs::path &fastaFile,
                                           dbg::BuildOptions params) {
  // parse data
  Colours colours;
  auto genome =
      parseFilterFNA(fastaFile, colours, params.threads, params.kmer_size);

  // construct graph
  return cdbg::construct(genome.data(), MetaColours(std::move(colours.ids)),
                         params);
}

// -----------------------------------------------------------------------------
// Fastq data
// -----------------------------------------------------------------------------

// BaseGraphFiles constructFqDbg(const reads::FastqDatasetChunkView &data,
//                               dbg::BuildOptions params) {}

WeightedGraphFiles
constructFqWeightedDbg(const reads::FastqDatasetChunkView &data,
                       dbg::BuildOptions params) {
  return wdbg::construct(data, params);
}
