
#include "maki/build/graph/construct_api.hpp"

// -----------------------------------------------------------------------------
// GFF data
// -----------------------------------------------------------------------------

ColouredGraphFiles constructGffColouredDbg(const GenomeManifest &im, dbg::BuildOptions params) {

  oneapi::tbb::global_control global_limit(
      oneapi::tbb::global_control::max_allowed_parallelism, params.threads);

  switch (im.type) {
  case InputFileType::Gff3FileType: {
    Colours colours;
    auto genomes =
        parse(im.files, parseGFF, colours, (std::size_t)params.kmer_size);
    auto view = toView(genomes);
    return construct(view, std::move(colours.ids), params);
  }
  case InputFileType::FastaFileType: {
    auto genomes = parse(im.files, parseFilterFNA, (std::size_t)params.threads,
                         (std::size_t)params.kmer_size);
    auto view = toView(genomes);
    return construct(view, params);
  }
  default:
    throw std::runtime_error(
        "Coloured graph construction only supported for FASTA or GFF3 files.");
  }
}

// -----------------------------------------------------------------------------
// Fasta data
// -----------------------------------------------------------------------------

ColouredGraphFiles constructFnaColouredDbg(const GenomeManifest &, dbg::BuildOptions params) {
  
}

// -----------------------------------------------------------------------------
// Fastq data
// -----------------------------------------------------------------------------

/*
OLD:

WeightedGraphFiles construct(const DataFilePair &fp, dbg::BuildOptions params) {
  oneapi::tbb::global_control global_limit(
      oneapi::tbb::global_control::max_allowed_parallelism, params.threads);

  auto chunks =
      chunkReads(detail::parsePairedFastq(fp, params.kmer_size,
                                          params.kmer_size, params.threads),
                 10 * params.threads);
  auto view = toView(chunks);
  return construct(view, params);
}

*/

BaseGraphFiles constructFqDbg(const FastqDataset &, dbg::BuildOptions params) {
  
}

WeightedGraphFiles constructFqWeightedDbg(const FastqDataset &, dbg::BuildOptions params) {

}

