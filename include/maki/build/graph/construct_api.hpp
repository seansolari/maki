#pragma once

// -----------------------------------------------------------------------------
// GFF data
// -----------------------------------------------------------------------------

ColouredGraphFiles constructGffColouredDbg(const GenomeManifest &, dbg::BuildOptions params);

// -----------------------------------------------------------------------------
// Fasta data
// -----------------------------------------------------------------------------

ColouredGraphFiles constructFnaColouredDbg(const GenomeManifest &, dbg::BuildOptions params);

// -----------------------------------------------------------------------------
// Fastq data
// -----------------------------------------------------------------------------

BaseGraphFiles constructFqDbg(const FastqDataset &, dbg::BuildOptions params);
WeightedGraphFiles constructFqWeightedDbg(const FastqDataset &, dbg::BuildOptions params);
