#pragma once

#include "maki/build/graph/construct_common.hpp"
#include "maki/build/io/fastq.hpp"
#include "maki/core/graph/base.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/graph/wdbg.hpp"

ColouredGraphFiles constructGffColouredDbg(const GenomeManifest &manifest,
                                           dbg::BuildOptions params);

ColouredGraphFiles constructFnaColouredDbg(const fs::path &fastaFile,
                                           dbg::BuildOptions params);

DeBruijnGraphFiles constructFqDbg(const reads::FastqDatasetChunkView &data,
                                  dbg::BuildOptions params);

WeightedGraphFiles
constructFqWeightedDbg(const reads::FastqDatasetChunkView &data,
                       dbg::BuildOptions params);
