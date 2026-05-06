
#include "maki/build/graph/construct_wdbg.hpp"
#include "maki/build/kmers/construct_terminals.hpp"
#include "maki/core/utils/logging.hpp"

namespace wdbg {

// -----------------------------------------------------------------------------
// Pipeline
// -----------------------------------------------------------------------------

std::unique_ptr<Bundle> SuffixwiseTerminals::operator()(uint64_t idx) const {
  auto sfx = ShortSuffix::fromIndex(idx, suffixSize);

  LOG_DEBUG() << "Processing suffix index " << idx << " (size=" << sfx.size()
              << ")";

  if (sfx.size() < suffixSize) {
    return extractPartialSuffix(idx, sfx);
  } else {
    return extractSuffix(idx, sfx);
  }
}

std::unique_ptr<Bundle>
SuffixwiseTerminals::extractPartialSuffix(uint64_t idx, ShortSuffix sfx) const {

}

std::unique_ptr<Bundle>
SuffixwiseTerminals::extractSuffix(uint64_t idx, ShortSuffix sfx) const {

}

// -----------------------------------------------------------------------------
// API
// -----------------------------------------------------------------------------

WeightedGraphFiles construct(const std::vector<const SequenceContainer *> &data,
                             dbg::BuildOptions params) {

  LOG_INFO() << "Starting WDBG construction";
  LOG_INFO() << "Input chunks: " << data.size();
  LOG_INFO() << "Output directory: " << params.out;

  TempBuffers outp{.files = {.edges = params.out / "temp-edges.sdsl",
                             .succ = params.out / "temp-succ.sdsl"},
                   .str = {}};
  
  // Graph edge counts
  CountBuffer rawCounts;

  LOG_INFO() << "Extracting terminal k-mers";
  auto terminals = extractTerminalsDense(data, params.suffix_size);
  LOG_INFO() << "Total terminal entries = " << terminals.size();

  
  Multi sinks{EdgeSink(outp.files.edges), SuccSink(outp.files.succ),
              CountSink(&rawCounts)};

  
}

} // namespace wdbg
