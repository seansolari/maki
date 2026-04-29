
#include "maki/build/graph/construct_wdbg.hpp"
#include "maki/build/kmers/construct_terminals.hpp"
#include "maki/core/utils/logging.hpp"

namespace wdbg {

// -----------------------------------------------------------------------------
// Intermediate data
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
// Pipeline
// -----------------------------------------------------------------------------

void SuffixwiseTerminals::setPool(std::shared_ptr<BundlePool> &p) { pool_ = p; }

std::unique_ptr<Bundle> SuffixwiseTerminals::operator()(uint64_t idx) const {
  auto sfx = ShortSuffix::fromIndex(idx, s_);

  LOG_DEBUG() << "Processing suffix index " << idx << " (size=" << sfx.size()
              << ")";

  if (sfx.size() < s_) {
    return _extractPartialKmers(idx, sfx);
  } else {
    return _extractKmers(idx, sfx);
  }
}

// -----------------------------------------------------------------------------
// API
// -----------------------------------------------------------------------------

WeightedGraphFiles construct(const std::vector<const SequenceContainer *> &data,
                             dbg::BuildOptions params) {

  LOG_INFO() << "Starting WDBG construction";
  LOG_INFO() << "Input chunks: " << data.size();
  LOG_INFO() << "Output directory: " << params.out;

  LOG_INFO() << "Extracting terminal k-mers";
  auto terminals = extractTerminalsDense(data, params.suffix_size);
  LOG_INFO() << "Total terminal entries = " << terminals.size();

  // suffix-wise k-mer processing
  // ...

  
}

}
