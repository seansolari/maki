
#include "maki/build/graph/construct_common.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/core/utils/logging.hpp"

namespace dbg {

std::size_t BuildOptions::chunks() const {
  auto n = ShortSuffix::numSuffixes(suffix_size);
  LOG_INFO() << "Configured build to use " << n
             << " suffix chunks (suffix size = " << suffix_size << ")";
  return n;
}

} // namespace dbg
