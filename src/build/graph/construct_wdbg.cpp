
#include "maki/build/graph/construct_wdbg.hpp"
#include "maki/core/graph/archive/counts.hpp"
#include "maki/core/utils/logging.hpp"
#include <oneapi/tbb/global_control.h>

namespace wdbg {

// -----------------------------------------------------------------------------
// Finalisation
// -----------------------------------------------------------------------------

WeightedGraphFiles finalise(TempBuffers inp, std::size_t k,
                            CountBuffer &&counts, const std::string &out) {

  LOG_INFO() << "Finalising weighted de Bruijn graph (k=" << k << ")";

  WeightedGraphFiles outp(out);
  dbg::detail::finaliseGraphBuffers(inp.files.edges, inp.files.succ, outp);

  WeightedGraph g;
  g.k = k;
  dbg::detail::finaliseGraphStructure(g, inp.str);

  LOG_INFO() << "Compressing k-mer counts buffer";
  const std::size_t orig_counts_size = detail::size_in_bytes(counts);
  g.occ = CompressedCountBuffer(std::move(counts));
  const std::size_t new_counts_size = detail::size_in_bytes(g.occ);
  LOG_INFO() << "Compression complete, uncmp=" << orig_counts_size
             << ", cmp=" << new_counts_size;

  dbg::detail::serialize(g, outp.meta);

  LOG_INFO() << "Graph construction complete";
  return outp;
}

} // namespace wdbg
