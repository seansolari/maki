
#include "maki/build/graph/construct_cdbg.hpp"
#include "maki/build/graph/build_colours.hpp"
#include "maki/build/graph/construct_common.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/utils/logging.hpp"
#include <filesystem>
#include <oneapi/tbb/global_control.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

namespace cdbg {

// -----------------------------------------------------------------------------
// Finalisation
// -----------------------------------------------------------------------------

ColouredGraphFiles finalise(TempBuffers inp, std::size_t k, MetaColours &&cols,
                            const std::string &out) {

  LOG_INFO() << "Finalising coloured de Bruijn graph (k=" << k << ")";

  ColouredGraphFiles outp(out);
  dbg::detail::finaliseGraphBuffers(inp.files.edges, inp.files.succ, outp);

  LOG_INFO() << "Moving colour archive";
  std::filesystem::rename(inp.files.colours, outp.archive);

  ColouredGraph g;
  g.k = k;
  dbg::detail::finaliseGraphStructure(g, inp.str);
  g.cmap = toRegistry(std::move(cols));

  dbg::detail::serialize(g, outp.meta);

  LOG_INFO() << "Graph construction complete";
  return outp;
}

} // namespace cdbg
