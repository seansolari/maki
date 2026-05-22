
#include "maki/build/graph/construct_common.hpp"
#include "maki/build/graph/wm/wm_construct.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/core/utils/logging.hpp"

namespace dbg {

BuildOptions::BuildOptions(std::size_t k, std::size_t s,
                           const std::filesystem::path &out,
                           std::size_t threads)
    : kmer_size(k), suffix_size(s), out(out), threads(threads) {}

std::size_t BuildOptions::chunks() const {
  auto n = ShortSuffix::numSuffixes(suffix_size);
  LOG_INFO() << "Configured build to use " << n
             << " suffix chunks (suffix size = " << suffix_size << ")";
  return n;
}

void detail::initSuccSupport(DeBruijnGraphFiles &outp) {
  LOG_INFO() << "Initialising rank/select support for successor bitvector";

  sdsl::bit_vector succ;
  sdsl::load_from_file(succ, outp.l);

  LOG_INFO() << "Successor bitvector length = " << succ.size();

  {
    sdsl::rank_support_v5<1, 1> lRnk;
    sdsl::util::init_support(lRnk, &succ);
    sdsl::store_to_file(std::move(lRnk), outp.lR);
  }

  {
    sdsl::select_support_mcl<1, 1> lSel;
    sdsl::util::init_support(lSel, &succ);
    sdsl::store_to_file(std::move(lSel), outp.lS);
  }
}

void detail::finaliseGraphBuffers(const std::string &edges,
                                  const std::string &succ,
                                  DeBruijnGraphFiles &outp) {
  LOG_INFO() << "Constructing edge wavelet matrix";
  initW(edges, outp.W, outp.base);

  LOG_INFO() << "Removing temporary edge file: " << edges;
  std::filesystem::remove(edges);

  LOG_INFO() << "Moving successor array to final location";
  std::filesystem::rename(succ, outp.l);

  detail::initSuccSupport(outp);
}

void detail::finaliseGraphStructure(DeBruijnGraph &g, push_summary &str) {
  g.F[0] = 0u;
  g.F[1] = str.F[0];
  g.F[2] = str.F[0] + str.F[1];
  g.F[3] = str.F[0] + str.F[1] + str.F[2];
  g.F[4] = str.F[0] + str.F[1] + str.F[2] + str.F[3];

  g.C[0] = 0u;
  g.C[1] = str.C[0];
  g.C[2] = str.C[0] + str.C[1];
  g.C[3] = str.C[0] + str.C[1] + str.C[2];
  g.C[4] = str.C[0] + str.C[1] + str.C[2] + str.C[3];
}

} // namespace dbg
