
#include "maki/build/graph/wm/wm_construct.hpp"
#include <cmath>
#include <cstdint>
#include <stdexcept>

// ---- pwm headers you will see in benchmark/register code ----
// These header names are representative; confirm exact names in
// `src/benchmark.cpp` and `src/register/parallel_prefix_counting.cpp`, then
// include those.
#include "pwm/external_memory/wx_ppc_ie.hpp"

// Prototype that matches the semi-external WM PPC call the benchmark uses.
// Find the exact symbol in `src/register/parallel_prefix_counting.cpp`.
//
// Typical shape (to confirm): a function that takes a file path (semi-external
// source), σ or logσ, #threads, and produces per-level bitvectors + Z[level].
//
// For clarity, we declare a function pointer we’ll resolve via a thin inline
// wrapper (included from the same headers the benchmark uses), so we keep names
// in one place.
namespace pwm_adapters {
// You will replace the body of `invoke_wm_ppc_semi_external` with the *exact*
// function call you locate in pwm's register/benchmark sources.
//
// It must:
//   - Read `path` as external text (semi-external streaming),
//   - Construct a WM with levels = ceil(log2(sigma)),
//   - For each level l, fill `out[l].words` with packed bits, and
//   `out[l].zeros` with Z[l].
//
// The exact iterator/sink type used by pwm can be adapted here.
static void invoke_wm_ppc_semi_external(const std::string &path, uint32_t sigma,
                                        int threads, std::vector<WMLevel> &out,
                                        uint64_t &n) {
  // -------------------------
  // BEGIN: glue to pwm call
  // -------------------------
  //
  // HOW TO FILL THIS IN:
  //   1) Open `external/pwm/src/benchmark.cpp` and find where `wm` + `ppc`
  //      + "semi-external" are selected.
  //   2) Note the *header* and the *function* used to perform construction.
  //      (You will likely see a call through a small registry wrapper.)
  //   3) Include that header above and call it here.
  //
  // PSEUDOCODE matching benchmark style (replace with actual symbol):
  //
  // auto log_sigma = static_cast<uint32_t>(std::ceil(std::log2(sigma)));
  // n = <get_input_length>(path);
  // out.resize(log_sigma);
  //
  // pwm::wm::ppc::semi_external::construct_from_file(
  //     path,
  //     sigma,
  //     threads,
  //     [&](uint/         out[level].words.assign(words, words + nwords);
  //         out[level].zeros = zeros;
  //     },
  //     n /* or return it */);
  //
  // -------------------------
  // END: glue to pwm call
  // -------------------------
  //
}

WMBuildResult build_wm_ppc_semi_external_from_file(const std::string &path,
                                                   uint32_t sigma,
                                                   int threads) {

  if (sigma == 0)
    throw std::invalid_argument("sigma must be >= 1");
  const uint32_t log_sigma = (sigma <= 1) ? 0u : 32u - __builtin_clz(sigma - 1);

  WMBuildResult res;
  res.log_sigma = log_sigma;
  res.levels.resize(log_sigma);
  pwm_adapters::invoke_wm_ppc_semi_external(path, sigma, threads, res.levels,
                                            res.n);
  return res;
}

} // namespace pwm_adapters
