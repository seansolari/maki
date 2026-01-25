
#include "maki/build/graph/wm/wm_construct.hpp"
#include <cmath>
#include <cstdint>
#include <omp.h>
#include <pwm/construction/wavelet_structure.hpp>
#include <pwm/external_memory/wx_ppc_ie.hpp>
#include <pwm/util/stats.hpp>

// Prototype that matches the semi-external WM PPC call the benchmark uses.
namespace pwm_adapters {
using wm_8_char = uint8_t;
using wm_ppc_8_ie = wx_ppc_ie<wm_8_char, false>;
using stats_type =
    statistics<wm_ppc_8_ie::external_in | wm_ppc_8_ie::external_out>;

/**
 * Given `n` the number of characters in the input file, and `log_sigma` the
 * bit-width of the alphabet, construct the wavelet matrix using PWM's
 * semi-external parallel prefix counting algorithm.
 */
static void invoke_wm_ppc_semi_external(const std::string &path, uint64_t n,
                                        uint64_t log_sigma, int threads,
                                        std::vector<WMLevel> &out) {
  omp_set_num_threads(threads);
  // call construction method - note by design the alphabet is minimal, so no
  // need for alphabet reduction stage
  stats_type stats;
  stats.start();
  wavelet_structure result = wm_ppc_8_ie::compute(path, n, log_sigma, stats);
  stats.finish();
  std::cerr << "Wavelet matrix constructed in " << stats.get_total_time()
            << "ms (" << path << ");\n";
  // wrap output
  
}

WMBuildResult construct_wm__pwm(const std::string &path, uint64_t n,
                                uint64_t log_sigma, int threads) {

  WMBuildResult res;
  res.log_sigma = log_sigma;
  res.levels.resize(log_sigma);

  pwm_adapters::invoke_wm_ppc_semi_external(path, n, log_sigma, threads, ...);
  return res;
}
} // namespace pwm_adapters
