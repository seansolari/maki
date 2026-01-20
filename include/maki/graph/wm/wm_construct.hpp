
#pragma once
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

// A minimal sink interface capturing WM output
struct WMLevel {
  std::vector<uint64_t> words; // packed bitvector (64-bit blocks), length = ceil(n/64)
  uint64_t zeros = 0;          // Z[level]
};

struct WMBuildResult {
  uint32_t log_sigma = 0;              // for σ=16 => 4
  std::vector<WMLevel> levels;         // size = log_sigma
  uint64_t n = 0;                      // length of the text
};

// Build using pwm's *parallel semi-external prefix-counting* WM constructor.
// `sigma` in [1..256]. For σ=16, log_sigma = 4.
// The call mirrors how `src/benchmark` invokes this algorithm.
WMBuildResult build_wm_ppc_semi_external_from_file(
    const std::string& path,
    uint32_t sigma,
    int threads /* 0=auto */);
