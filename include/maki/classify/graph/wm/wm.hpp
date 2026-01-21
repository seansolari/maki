
#pragma once
#include <sdsl/bit_vectors.hpp>
#include <string>
#include <vector>
#include <memory>

struct WMLevelSDSL
{
  sdsl::bit_vector bv;
  sdsl::rank_support_v5<> r1; // 6.25% overhead, faster than v
  sdsl::select_support_mcl<1> s1;
  sdsl::select_support_mcl<0> s0;
};

struct WaveletMatrixSDSL
{
  uint64_t n = 0;
  uint32_t log_sigma = 0;
  std::vector<uint64_t> Z; // zeros per level
  std::vector<WMLevelSDSL> L;

  uint8_t access(uint64_t i) const;
  uint64_t rank(uint8_t c, uint64_t i) const;
  uint64_t select(uint8_t c, uint64_t k) const;
};

// Load per-level SDSL files and build supports
WaveletMatrixSDSL load_wm_levels_from_files(const std::string &base_path,
                                            const std::vector<uint64_t> &Z,
                                            uint32_t log_sigma, uint64_t n);
