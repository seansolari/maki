
#include "maki/graph/wm/wm.hpp"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>

static std::string level_filename(const std::string &base, uint32_t level)
{
  return base + ".wmlevel" + std::to_string(level) + ".sdsl";
}

WaveletMatrixSDSL load_wm_levels_from_files(const std::string &base_path,
                                            const std::vector<uint64_t> &Z,
                                            uint32_t log_sigma, uint64_t n)
{
  WaveletMatrixSDSL wm;
  wm.n = n;
  wm.log_sigma = log_sigma;
  wm.Z = Z;
  wm.L.resize(log_sigma);

  for (uint32_t l = 0; l < log_sigma; ++l)
  {
    auto f = level_filename(base_path, l);
    sdsl::load_from_file(wm.L[l].bv, f);                 // linear read
    sdsl::util::init_support(wm.L[l].r1, &(wm.L[l].bv)); // rank
    sdsl::util::init_support(wm.L[l].s1, &(wm.L[l].bv)); // select1
    sdsl::util::init_support(wm.L[l].s0, &(wm.L[l].bv)); // select0
  }
  return wm;
}

// WM navigation using per-level Z and SDSL rank/select (σ<=16 -> 4 levels)
static inline uint32_t bit_at(uint8_t c, uint32_t level, uint32_t log_sigma)
{
  return (c >> ((log_sigma - 1) - level)) & 1u;
}

uint8_t WaveletMatrixSDSL::access(uint64_t i) const
{
  uint8_t c = 0;
  uint64_t pos = i;
  for (uint32_t l = 0; l < log_sigma; ++l)
  {
    auto &Lvl = L[l];
    uint8_t b = Lvl.bv[pos];
    c = (uint8_t)((c << 1) | b);
    if (b == 0)
    {
      uint64_t ones = Lvl.r1.rank(pos);
      pos = pos - ones; // zeros before pos
    }
    else
    {
      uint64_t ones = Lvl.r1.rank(pos);
      pos = Z[l] + ones;
    }
  }
  return c;
}

uint64_t WaveletMatrixSDSL::rank(uint8_t c, uint64_t i) const
{
  if (i == 0)
    return 0;
  uint64_t Lp = 0, Rp = i;
  for (uint32_t l = 0; l < log_sigma; ++l)
  {
    auto &Lvl = L[l];
    uint8_t b = bit_at(c, l, log_sigma);
    uint64_t L1 = Lvl.r1.rank(Lp);
    uint64_t R1 = Lvl.r1.rank(Rp);
    if (b == 0)
    {
      Lp = Lp - L1;
      Rp = Rp - R1;
    }
    else
    {
      Lp = Z[l] + L1;
      Rp = Z[l] + R1;
    }
  }
  return Rp - Lp;
}

uint64_t WaveletMatrixSDSL::select(uint8_t c, uint64_t k) const
{
  if (k == 0)
    return uint64_t(-1);

  // Descend: find [L,R) range for c on the last level
  uint64_t Lp = 0, Rp = n;
  for (uint32_t l = 0; l < log_sigma; ++l)
  {
    auto &Lvl = L[l];
    uint8_t b = bit_at(c, l, log_sigma);
    uint64_t L1 = Lvl.r1.rank(Lp);
    uint64_t R1 = Lvl.r1.rank(Rp);
    if (b == 0)
    {
      Lp = Lp - L1;
      Rp = Rp - R1;
    }
    else
    {
      Lp = Z[l] + L1;
      Rp = Z[l] + R1;
    }
  }
  if (k > (Rp - Lp))
    return uint64_t(-1);

  uint64_t idx = Lp + (k - 1); // bottom index among b-bucket
  // Ascend: invert per level using select0 / select1
  for (int l = (int)log_sigma - 1; l >= 0; --l)
  {
    auto &Lvl = L[(size_t)l];
    uint8_t b = bit_at(c, (uint32_t)l, log_sigma);
    idx = (b == 0) ? Lvl.s0.select(idx + 1)
                   : Lvl.s1.select(idx + 1);
  }
  return idx;
}
