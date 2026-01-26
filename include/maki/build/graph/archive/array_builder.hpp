
#pragma once
#include <algorithm>
#include <cstdint>
#include <vector>

class ArrayBuilder {
public:
  void reserve(std::size_t size_) {
    values_.reserve(size_);
  }

  uint64_t get(std::size_t i) {
    return values_[i];
  }

  void push(uint64_t v) {
    values_.push_back(v);
    max_ = std::max(max_, v);
  }

  void clear() {
    values_.clear();
    max_ = 0;
  }

  uint64_t size() const { return values_.size(); }

  uint8_t required_bit_width() const {
    if (max_ == 0)
      return 1;
    return 64 - __builtin_clzll(max_);
  }

  // Produce RAW_PACKED byte vector
  void finalize_packed(std::vector<uint8_t> &out) const {
    const uint64_t n = values_.size();
    const uint8_t w = required_bit_width();
    const __uint128_t total_bits = (__uint128_t)n * w;
    const size_t bytes_needed = (size_t)((total_bits + 7) / 8);
    out.assign(bytes_needed, 0);

    __uint128_t bitpos = 0;
    for (uint64_t i = 0; i < n; ++i) {
      uint64_t v = values_[i];
      if (w < 64)
        v &= ((uint64_t(1) << w) - 1u);
      const size_t byte_idx = (size_t)(bitpos >> 3);
      const unsigned bit_in_byte = (unsigned)(bitpos & 7);

      __uint128_t lane = (__uint128_t)v << bit_in_byte;
      for (unsigned b = 0; b < 16; ++b) {
        uint64_t shifted = (uint64_t)((lane >> (8 * b)) & 0xFF);
        if (!shifted)
          continue;
        size_t p = byte_idx + b;
        if (p < out.size())
          out[p] |= (uint8_t)shifted;
      }
      bitpos += w;
    }
  }

private:
  std::vector<uint64_t> values_;
  uint64_t max_ = 0;
};
