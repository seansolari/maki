
#pragma once
#include <algorithm>
#include <cstdint>
#include <span>
#include <vector>

/**
 * @brief Create bit-packed vectors.
 * @note Holds raw values until `clear()` is called.
 */
class ArrayBuilder {
public:
  /**
   * @brief Reserve slots.
   * @param size_ Number of slots to reserve.
   */
  void reserve(std::size_t size_) {
    values_.reserve(size_);
  }

  /**
   * @brief `true` if there are no values to pack.
   * @return true 
   * @return false 
   */
  bool empty() const noexcept {
    return values_.empty();
  }

  /**
   * @brief Get value at slot `i`.
   * @param i
   * @return uint64_t 
   */
  uint64_t get(std::size_t i) {
    return values_[i];
  }

  /**
   * @brief Push value.
   * @param v Value to append.
   */
  void push(uint64_t v) {
    values_.push_back(v);
    max_ = std::max(max_, v);
  }

  /**
   * @brief Reset container (calls `clear()` on values).
   */
  void clear() {
    values_.clear();
    max_ = 0;
  }

  /**
   * @brief Number of values in container.
   * @return uint64_t 
   */
  uint64_t size() const {
    return values_.size();
  }

  /**
   * @brief Calculate bit width based on current values in container.
   * @note Calculated based on bits required to represent largest current value.
   * @return uint8_t 
   */
  uint8_t required_bit_width() const {
    if (max_ == 0)
      return 1;
    return 64 - __builtin_clzll(max_);
  }

  /**
   * @brief Write bit-packed values to external container.
   * @note Resets size and contents of the external container.
   * @param out 
   */
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

  /**
   * @brief Get view on current values.
   * @return std::span<const uint64_t> 
   */
  std::span<const uint64_t> view() const {
    return std::span<const uint64_t>(values_.begin(), values_.end());
  }

private:
  std::vector<uint64_t> values_;
  uint64_t max_ = 0;
};
