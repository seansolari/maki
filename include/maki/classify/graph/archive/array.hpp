
#pragma once
#include <cstdint>
#include <stdexcept>

class RawPackedArray
{
public:
  RawPackedArray(const uint8_t *ptr,
                 size_t len,
                 uint64_t elem_count,
                 uint8_t bit_width)
      : base_(ptr), len_(len), n_(elem_count), w_(bit_width)
  {
    if (!ptr)
      throw std::invalid_argument("null ptr");
  }

  uint64_t size() const { return n_; }
  uint64_t bit_width() const { return w_; }

  uint64_t get(uint64_t i) const
  {
    if (i >= n_)
      throw std::out_of_range("index");
    __uint128_t bitpos = (__uint128_t)i * w_;
    size_t byteidx = (size_t)(bitpos >> 3);
    unsigned b = (unsigned)(bitpos & 7);

    __uint128_t lane = 0;
    const size_t avail = len_ - (byteidx < len_ ? byteidx : len_);
    const size_t take = avail >= 16 ? 16 : avail;

    for (size_t j = 0; j < take; ++j)
      lane |= (__uint128_t)base_[byteidx + j] << (8 * j);

    lane >>= b;

    if (w_ == 64)
      return (uint64_t)lane;
    uint64_t mask = ((uint64_t)1 << w_) - 1;
    return (uint64_t)lane & mask;
  }

private:
  const uint8_t *base_;
  size_t len_;
  uint64_t n_;
  uint8_t w_;
};
