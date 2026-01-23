
#pragma once
#include <cstdint>
#include <bit>

template <typename T>
inline constexpr uint8_t floor_log2(T __x)
{
  if (__x == 0)
    return 0;
  return uint8_t(sizeof(T) * 8 - 1) - static_cast<uint8_t>(std::countl_zero(__x));
}

// number of bits required to represent a number
template <typename T>
inline constexpr uint8_t required_bits(T __x)
{
  return floor_log2(__x) + 1u;
}

template <typename T>
inline constexpr uint8_t ceil_log2(T __x)
{
  if (__x == 0)
    return 0;
  return floor_log2(__x - 1u) + 1u;
}

// number of bytes required to represent a number
template <typename T>
inline constexpr uint8_t required_bytes(T __x)
{
  return (ceil_log2(__x) + 7) / 8;
}

template <typename T>
inline constexpr long double floor_log2_l(T __x)
{
  return static_cast<long double>(floor_log2(__x));
}

template <typename T>
inline constexpr long double ceil_log2_l(T __x)
{
  return static_cast<long double>(ceil_log2(__x));
}

// number of bytes per key
template <uint8_t bitsPerNucl = 2>
inline constexpr uint8_t key_size(uint8_t k)
{
  return ((bitsPerNucl * k) + 7) / 8;
}

// number of bytes per value
inline constexpr uint8_t value_size(uint64_t colourBits)
{
  return (colourBits + /* 3 bits per edge */ 3u + 7u) / 8u;
}

// number of bytes per record
inline constexpr uint32_t record_size(uint8_t keyBytes, uint8_t valueBytes)
{
  return static_cast<uint32_t>(keyBytes + valueBytes);
}

inline void ser(std::size_t src, uint8_t *dst, uint32_t bytes) {
  static constexpr std::size_t msk = 0xFF;
  for (std::size_t i = 0; i < bytes; ++i)
    *(dst + i) = static_cast<uint8_t>((src & (msk << (8 * i))) >> (8 * i));
}

inline std::size_t deser(uint8_t const *src, uint32_t bytes) {
  std::size_t dest = 0;
  for (std::size_t i = 0; i < bytes; ++i)
    dest |= static_cast<std::size_t>(*(src + i)) << (i * 8);
  return dest;
}
