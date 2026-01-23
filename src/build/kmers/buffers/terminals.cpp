
#include "maki/build/kmers/buffers/terminals.hpp"

namespace {

void ser(std::size_t src, uint8_t *dst, uint32_t bytes) {
  static constexpr std::size_t msk = 0xFF;
  for (std::size_t i = 0; i < bytes; ++i)
    *(dst + i) = static_cast<uint8_t>((src & (msk << (8 * i))) >> (8 * i));
}

std::size_t deser(uint8_t const *src, uint32_t bytes) {
  std::size_t dest = 0;
  for (std::size_t i = 0; i < bytes; ++i)
    dest |= static_cast<std::size_t>(*(src + i)) << (i * 8);
  return dest;
}

} // namespace
