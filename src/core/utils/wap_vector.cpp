
#include "maki/core/utils/wap_vector.hpp"
#include "maki/core/utils/algo.hpp"
#include <limits>
#include <sdsl/io.hpp>

wo_wap_vector::wo_wap_vector(std::size_t reserve_) {
  is_u8.reserve(reserve_);
  is_u16.reserve(reserve_);
  is_u32.reserve(reserve_);
  is_u64.reserve(reserve_);
}

wo_wap_vector::wo_wap_vector(const std::vector<uint64_t> &data)
    : wo_wap_vector::wo_wap_vector(data.size()) {
  for (const auto &x : data)
    push_back(x);
}

std::size_t wo_wap_vector::size() const {
  assert(is_u16.size() == is_u8.size());
  assert(is_u32.size() == is_u8.size());
  assert(is_u64.size() == is_u8.size());
  return is_u8.size();
}

void wo_wap_vector::push_back(uint64_t value_) {
  is_u8.push_back(0);
  is_u16.push_back(0);
  is_u32.push_back(0);
  is_u64.push_back(0);

  if (value_ <= std::numeric_limits<uint8_t>::max()) {
    is_u8.back() = 1;
    u8.push_back(value_);
  } else if (value_ <= std::numeric_limits<uint16_t>::max()) {
    is_u16.back() = 1;
    u16.push_back(value_);
  } else if (value_ <= std::numeric_limits<uint32_t>::max()) {
    is_u32.back() = 1;
    u32.push_back(value_);
  } else {
    is_u64.back() = 1;
    u64.push_back(value_);
  }
}

void wo_wap_vector::clear() noexcept {
  is_u8.clear();
  is_u16.clear();
  is_u32.clear();
  is_u64.clear();
  u8.clear();
  u16.clear();
  u32.clear();
  u64.clear();
}

void wo_wap_vector::append(const wo_wap_vector &rhs) {
  vector_append(is_u8, rhs.is_u8);
  vector_append(is_u16, rhs.is_u16);
  vector_append(is_u32, rhs.is_u32);
  vector_append(is_u64, rhs.is_u64);
  vector_append(u8, rhs.u8);
  vector_append(u16, rhs.u16);
  vector_append(u32, rhs.u32);
  vector_append(u64, rhs.u64);
}

ro_wap_vector::ro_wap_vector(const wo_wap_vector &v)
    : is_u8(v.is_u8), is_u16(v.is_u16), is_u32(v.is_u32), is_u64(v.is_u64),
      rs_u8(), rs_u16(), rs_u32(), rs_u64(), u8(v.u8), u16(v.u16), u32(v.u32),
      u64(v.u64) {
  init_support();
}

ro_wap_vector::ro_wap_vector(wo_wap_vector &&v)
    : is_u8(v.is_u8), is_u16(v.is_u16), is_u32(v.is_u32), is_u64(v.is_u64),
      rs_u8(), rs_u16(), rs_u32(), rs_u64(), u8(std::move(v.u8)),
      u16(std::move(v.u16)), u32(std::move(v.u32)), u64(std::move(v.u64)) {
  init_support();
}

uint64_t ro_wap_vector::operator[](std::size_t i) const {
  if (is_u8[i])
    return u8[rs_u8(i)];
  if (is_u16[i])
    return u16[rs_u16(i)];
  if (is_u32[i])
    return u32[rs_u32(i)];
  if (is_u64[i])
    return u64[rs_u64(i)];
  return 0;
}

std::size_t ro_wap_vector::size() const {
  assert(is_u16.size() == is_u8.size());
  assert(is_u32.size() == is_u8.size());
  assert(is_u64.size() == is_u8.size());
  return is_u8.size();
}

void ro_wap_vector::init_support() {
  sdsl::util::init_support(rs_u8, &is_u8);
  sdsl::util::init_support(rs_u16, &is_u16);
  sdsl::util::init_support(rs_u32, &is_u32);
  sdsl::util::init_support(rs_u64, &is_u64);
}

void ro_wap_vector::stabilize_support() {
  rs_u8.set_vector(&is_u8);
  rs_u16.set_vector(&is_u16);
  rs_u32.set_vector(&is_u32);
  rs_u64.set_vector(&is_u64);
}

std::size_t detail::size_in_bytes(const wo_wap_vector &vec) {
  return sdsl::size_in_bytes(vec.is_u8) + sdsl::size_in_bytes(vec.is_u16) +
         sdsl::size_in_bytes(vec.is_u32) + sdsl::size_in_bytes(vec.is_u64) +
         (sizeof(uint8_t) * vec.u8.capacity()) +
         (sizeof(uint16_t) * vec.u16.capacity()) +
         (sizeof(uint32_t) * vec.u32.capacity()) +
         (sizeof(uint64_t) * vec.u64.capacity());
}

std::size_t detail::size_in_bytes(const ro_wap_vector &vec) {
  return sdsl::size_in_bytes(vec.is_u8) + sdsl::size_in_bytes(vec.is_u16) +
         sdsl::size_in_bytes(vec.is_u32) + sdsl::size_in_bytes(vec.is_u64) +
         sdsl::size_in_bytes(vec.rs_u8) + sdsl::size_in_bytes(vec.rs_u16) +
         sdsl::size_in_bytes(vec.rs_u32) + sdsl::size_in_bytes(vec.rs_u64) +
         (sizeof(uint8_t) * vec.u8.capacity()) +
         (sizeof(uint16_t) * vec.u16.capacity()) +
         (sizeof(uint32_t) * vec.u32.capacity()) +
         (sizeof(uint64_t) * vec.u64.capacity());
}
