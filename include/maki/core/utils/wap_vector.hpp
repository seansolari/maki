
#pragma once
#include <cassert>
#include <cstdint>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rrr_vector.hpp>
#include <vector>

class wo_wap_vector;
class ro_wap_vector;

namespace detail {

std::size_t size_in_bytes(const wo_wap_vector &);
std::size_t size_in_bytes(const ro_wap_vector &);

} // namespace detail

// Width-adaptive packing vector
class wo_wap_vector {
  friend class ro_wap_vector;
  friend std::size_t detail::size_in_bytes(const wo_wap_vector &);

public:
  wo_wap_vector() = default;
  explicit wo_wap_vector(std::size_t reserve_);
  wo_wap_vector(const std::vector<uint64_t> &);
  std::size_t size() const;
  void push_back(uint64_t value_);
  void clear() noexcept;
  void append(const wo_wap_vector &);

private:
  sdsl::bit_vector is_u8, is_u16, is_u32, is_u64;
  std::vector<uint8_t> u8;
  std::vector<uint16_t> u16;
  std::vector<uint32_t> u32;
  std::vector<uint64_t> u64;
};

class ro_wap_vector {
  friend std::size_t detail::size_in_bytes(const ro_wap_vector &);

public:
  ro_wap_vector(const wo_wap_vector &v);
  ro_wap_vector(wo_wap_vector &&v);
  uint64_t operator[](std::size_t i) const;
  std::size_t size() const;

protected:
  sdsl::rrr_vector<> is_u8, is_u16, is_u32, is_u64;
  sdsl::rrr_vector<>::rank_1_type rs_u8, rs_u16, rs_u32, rs_u64;
  std::vector<uint8_t> u8;
  std::vector<uint16_t> u16;
  std::vector<uint32_t> u32;
  std::vector<uint64_t> u64;

  void init_support();
};
