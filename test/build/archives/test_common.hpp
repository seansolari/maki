
#pragma once
#include "maki/build/graph/archive/array_builder.hpp"
#include "maki/core/graph/archive/array.hpp"
#include "maki/core/utils/tempfile.hpp"
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <random>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

inline std::string MakeTempPath(const char *ext) {
  return tempio::create_temporary_file(std::filesystem::temp_directory_path(),
                                       ext);
}

inline std::vector<uint8_t> ReadFile(const std::string &path) {
  std::ifstream ifs(path, std::ios::binary);
  if (!ifs)
    throw std::runtime_error("ReadFile open failed: " + path);
  return std::vector<uint8_t>((std::istreambuf_iterator<char>(ifs)),
                              std::istreambuf_iterator<char>());
}

inline void WriteFile(const std::string &path,
                      const std::vector<uint8_t> &bytes) {
  std::ofstream ofs(path, std::ios::binary);
  if (!ofs)
    throw std::runtime_error("WriteFile open failed: " + path);
  ofs.write(reinterpret_cast<const char *>(bytes.data()), bytes.size());
  ofs.flush();
}

inline ArrayBuilder GenValues(size_t n, uint64_t mask) {
  ArrayBuilder v;
  v.reserve(n);
  std::mt19937_64 rng(12345);
  for (size_t i = 0; i < n; ++i)
    v.push(rng() & mask);
  return v;
}

template <typename Int = uint64_t>
inline std::vector<Int> GenVector(size_t n, Int mask) {
  std::vector<Int> v(n);
  std::mt19937_64 rng(12345);
  for (size_t i = 0; i < n; ++i)
    v[i] = static_cast<Int>(rng()) & mask;
  return v;
}

struct PackedResult {
  std::vector<uint8_t> data;
  uint64_t elem_count;
  uint8_t bit_width;

  RawPackedArray view() const {
    return RawPackedArray{data.data(), data.size(), elem_count, bit_width};
  }
};

inline PackedResult PackRaw(const ArrayBuilder &arr) {
  PackedResult res;
  arr.finalize_packed(res.data);
  res.elem_count = arr.size();
  res.bit_width = arr.required_bit_width();
  return res;
}