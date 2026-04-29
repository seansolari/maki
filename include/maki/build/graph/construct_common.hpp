#pragma once

#include <cstddef>
#include <filesystem>

namespace dbg {

struct BuildOptions {
  // Algorithm parameters
  std::size_t kmer_size = 31;
  std::size_t suffix_size = 8;
  // Output parameters
  std::filesystem::path out;
  // Space parameters
  std::size_t pool_size = 16;
  std::size_t reserve_per_chunk = 0;
  std::size_t chunks() const;
};

} // namespace dbg
