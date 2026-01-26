
#pragma once
#include <cstdint>
#include <string>
#include <vector>

struct ColourRegistry {
  std::vector<std::string> seeds; // seed names
  std::vector<uint64_t> occs; // occurrences of each seed
  std::vector<std::vector<uint64_t>> metas; // meta colours
};
