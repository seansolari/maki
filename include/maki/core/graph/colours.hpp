
#pragma once
#include <cstdint>
#include <string>
#include <vector>

using colour_t = uint32_t;

struct ColourRegistry {
  std::vector<std::string> seeds; // seed names
  std::vector<uint64_t> occs; // occurrences of each seed
  std::vector<std::vector<colour_t>> metas; // meta colours
};
