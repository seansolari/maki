
#pragma once
#include <cassert>
#include <cstdint>
#include <functional>
#include <string>
#include <variant>
#include <vector>

using colour_t = uint32_t;

struct ColourRegistry {
  using vector_ref = std::reference_wrapper<const std::vector<colour_t>>;
  using id_result = std::variant<colour_t, vector_ref>;

  // Get seed name for a colour id.
  inline const std::string &seed(colour_t c) const {
    assert(c < seeds.size());
    return seeds[c];
  }

  // prevalence of colour
  inline uint64_t count(colour_t c) const {
    assert(c < occs.size());
    return occs[c];
  }

  // Get colours for a meta colour id
  inline id_result colours(uint64_t id) const {
    if (id < seeds.size()) {
      return static_cast<colour_t>(id);
    } else {
      return std::cref(metas[id - seeds.size()]);
    }
  }

  std::vector<std::string> seeds;           // seed names
  std::vector<uint64_t> occs;               // occurrences of each seed
  std::vector<std::vector<colour_t>> metas; // meta colours
};
