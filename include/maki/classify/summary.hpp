
#pragma once

#include "maki/core/graph/cdbg.hpp"
#include "maki/core/graph/colours.hpp"
#include "maki/core/graph/wdbg.hpp"
#include <cassert>
#include <cstdint>
#include <vector>

using id_result = ColourRegistry::id_result;

class ClassificationSummary {
  std::unordered_map<std::pair<std::size_t,colour_t>,uint32_t> _data;

public:
  void count(std::size_t edge, uint64_t depth, id_result colours);
};

ClassificationSummary Accumulate(const std::vector<int64_t> &edges,
                                 const WeightedGraph &qry,
                                 const ColouredGraph &ref);
