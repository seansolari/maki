
#pragma once

#include "maki/core/graph/cdbg.hpp"
#include "maki/core/graph/colours.hpp"
#include "maki/core/graph/wdbg.hpp"
#include <cassert>
#include <cstdint>
#include <vector>

using id_result = ColourRegistry::id_result;

class ClassificationSummary {
public:
  void count(std::size_t edge, uint64_t depth, id_result colours);
};

ClassificationSummary Accumulate(const std::vector<int64_t> &edges,
                                 const WeightedGraph &qry,
                                 const ColouredGraph &ref);
