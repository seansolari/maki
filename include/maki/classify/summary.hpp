
#pragma once

#include "maki/core/graph/cdbg.hpp"
#include "maki/core/graph/colours.hpp"
#include "maki/core/graph/wdbg.hpp"
#include <cassert>
#include <cstdint>
#include <tuple>
#include <vector>

using vector_ref = ColourRegistry::vector_ref;
using id_result = ColourRegistry::id_result;

template<class... Ts>
struct overloaded : Ts... { using Ts::operator()...; };

class ClassificationSummary {
  std::vector<std::tuple<std::size_t,colour_t,uint64_t>> _data;

public:
  void count(std::size_t edge, uint64_t depth, id_result colours);
};

ClassificationSummary Accumulate(const std::vector<int64_t> &edges,
                                 const WeightedGraph &qry,
                                 const ColouredGraph &ref);
