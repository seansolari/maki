
#include "maki/classify/summary.hpp"
#include <stdexcept>

ClassificationSummary Accumulate(const std::vector<int64_t> &edges,
                                 const WeightedGraph &qry,
                                 const ColouredGraph &ref) {
  ClassificationSummary results;

  std::size_t i = 0, numEdges = edges.size();
  while (i < numEdges && edges[i] == -1)
    ++i;

  while (i < numEdges) {
    // get region of colour archive
    auto region = ref.carch->find_chunk_by_global_index(edges[i]);
    if (region == -1) {
      throw std::runtime_error("Reached invalid reference edge: " +
                               std::to_string(edges[i]));
    }

    // get colour array
    auto cols = ref.carch->view(region);

    // get range of edge colours covered by this array
    const auto &meta = ref.carch->meta(region);
    auto regionEnd = static_cast<int64_t>(meta.start_index + meta.elem_count);
    assert(static_cast<int64_t>(meta.start_index) <= edges[i]);
    assert(edges[i] < regionEnd);

    // count through edges in this range
    while (i < numEdges) {
      if (edges[i] != -1) {
        if (edges[i] < regionEnd) {
          results.count(
              i, qry.edge_count(i),
              ref.cmap.colours(cols.get(edges[i] - meta.start_index)));
        } else {
          break;
        }
      }

      ++i;
    }
  }

  return results;
}
