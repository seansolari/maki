
#include "maki/build/graph/interleave_buffers.hpp"

#include <algorithm>

void pushNode(packet &pkt, sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
              uint8_t msb_dna5, ArchivePayload &carch, MetaColours &cmap) {
  if (pkt.data.empty()) {
    return;
  }

  // sort temp data and make it unique
  std::sort(pkt.data.begin(), pkt.data.end(), value_comp{});
  auto it = pkt.data.begin(), end = std::unique(it, pkt.data.end());
  // if there are valid edges other than `$`, then `$` is not required
  {
    auto valid = std::find_if(
        it, end, [](const BufferValue &v) { return v.edge() != 0b000; });
    if (valid != end)
      it = valid;
  }

  while (it != end) {
    // collect all colours for edge
    uint64_t edge = it->edge();
    do {
      pkt.colours.push_back(it->colour());
      ++it;
    } while ((it != end) && (it->edge() == edge));
    // insert new edge
    ++pkt.str.F[msb_dna5];
    if (pkt.last[edge] != pkt.block) {
      edges.push_back(edge | 0b1000);
      pkt.last[edge] = pkt.block;
    } else {
      edges.push_back(edge);
    }
    // attach edge to node
    succ.push_back(0);
    // push colour
    auto ccode = cmap.insert(pkt.colours);
    carch.raw.push(ccode);
    pkt.colours.clear();
  }

  // finalise node
  succ.back() = 1;
  ++pkt.str.C[msb_dna5];
  pkt.data.clear();
}

void pushNode(packet &pkt, sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
              uint8_t msb_dna5, CountBuffer &c) {
  if (pkt.data.empty()) {
    return;
  }

  // sort temp values but do not remove duplicates
  std::sort(pkt.data.begin(), pkt.data.end(), value_comp{});
  auto it = pkt.data.begin(), end = pkt.data.end();
  // if there are valid edges other than `$`, then `$` is not required
  {
    auto valid = std::find_if(
        it, end, [](const BufferValue &v) { return v.edge() != 0b000; });
    if (valid != end)
      it = valid;
  }

  // count occurrence of each edge
  std::array<uint64_t, 5> counts = {0, 0, 0, 0, 0};
  uint8_t outdegree = 0;
  while (it != end) {
    uint64_t edge = it->edge();
    do {
      ++counts[outdegree];
      ++it;
    } while ((it != end) && (it->edge() == edge));
    ++outdegree;
    // insert structural elements for this edge
    ++pkt.str.F[msb_dna5];
    if (pkt.last[edge] != pkt.block) {
      edges.push_back(edge | 0b1000);
      pkt.last[edge] = pkt.block;
    } else {
      edges.push_back(edge);
    }
    // attach edge to node
    succ.push_back(0);
  }

  // encode edge counts for this node
  c.insert_node({.outdegree = outdegree, .edge_counts = counts.data()});
  // finalise node
  succ.back() = 1;
  ++pkt.str.C[msb_dna5];
  pkt.data.clear();
}
