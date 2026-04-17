
#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

#include <sdsl/int_vector.hpp>

#include "archive/archive_writer.hpp"
#include "maki/build/kmers/buffers/kmers.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/core/seq/seq_io.hpp"

struct value_comp {
  inline constexpr bool operator()(const BufferValue &lhs,
                                   const BufferValue &rhs) const {
    return ((lhs.edge() << 61) | lhs.colour()) <
           ((rhs.edge() << 61) | rhs.colour());
  }
};

struct push_summary {
  std::array<std::size_t, 5> F = {0, 0, 0, 0, 0}, C = {0, 0, 0, 0, 0};
};

struct packet {
  push_summary str;
  std::vector<BufferValue> data;
  ColourVector colours;
  int64_t block = -1;
  std::array<int64_t, 5> last = {-1, -1, -1, -1, -1};

  // Push a coloured edge into the current node's buffer
  void emplace(BufferValue &&dna4) {
    data.emplace_back(dna4.colour(), parsing::dna4ToDna5(dna4.edge()));
  }
};

/**
 * Push a node into the graph buffers, comprising graph
 * structure (edges, succ) as well as colour data.
 */
void pushNode(packet &pkt, sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
              ArchivePayload &carch, MetaColours &cmap, uint8_t msb_dna5);

/**
 * Iterate over k-mers and push structure into graph buffers.
 */
template <typename It>
void pushNodes(packet &pkt, It it, It end, sdsl::int_vector<2>::iterator b,
               sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
               ArchivePayload &carch, MetaColours &cmap, uint8_t msb_dna4) {
  while (it != end) {
    if (*b == BW_0_K)
      ++pkt.block;

    do {
      pkt.emplace(it.readValue());
      ++it;
      ++b;
    } while ((it != end) && (*b == IS_0));

    pushNode(pkt, edges, succ, carch, cmap, parsing::dna4ToDna5(msb_dna4));
  }
}

/**
 * Interleave k-mers with terminals, writing edges, succ and colours.
 */
push_summary interleave(KmerBuffer &, sdsl::int_vector<2>::iterator,
                        TerminalRange &, sdsl::int_vector<2>::iterator,
                        sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
                        ArchivePayload &carch, MetaColours &cmap, uint8_t msb_dna4);

/**
 * Push terminals to output.
 */
push_summary pushRange(TerminalRange &, sdsl::int_vector<2>::iterator,
                       sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
                       ArchivePayload &carch, MetaColours &cmap, uint8_t msb_dna4);
