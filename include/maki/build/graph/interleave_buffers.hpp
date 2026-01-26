
#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

#include <sdsl/int_vector.hpp>

#include "archive/archive_writer.hpp"
#include "maki/build/kmers/buffers/kmers.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"

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
};

/**
 * Interleave k-mers with terminals, writing edges, succ and colours.
 */
push_summary interleave(KmerBuffer &, sdsl::int_vector<2>::iterator,
                        TerminalRange &, sdsl::int_vector<2>::iterator,
                        sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
                        ArchivePayload &carch, MetaColours &cmap, uint8_t msb);

/**
 * Push terminals to output.
 */
push_summary pushRange(TerminalRange &, sdsl::int_vector<2>::iterator,
                       sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
                       ArchivePayload &carch, MetaColours &cmap, uint8_t msb);
