
#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

#include <sdsl/int_vector.hpp>

#include "archive/archive_writer.hpp"
#include "maki/build/kmers/buffers/kmers.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/core/graph/archive/counts.hpp"
#include "maki/core/seq/io.hpp"

struct value_comp {
  inline constexpr bool operator()(const BufferValue &lhs,
                                   const BufferValue &rhs) const {
    return ((lhs.edge() << 61) | lhs.colour()) <
           ((rhs.edge() << 61) | rhs.colour());
  }
};

struct edge_comp {
  inline constexpr bool operator()(const BufferValue &lhs,
                                   const BufferValue &rhs) const {
    return lhs.edge() < rhs.edge();
  }
};

struct edge_eq {
  inline constexpr bool operator()(const BufferValue &lhs,
                                   const BufferValue &rhs) const {
    return lhs.edge() == rhs.edge();
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
              uint8_t msb_dna5, ArchivePayload &carch, MetaColours &cmap);

/**
 * Push a node into the graph buffers, comprising graph
 * structure (edges, succ) as well as edge count data.
 */
void pushNode(packet &pkt, sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
              uint8_t msb_dna5, CountBuffer &c);

/**
 * Push a node into the graph buffers, comprising graph
 * structure (edges, succ).
 */
void pushNode(packet &pkt, sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
              uint8_t msb_dna5);

/**
 * Iterate over k-mers and push structure into graph buffers.
 */
template <typename It, class... Args>
void pushNodes(packet &pkt, It it, It end, sdsl::int_vector<2>::iterator b,
               sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
               uint8_t msb_dna4, Args &...args) {
  while (it != end) {
    if (*b == BW_0_K)
      ++pkt.block;

    do {
      pkt.emplace(it.readValue());
      ++it;
      ++b;
    } while ((it != end) && (*b == IS_0));

    pushNode(pkt, edges, succ, parsing::dna4ToDna5(msb_dna4), args...);
  }
}

/**
 * Push terminals to output.
 */
template <class Buffer, class... Args>
push_summary pushRange(Buffer &data, sdsl::int_vector<2>::iterator to,
                       sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
                       uint8_t msb_dna4, Args &...args) {
  packet pkt;
  pushNodes(pkt, data.begin(), data.end(), to, edges, succ, msb_dna4, args...);
  return pkt.str;
}

/**
 * Interleave k-mers and terminals (sized k-mers), pushing co-lex ordered
 * nodes into the graph buffers, encoding colour information with basic
 * bit-packing and encoding colour tuples.
 *
 * Returns counts for the number of nodes and edges pushed, for each
 * edge label.
 */
template <class... Args>
push_summary interleave(KmerBuffer &kmers, sdsl::int_vector<2>::iterator ko,
                        TerminalRange &terminals,
                        sdsl::int_vector<2>::iterator to,
                        sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
                        uint8_t msb_dna4, Args &...args) {
  std::size_t k = kmers.getK(), km1 = k - 1, keff = kmers.getEffK();
  std::size_t prevTerminalSize = 0;
  assert(km1 > prevTerminalSize);

  MaskedBytesDiff Diff(keff);
  MaskedAlignedBytesLessThan<RHS> LessThan(keff, terminals.lengthBytes());

  auto km_it = kmers.begin(), km_end = kmers.end();
  auto tm_it = terminals.constBegin(), tm_end = terminals.constEnd();

  packet pkt;
  while ((km_it != km_end) || (tm_it != tm_end)) {
    // insert terminals (they always start a new block)

    auto tm_pivot = km_it == km_end
                        ? tm_end
                        : std::upper_bound(tm_it, tm_end, *km_it, LessThan);

    if (tm_it != tm_pivot) {
      pushNodes(pkt, tm_it, tm_pivot, to, edges, succ, msb_dna4, args...);
      size_t insertedTerminals = tm_pivot - tm_it;
      assert(insertedTerminals > 0);
      to += insertedTerminals;

      // capture last terminal size, used for
      // overlap correction when comparing against k-mers
      auto tm_penul = tm_it + insertedTerminals - 1;
      prevTerminalSize = tm_penul.readSize();
      tm_it = tm_pivot;
    }

    // insert k-mers
    auto km_pivot = tm_it == tm_end
                        ? km_end
                        : std::lower_bound(km_it, km_end, *tm_it, LessThan);

    if (km_it != km_pivot) {
      // correct border overlap value
      if (prevTerminalSize == km1) {
        // if their bytes are equal and the terminal size is `k-1`, then the
        // only difference is the final `$` char

        auto dtk = Diff(std::to_address(/* guaranteed */ tm_it - 1) +
                            tm_it.lengthBytes(),
                        std::to_address(km_it));
        *ko = dtk == IS_0 ? IS_K : dtk;
      }

      // insert nodes
      pushNodes(pkt, km_it, km_pivot, ko, edges, succ, msb_dna4, args...);
      size_t insertedKmers = km_pivot - km_it;
      assert(insertedKmers > 0);
      ko += insertedKmers;
      km_it = km_pivot;
    }
  }

  return pkt.str;
}
