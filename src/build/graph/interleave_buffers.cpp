
#include "maki/build/graph/interleave_buffers.hpp"
#include "maki/build/kmers/buffers/kmers.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"

#include <algorithm>
#include <cassert>
#include <memory>

void pushNode(packet &pkt, sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
              ArchivePayload &carch, MetaColours &cmap, uint8_t msb) {
  if (pkt.data.empty())
    return;

  // sort temp data, but don't make it unique
  std::sort(pkt.data.begin(), pkt.data.end(), value_comp{});
  auto it = pkt.data.begin(), end = std::unique(it, pkt.data.end());
  // if there are valid edges other than `$`, then `$` is not required
  {
    auto valid = std::find_if(
        it, end, [](const BufferValue &v) { return v.edge() != 0b000; });
    if (valid != end)
      it = valid;
  }

  // iterate through each edge
  while (it != end) {
    // collect all colours for edge
    uint64_t edge = it->edge();
    do {
      pkt.colours.push_back(it->colour());
      ++it;
    } while ((it != end) && (it->edge() == edge));
    // insert new edge
    ++pkt.str.F[msb];
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
  ++pkt.str.C[msb];
  pkt.data.clear();
}

/**
 * Interleave k-mers and terminals (sized k-mers), pushing colex ordered
 * nodes into the graph buffers, encoding colour information with basic
 * bit-packing and encoding colour tuples.
 *
 * Returns counts for the number of nodes and edges pushed, for each
 * edge label.
 */
push_summary interleave(KmerBuffer &kmers, sdsl::int_vector<2>::iterator ko,
                        TerminalRange &terminals,
                        sdsl::int_vector<2>::iterator to,
                        sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
                        ArchivePayload &carch, MetaColours &cmap, uint8_t msb) {
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
      pushNodes(pkt, tm_it, tm_pivot, to, edges, succ, carch, cmap, msb);
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
      pushNodes(pkt, km_it, km_pivot, ko, edges, succ, carch, cmap, msb);
      size_t insertedKmers = km_pivot - km_it;
      assert(insertedKmers > 0);
      ko += insertedKmers;
      km_it = km_pivot;
    }
  }

  return pkt.str;
}

push_summary pushRange(TerminalRange &terminals,
                       sdsl::int_vector<2>::iterator to,
                       sdsl::int_vector<4> &edges, sdsl::bit_vector &succ,
                       ArchivePayload &carch, MetaColours &cmap, uint8_t msb) {
  packet pkt;
  pushNodes(pkt, terminals.begin(), terminals.end(), to, edges, succ, carch,
            cmap, msb);
  return pkt.str;
}
