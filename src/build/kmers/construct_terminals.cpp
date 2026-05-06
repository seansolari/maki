
#include "maki/build/kmers/construct_terminals.hpp"

#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "oneapi/tbb/parallel_for.h"
#include <numeric>

std::vector<std::size_t>
planTerminalRanges(const std::vector<const SequenceContainer *> &data_,
                   std::size_t k_) {
  std::vector<std::size_t> out(data_.size());
  oneapi::tbb::parallel_for(
      (std::size_t)0, data_.size(), (std::size_t)1,
      [&](std::size_t i_) { out[i_] = data_[i_]->numTerminals(k_); });
  std::partial_sum(out.cbegin(), out.cend(), out.begin());
  return out;
}

/**
 * Calculate space requried to store terminal sequences for each input
 * chunk, and then parallel extract sequences into disjoint regions of
 * a single buffer.
 */
TerminalBuffer
extractTerminalsSparse(const std::vector<const SequenceContainer *> &data_,
                       std::size_t k_) {
  auto blocks = planTerminalRanges(data_, k_);
  TerminalBuffer buffer(blocks.back(), k_, TerminalBuffer::autofit_tag);
  buffer.fill(data_, blocks);
  return buffer.OOPsort();
}

/**
 * Preallocate maximum possible space for a given suffix size, and insert
 * terminal sequences using a lock-free gating system that checks if a
 * given terminal has been observed before.
 */
TerminalBuffer
extractTerminalsDense(const std::vector<const SequenceContainer *> &data_,
                      std::size_t s_) {
  assert(s_ > 0);
  TerminalBuffer buffer(5 * ShortSuffix::numSuffixes(s_ - 1), s_,
                        TerminalBuffer::autofit_tag);
  LongSuffixGate lock(s_);
  buffer.fill(data_, lock);
  buffer.shrink(lock.count());
  return buffer.OOPsort(false);
}
