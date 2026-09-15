
#pragma once
#include <oneapi/tbb/parallel_for.h>
#include <vector>

#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/core/seq/concepts.hpp"

/**
 * Count number of terminals in each `SequenceContainer` and accumulate output
 * positions.
 */
std::vector<std::size_t> planTerminalRanges(const container_span auto &data_,
                                            std::size_t k_) {
  std::vector<std::size_t> out(data_.size());

  oneapi::tbb::parallel_for(
      (std::size_t)0, data_.size(), (std::size_t)1,
      [&](std::size_t i_) { out[i_] = data_[i_].numTerminals(k_); });

  std::partial_sum(out.cbegin(), out.cend(), out.begin());

  return out;
}

/**
 * Extract all terminal sequences into a `TerminalBuffer`. Result is sorted
 * and made unique.
 */
TerminalBuffer extractTerminalsSparse(const container_span auto &data_,
                                      std::size_t k_) {
  auto blocks = planTerminalRanges(data_, k_);

  TerminalBuffer buffer(blocks.back(), k_, TerminalBuffer::autofit_tag);
  buffer.fillTerminals(data_, blocks);

  return buffer.OOPsort();
}

/**
 * Extract unique terminal sequences into a `TerminalBuffer` by
 * competitive multi-threaded insertion of sequences and an atomic
 * pressence-absence pre-insertion check. Results is sorted and
 * contains unique values.
 *
 * Preallocate maximum possible space for a given suffix size, and insert
 * terminal sequences using a lock-free gating system that checks if a
 * given terminal has been observed before.
 */
TerminalBuffer extractTerminalsDense(const container_span auto &data_,
                                     std::size_t s_) {
  assert(s_ > 0);
  TerminalBuffer buffer(5 * ShortSuffix::numSuffixes(s_ - 1), s_,
                        TerminalBuffer::autofit_tag);

  LongSuffixGate lock(s_);
  buffer.fillTerminals(data_, lock);
  buffer.shrink(lock.count());

  return buffer.OOPsort(false);
}
