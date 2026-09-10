
#pragma once
#include <oneapi/tbb/parallel_for.h>
#include <vector>

#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/core/seq/concepts.hpp"

/**
 * Count number of terminals in each `SequenceContainer` and accumulate output
 * positions.
 */
template <countable_container T>
std::vector<std::size_t>
planTerminalRanges(std::span<T> data_, std::size_t k_) {
  std::vector<std::size_t> out(data_.size());

  oneapi::tbb::parallel_for((std::size_t)0, data_.size(), (std::size_t)1,
      [&](std::size_t i_) { out[i_] = data_[i_]->numTerminals(k_); });

  std::partial_sum(out.cbegin(), out.cend(), out.begin());

  return out;
}

/**
 * Extract all terminal sequences into a `TerminalBuffer`. Result is sorted
 * and made unique.
 */
template <sequence_container_like T>
TerminalBuffer
extractTerminalsSparse(std::span<T> data_, std::size_t k_) {
  auto blocks = planTerminalRanges(data_, k_);

  TerminalBuffer buffer(blocks.back(), k_, TerminalBuffer::autofit_tag);
  buffer.fill(data_, blocks);
}

/**
 * Extract unique terminal sequences into a `TerminalBuffer` by
 * competitive multi-threaded insertion of sequences and an atomic
 * pressence-absence pre-insertion check. Results is sorted and
 * made unique.
 */
TerminalBuffer
extractTerminalsDense(const std::vector<const SequenceContainer *> &data_,
                      std::size_t s_);
