
#pragma once
#include <vector>

#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/core/seq/concepts.hpp"

/**
 * Count number of terminals in each `SequenceContainer` and accumulate output
 * positions.
 */
std::vector<std::size_t>
planTerminalRanges(const std::vector<const SequenceContainer *> &data_,
                   std::size_t k_);

/**
 * Extract all terminal sequences into a `TerminalBuffer`. Result is sorted
 * and made unique.
 */
TerminalBuffer
extractTerminalsSparse(const std::vector<const SequenceContainer *> &data_,
                       std::size_t k_);

/**
 * Extract unique terminal sequences into a `TerminalBuffer` by
 * competitive multi-threaded insertion of sequences and an atomic
 * pressence-absence pre-insertion check. Results is sorted and
 * made unique.
 */
TerminalBuffer
extractTerminalsDense(const std::vector<const SequenceContainer *> &data_,
                      std::size_t s_);
