
#pragma once
#include <vector>

#include "./buffers/terminals.hpp"
#include "maki/core/seq/seq_concepts.hpp"

/**
 * Count number of terminals in each `SequenceContainer` and accumulate output
 * positions.
 */
std::vector<std::size_t>
planTerminalRanges(const std::vector<const SequenceContainer *> &data_,
                   std::size_t k_);

/**
 * Extract all terminal sequences into a `TerminalBuffer`. Result is neither
 * unique nor sorted.
 */
TerminalBuffer
extractTerminalsSparse(const std::vector<const SequenceContainer *> &data_,
                       std::size_t k_);

/**
 * Extract unique terminal sequences into a `TerminalBuffer` by
 * competetive multi-threaded insertion of sequences and an atomic
 * pressence-absence pre-insertion check. Results are unique but not sorted.
 */
TerminalBuffer
extractTerminalsDense(const std::vector<const SequenceContainer *> &data_,
                      std::size_t k_);
