
#include "maki/build/kmers/construct_terminals.hpp"

#include <numeric>
#include "oneapi/tbb/parallel_for.h"

std::vector<std::size_t> planTerminalRanges(const std::vector<const SequenceContainer*> &data_, std::size_t k_) {
  std::vector<std::size_t> out(data_.size());
  oneapi::tbb::parallel_for((std::size_t)0, data_.size(), (std::size_t)1, [&](std::size_t i_) {
    out[i_] = data_[i_]->numTerminals(k_);
  });
  std::partial_sum(out.cbegin(), out.cend(), out.begin());
  return out;
}

TerminalBuffer extractTerminalsSparse(const std::vector<const SequenceContainer*> &data_, std::size_t k_) {
  auto blocks = planTerminalRanges(data_, k_);
  TerminalBuffer buffer(blocks.back(), k_, TerminalBuffer::autofit_tag);
  buffer.fill(data_, blocks);
  return buffer;
}
