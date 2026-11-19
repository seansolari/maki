#pragma once

#include "maki/core/seq/concepts.hpp"
#include "maki/build/io/gff.hpp"
#include <stdexcept>
#include <string>
#include <vector>

inline std::string revcomp(std::string const & seq)
{
    std::string result;
    result.reserve(seq.size());

    for (auto it = seq.rbegin(); it != seq.rend(); ++it)
    {
        switch (*it)
        {
            case 'A': result.push_back('T'); break;
            case 'C': result.push_back('G'); break;
            case 'G': result.push_back('C'); break;
            case 'T': result.push_back('A'); break;
            default:
                throw std::runtime_error(
                    "Invalid dna4 character");
        }
    }

    return result;
}

std::vector<std::string> collectSequences(sequence_fragment_input_range auto &&range) {
  std::vector<std::string> result;

  for (const auto &fragment : range)
    result.push_back(toString(Dna4Sequence(fragment.data())));

  return result;
}

std::vector<uint64_t> collectIds(sequence_fragment_input_range auto &&range) {
  std::vector<uint64_t> ids;

  for (auto const &fragment : range)
    ids.push_back(fragment.id());

  return ids;
}
