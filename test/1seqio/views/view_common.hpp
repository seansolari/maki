#pragma once

#include <stdexcept>
#include <string>

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
