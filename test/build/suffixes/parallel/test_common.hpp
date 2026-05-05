#pragma once

#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include <array>
#include <cstdint>

constexpr uint8_t A = 0;
constexpr uint8_t C = 1;
constexpr uint8_t G = 2;
constexpr uint8_t T = 3;

static constexpr std::array<uint8_t, 5> AllEdges{A, C, G, T, TerminalEdge};

static inline LongSuffix makeLongSuffixFromDNA(std::string_view dna) {
  LongSuffix sfx(dna.size());
  for (const auto &c : dna) {
    switch (c) {
    case 'A':
      sfx.push(0);
      break;
    case 'C':
      sfx.push(1);
      break;
    case 'G':
      sfx.push(2);
      break;
    case 'T':
      sfx.push(3);
      break;
    }
  }
  return sfx;
}

static inline void enumerateDNA(size_t maxLen, std::vector<std::string> &out,
                                std::string current = "") {
  out.push_back(current); // includes empty suffix

  if (current.size() == maxLen)
    return;

  static const char *bases = "ACGT";
  for (int i = 0; i < 4; ++i) {
    enumerateDNA(maxLen, out, current + bases[i]);
  }
}
