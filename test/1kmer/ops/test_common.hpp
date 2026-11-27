#pragma once

#include "maki/build/kmers/buffers/terminals.hpp"

inline LongSuffix makeTerminal(uint8_t k, std::string_view seq) {
  LongSuffix t(k);

  for (char c : seq) {
    t.push(seqan3::to_rank(seqan3::dna4{}.assign_char(c)));
  }

  return t;
}

inline void writeRecord(TerminalBuffer::iterator it, LongSuffix const &t,
                 unsigned char edge_) {
  it.writeKey(t);
  it.writeEdge(edge_);
}

inline std::vector<std::string> exportKmers(const TerminalBuffer &buffer) {
  std::vector<std::string> expected;

  for (auto it = buffer.begin(), end = buffer.end(); it != end; ++it) {
    expected.push_back(toString(it.readKmer()));
  }

  return expected;
}

inline std::vector<unsigned char> exportEdges(const TerminalBuffer &buffer) {
  std::vector<unsigned char> expected;

  for (auto it = buffer.begin(), end = buffer.end(); it != end; ++it) {
    expected.push_back(it.readEdge());
  }

  return expected;
}
