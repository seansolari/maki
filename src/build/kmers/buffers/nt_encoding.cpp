#include "maki/build/kmers/buffers/nt_encoding.hpp"

#include <iterator>

ShortSuffix::ShortSuffix(std::size_t _size, uint64_t init)
    : _s(_size), _data(init) {
  assert((sizeof(uint64_t) * 8) / 2 >= _s);
}
ShortSuffix::ShortSuffix(std::size_t _size) : ShortSuffix(_size, 0) {}
ShortSuffix::ShortSuffix(std::size_t _size, Dna4SequenceConstIter _it)
    : ShortSuffix(_size, 0) {
  for (std::size_t i = 0; i < _s; ++i)
    _data |= parsing::dna4ToLong(*_it++) << (2 * i);
}
ShortSuffix::ShortSuffix(Dna4SequenceConstIter begin, Dna4SequenceConstIter end)
    : ShortSuffix(std::distance(begin, end), 0) {
  for (std::size_t i = 0; i < _s; ++i)
    _data |= parsing::dna4ToLong(*begin++) << (2 * i);
}

std::string Kmer::toString() const {
  std::string sequence;
  sequence.reserve(_num_nts);
  for (size_t i = 0; i < _num_nts; ++i) {
    uint8_t rank = (_data[i / 4] >> ((i % 4) * 2)) & 0b11;
    switch (rank) {
    case 0:
      sequence.push_back('A');
      break;
    case 1:
      sequence.push_back('C');
      break;
    case 2:
      sequence.push_back('G');
      break;
    case 3:
      sequence.push_back('T');
      break;
    default:
      break;
    }
  }
  return sequence;
}

void Kmer::roll(uint8_t nt) {
  // proceed whole bytes
  uint8_t i = 0, _e;
  for (; i < _num_cells - 1; ++i) {
    _e = _data[i + 1] & 0b11;
    _data[i] = (_data[i] >> 2) | (_e << 6);
  }

  // finalise partial byte
  _data[i] = (_data[i] >> 2) | (nt << _edge_bit_offset);
}
