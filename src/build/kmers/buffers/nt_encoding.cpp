#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include <iterator>
#include <oneapi/tbb/parallel_for.h>

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

std::size_t ShortSuffix::numSuffixes(std::size_t s_) {
  /*
    1 + 4 + 4^2 + ... + 4^s = 4 * (4^(s + 1) - 1) / 3 (geometric series)
  */
  return (((std::size_t)1 << (2 * (s_ + 1))) - 1u) / 3u;
}

ShortSuffix ShortSuffix::fromRank(std::size_t r_) {
  std::size_t l = 0u;
  while (numSuffixes(l) < r_)
    ++l;
  if (l == 0u)
    return ShortSuffix(0u);
  else
    return ShortSuffix(l, r_ - numSuffixes(l - 1u) - 1u);
}

std::string ShortSuffix::toString() const {
  if (_s == 0)
    return "NULL";
  std::string nucleotideSequence;
  nucleotideSequence.reserve(_s);
  for (size_t i = 0; i < _s; ++i) {
    uint64_t rank = (_data >> (2 * i)) & 0b11;
    switch (rank) {
    case 0:
      nucleotideSequence.push_back('A');
      break;
    case 1:
      nucleotideSequence.push_back('C');
      break;
    case 2:
      nucleotideSequence.push_back('G');
      break;
    case 3:
      nucleotideSequence.push_back('T');
      break;
    default:
      break;
    }
  }
  return nucleotideSequence;
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

void SuffixTable::count(Dna4SequenceConstIter it, Dna4SequenceConstIter end) {
  ShortSuffix suffix{s};

  // initialise suffix
  for (size_t i = 0; i < s; ++i)
    suffix.roll(parsing::dna4ToLong(*it++));
  ++(*this)[suffix];

  // count rest of sequence
  while (it != end) {
    suffix.roll(parsing::dna4ToLong(*it++));
    ++(*this)[suffix];
  }
}

size_t SuffixTable::maxValue() const {
  return *std::max_element(_data.cbegin(), _data.cend());
}

SuffixTable &operator+=(SuffixTable &lhs, const SuffixTable &rhs) {
  assert(lhs.suffixSize() == rhs.suffixSize());
  auto it = lhs.begin();
  for (uint64_t v : rhs) {
    *it += v;
    ++it;
  }
  return lhs;
}

std::vector<SuffixTable>
createSuffixPlan(const std::vector<const SequenceContainer *> &data,
                 std::size_t k, std::size_t s, bool accumulate) {
  std::vector<SuffixTable> tables(data.size());

  // count suffixes
  oneapi::tbb::parallel_for((std::size_t)0, data.size(), (std::size_t)1,
                            [&](std::size_t i) {
                              tables[i].resize(s);
                              for (auto fmt : data[i]->fragments(k)) {
                                auto it = fmt.begin(), end = fmt.end();
                                if (it != end && !fmt.endIsTerminal())
                                  --end;
                                std::size_t size = end - it;
                                if (size >= k)
                                  tables[i].count(it + (k - s), end);
                              }
                            });

  if (accumulate) {
    for (std::size_t i = 1; i < tables.size(); ++i) {
      tables[i] += tables[i - 1];
    }
  }

  return tables;
}
