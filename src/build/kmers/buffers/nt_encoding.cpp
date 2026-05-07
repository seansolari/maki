#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/core/seq/io.hpp"
#include <iterator>
#include <oneapi/tbb/parallel_for.h>
#include <optional>

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
    1 + 4 + 4^2 + ... + 4^s = 1 * (4^(s + 1) - 1) / 3 (geometric series)
  */
  return (((std::size_t)1 << (2 * (s_ + 1))) - 1u) / 3u;
}

ShortSuffix ShortSuffix::fromIndex(std::size_t r_, std::size_t s_) {
  std::size_t x = 0u, l = s_, ls, f;
  while (r_ > 0u) {
    assert(l > 0u);
    ls = numSuffixes(l - 1u), f = (r_ - 1u) / ls;
    x = (x << 2u) | f;
    r_ -= f * ls + 1u;
    l -= 1;
  }
  return ShortSuffix(s_ - l, x);
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

size_t LongSuffix::rank() const noexcept {
  size_t res = 0;
  for (size_t i = 0; 4 * i < _size; ++i) {
    res |= static_cast<size_t>(_data[i]) << (size_t(8) * i);
  }
  return res;
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
                 std::size_t k, std::size_t s, std::size_t offset,
                 bool accumulate) {
  std::vector<SuffixTable> tables(data.size());

  // count suffixes
  oneapi::tbb::parallel_for((std::size_t)0, data.size(), (std::size_t)1,
                            [&, k, s](std::size_t i) {
                              tables[i].resize(s);
                              for (auto fmt : data[i]->fragments(k)) {
                                auto it = fmt.begin(), end = fmt.end();
                                if (it != end && !fmt.endIsTerminal())
                                  --end;
                                std::size_t size = end - it;
                                if (size >= k)
                                  tables[i].count(it + offset, end);
                              }
                            });

  if (accumulate) {
    for (std::size_t i = 1; i < tables.size(); ++i) {
      tables[i] += tables[i - 1];
    }
  }

  return tables;
}

LongSuffixGate::LongSuffixGate(std::size_t _length)
    : _words((ShortSuffix::numSuffixes(_length) * 5 + 7) / 8),
      _locks(std::make_unique<std::atomic_uint8_t[]>(_words)), _elems(0) {}

std::optional<std::size_t> LongSuffixGate::trySet(const LongSuffix &sfx,
                                                  uint8_t dna4_edge) {
  if (trySetBit(sfx, dna4_edge)) {
    return std::make_optional(_elems++);
  } else
    return std::nullopt;
}

std::size_t LongSuffixGate::count() const {
  return _elems.load(std::memory_order_relaxed);
}

bool LongSuffixGate::trySetBit(const LongSuffix &sfx, uint8_t dna4_edge) {
  std::size_t suffix_index = calculateBitPosition(sfx, dna4_edge),
              w = suffix_index / 8, b = suffix_index % 8;
  assert(w < _words);
  uint8_t mask = 1U << b, expected = _locks[w].load(std::memory_order_relaxed);
  while (!(expected & mask)) {
    if (_locks[w].compare_exchange_weak(expected, expected | mask,
                                        std::memory_order_acquire,
                                        std::memory_order_relaxed)) {
      return true;
    }
  }
  return false;
}

std::size_t LongSuffixGate::calculateBitPosition(const LongSuffix &sfx,
                                                 uint8_t dna4_edge) {
  auto dna5_edge = parsing::dna4ToDna5(dna4_edge);
  if (sfx.size() == 0)
    return dna5_edge;
  else {
    return 5 * (ShortSuffix::numSuffixes(sfx.size() - 1) + sfx.rank()) +
           dna5_edge;
  }
}
