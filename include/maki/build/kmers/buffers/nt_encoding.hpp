
#pragma once
#include <cassert>
#include <cstdint>
#include <iterator>
#include <string>
#include <vector>

#include "base_buffer.hpp"
#include "maki/build/utils/bits.hpp"
#include "maki/core/seq/concepts.hpp"
#include "maki/core/seq/io.hpp"

/**
 * Represent nt using the 2-bit encoding: `[A, C, G, T] <=> [0, 1, 2, 3]`.
 * These are used for efficiently representing suffixes in a single cell
 * so that they can be quickly counted or compared. Most significant
 * base-pair of each window is the one most recently inserted.
 */
class ShortSuffix {
  friend class SuffixTable;

public:
  ShortSuffix(std::size_t _size, uint64_t init);
  ShortSuffix(std::size_t _size);

  ShortSuffix(std::size_t _size, random_dna4_iter auto _it)
      : ShortSuffix(_size, 0) {
    for (std::size_t i = 0; i < _s; ++i)
      _data |= parsing::dna4ToLong(*_it++) << (2 * i);
  }

  ShortSuffix(random_dna4_range auto &&rng_)
      : ShortSuffix(std::ranges::size(rng_), 0) {
    std::size_t i = 0;
    for (uint64_t nt : rng_ | std::views::transform([](seqan3::dna4 &&nt) {
                         return parsing::dna4ToLong(std::move(nt));
                       })) {
      _data |= nt << (2 * i++);
    }
  }

  static std::size_t numSuffixes(
      std::size_t s_); // number of suffixes up to and including size `s_`
  static ShortSuffix fromIndex(std::size_t r_, std::size_t s_);

  ShortSuffix(const ShortSuffix &other) = default;
  ShortSuffix(ShortSuffix &&other) = default;

  ShortSuffix &operator=(const ShortSuffix &other) = default;
  ShortSuffix &operator=(ShortSuffix &&other) = default;

  bool operator==(const ShortSuffix &other) const = default;
  bool operator!=(const ShortSuffix &other) const = default;

  inline uint64_t data() const { return _data; }
  inline std::size_t size() const { return _s; }
  std::string toString() const;

  // Return dna4 representation of most significant nucleotide, where
  // `TerminalEdge` represents the empty suffix
  inline uint8_t msb() const noexcept {
    return _s > 0 ? (_data >> (2 * (_s - 1))) : TerminalEdge;
  }

  /**
   * Rolls back base-pair to most significant position, removing least
   * significant nucleotide (does not increase size).
   */
  inline void roll(uint64_t nt) {
    _data = (_data >> 2) | (nt << (2 * (_s - 1)));
  }

  /**
   * In-place append a character to the least significant position. Increase
   * size by 1.
   */
  inline void push(uint64_t nt) {
    _data = (_data << 2) | (nt & 0b11);
    ++_s;
  }

  /**
   * Append a character to the least significant position. Increase size by 1.
   */
  inline ShortSuffix operator+(uint64_t nt) const {
    return ShortSuffix(_s + 1, (_data << 2) | (nt & 0b11));
  }

  /**
   * Shift in units of base-pairs. Equivalent to appending `w_` 'A's to the
   * least significant end.
   */
  inline ShortSuffix &operator<<=(std::size_t w_) {
    _s += w_;
    _data <<= 2 * w_;
    return *this;
  }

  /**
   * Shift in units of base-pairs. Equivalent to appending `w_` 'A's to the
   * least significant end.
   */
  inline ShortSuffix operator<<(std::size_t w_) const {
    return ShortSuffix(_s + w_, _data << (2 * w_));
  }

  /**
   * Shift in units of base-pairs. Removes `w_` least significant nucleotides.
   */
  inline ShortSuffix &operator>>=(std::size_t w_) {
    assert(w_ <= _s);
    _s -= w_;
    _data >>= 2 * w_;
    return *this;
  }

  /**
   * Shift in units of base-pairs. Removes `w_` least significant nucleotides.
   */
  inline ShortSuffix operator>>(std::size_t w_) const {
    return ShortSuffix(_s - w_, _data >> (2 * w_));
  }

  /**
   * Iterates from most-to-least significant byte of small sequence.
   */
  class Iterator {
  public:
    // initialise from value and width in bits
    Iterator(uint64_t &v_, uint8_t w_)
        : _value(v_), i(static_cast<int>((w_ - 1) / 8)) {
      assert(w_ > 0);
    }

    Iterator(const Iterator &rhs) : _value(rhs._value), i(rhs.i) {}
    Iterator(Iterator &&rhs) : _value(rhs._value), i(rhs.i) {}

  protected:
    uint64_t &_value;
    int i;

  public:
    using difference_type = int;
    using value_type = uint8_t;

    inline Iterator &operator=(const Iterator &rhs_) {
      _value = rhs_._value;
      i = rhs_.i;
      return *this;
    }

    inline Iterator &operator=(Iterator &&rhs_) {
      _value = rhs_._value;
      i = rhs_.i;
      return *this;
    }

    inline uint8_t operator*() const {
      assert(i >= 0);
      return static_cast<uint8_t>((_value >> (8 * i)) & 0xFF);
    }

    inline Iterator &operator++() {
      --i;
      return *this;
    }

    void operator++(int) { ++*this; }
  };

  // return iterator to most significant byte
  inline Iterator bytes() { return Iterator(_data, _s * 2); }

protected:
  std::size_t _s;
  uint64_t _data;
};

/**
 * Represent nt using the 2-bit encoding `[A, C, G, T] <=> [0, 1, 2, 3]`,
 * also allowing for a terminal edge `$ < A`. Bytes are ordered in vector
 * from least significant to most significant. This is useful for
 * efficiently encoding consecutive k-mers without needing to re-encode
 * the previous letters.
 */
class Kmer {
public:
  using Span = ndim::Span<const uint8_t *>;
  using Vector = std::vector<uint8_t>;
  using ReverseIterator = Vector::const_reverse_iterator;

public:
  Kmer(size_t num_nts)
      : _num_nts(num_nts), _num_cells(key_size(num_nts)),
        _edge_bit_offset(((num_nts - 1) % 4) * 2), _data(_num_cells, 0) {}

  Kmer(size_t _length, random_dna4_iter auto _it) : Kmer(_length) {
    writeMerTo(_it, _data.data(), _length);
  }

  Kmer(random_dna4_range auto &&rng_) : Kmer(std::ranges::distance(rng_)) {
    writeMerTo(std::ranges::begin(rng_), _data.data(), _num_nts);
  }

public:
  size_t _num_nts;
  uint8_t _num_cells, _edge_bit_offset;
  Vector _data;

public:
  inline size_t size() const noexcept { return _num_nts; }
  inline const uint8_t *data() const { return _data.data(); }
  inline ReverseIterator rbegin() const { return _data.rbegin(); }
  inline ReverseIterator rend() const { return _data.rend(); }
  std::string toString() const;
  void roll(uint8_t nt);
  inline Span view() const {
    return Span(_data.data(), static_cast<size_t>(_num_cells));
  }
};

/**
 * Subclass on contiguous rolling sequence that tracks how many
 * base-pairs have been inserted so far. This is used to represent
 * terminals at the beginning of contigs.
 */
struct LongSuffix : public Kmer {
  LongSuffix(size_t _length) : Kmer(_length), _size(0) {}

  /**
   * Adds nucleotide into most significant position, increasing size by 1.
   */
  inline void push(uint8_t nt) {
    Kmer::roll(nt);
    if (_size < _num_nts)
      ++_size;
  }

  size_t rank() const noexcept;

  // Returns number of base-pairs currently in terminal (starts with 0).
  inline size_t terminalLength() const noexcept { return _size; }

public:
  size_t _size;
};

class SuffixTable {
public:
  using iterator = std::vector<uint64_t>::iterator;
  using const_iterator = typename std::vector<uint64_t>::const_iterator;
  using pointer = uint64_t *;
  using const_pointer = const uint64_t *;

  inline void resize(size_t s_) {
    s = s_;
    _data.resize((size_t)1u << (2 * s), 0);
  }

  inline SuffixTable &operator=(const SuffixTable &other) = default;
  inline SuffixTable &operator=(SuffixTable &&other) = default;

  inline size_t suffixSize() const noexcept { return s; }

  const std::vector<uint64_t> &cdata() const { return _data; }
  inline iterator begin() { return _data.begin(); }
  inline iterator end() { return _data.end(); }
  inline const_iterator begin() const { return _data.begin(); }
  inline const_iterator end() const { return _data.end(); }
  inline const_iterator cbegin() const { return _data.cbegin(); }
  inline const_iterator cend() const { return _data.cend(); }

  inline uint64_t &operator[](const ShortSuffix &sfx) {
    return _data.operator[](sfx._data);
  }
  inline const uint64_t &operator[](uint64_t sfx) const {
    return _data.operator[](sfx);
  }
  inline const uint64_t &operator[](const ShortSuffix &sfx) const {
    return _data.operator[](sfx._data);
  }

  size_t maxValue() const;
  inline size_t size() const noexcept { return _data.size(); }

  friend SuffixTable &operator+=(SuffixTable &lhs, const SuffixTable &rhs);

public:
  // --------------
  // Insert methods
  // --------------

  void countRegion(const random_dna4_range auto &sequence) {
    ShortSuffix suffix{s};

    auto it = std::ranges::cbegin(sequence);
    auto end = std::ranges::cend(sequence);

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

  template <random_dna4_range R>
  void countKmerSuffixes(SequenceFragment<R> &&fmt, std::size_t k,
                         std::size_t offset) {
    auto &sequence = fmt.data();
    std::size_t size = fmt.size();

    if (fmt.endIsTerminal())
      sequence = sequence | std::views::take(--size);

    if (size >= k)
      countRegion(sequence | std::views::drop(offset));
  }

  // --------------------
  // Container interfaces
  // --------------------

  void count(const sequence_like auto &data, std::size_t k,
             std::size_t offset) {
    countKmerSuffixes(data.fragments(k), k, offset);
  }

  void count(const sequence_container_like auto &data, std::size_t k,
             std::size_t offset) {
    std::ranges::for_each(data.fragments(k), [&](auto &&fmt) {
      countKmerSuffixes(std::forward<decltype(fmt)>(fmt), k, offset);
    });
  }

  void count(const stranded_sequence_container_like auto &data, std::size_t k, std::size_t offset) {
    std::ranges::for_each(data.forwardFragments(k), [&](auto &&fmt) {
      countKmerSuffixes(std::forward<decltype(fmt)>(fmt), k, offset);
    });

    std::ranges::for_each(data.reverseFragments(k), [&](auto &&fmt) {
      countKmerSuffixes(std::forward<decltype(fmt)>(fmt), k, offset);
    });
  }

private:
  size_t s;
  std::vector<uint64_t> _data;
};

class LongSuffixGate {
public:
  LongSuffixGate(std::size_t _length);
  std::optional<std::size_t> trySet(const LongSuffix &sfx, uint8_t dna4_edge);
  std::size_t count() const;

private:
  std::size_t _words;
  std::unique_ptr<std::atomic_uint8_t[]> _locks;
  std::atomic_size_t _elems;

  std::size_t calculateBitPosition(const LongSuffix &sfx, uint8_t dna4_edge);
  bool trySetBit(const LongSuffix &sfx, uint8_t dna4_edge);
};

/**
 * Count suffixes in each genome, and accumulate. The result is a 2D matrix
 * where `rows ~ genome` and `columns ~ suffix`, and each cell contains the
 * cumulative number of k-mers with that suffix including all genomes before it
 * in the container.
 */
template <sequence_fragment_container T>
std::vector<SuffixTable> createSuffixPlan(std::span<const T> data_, std::size_t k,
                                          std::size_t s, std::size_t offset = 0,
                                          bool accumulate = true) {
  std::vector<SuffixTable> tables(data_.size());

  // count suffixes
  oneapi::tbb::parallel_for((std::size_t)0, data_.size(), (std::size_t)1,
                            [&, k, s](std::size_t i) {
                              tables[i].resize(s);
                              tables[i].count(data_[i], k, offset);
                            });

  if (accumulate) {
    for (std::size_t i = 1; i < tables.size(); ++i) {
      tables[i] += tables[i - 1];
    }
  }

  return tables;
}
