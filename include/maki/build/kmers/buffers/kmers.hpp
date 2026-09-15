
#pragma once
#include "base_buffer.hpp"
#include "maki/build/graph/build_colours.hpp"
#include "maki/build/utils/bits.hpp"
#include "maki/core/seq/io.hpp"
#include "nt_encoding.hpp"
#include <cstdint>
#include <sdsl/int_vector.hpp>

template <typename T, typename Derived>
class BaseKmerRandomIterator : public ndim::StridedIteratorBase<T, Derived> {
  using BaseIterType = ndim::StridedIteratorBase<T, Derived>;

public:
  using difference_type = std::ptrdiff_t;

protected:
  uint32_t _key_bytes, _value_offset;

public:
  BaseKmerRandomIterator(uint32_t key_bytes, uint32_t value_offset)
      : BaseIterType(), _key_bytes(key_bytes), _value_offset(value_offset) {}
  BaseKmerRandomIterator(T *p, size_t w, uint32_t key_bytes,
                         uint32_t value_offset)
      : BaseIterType(p, w), _key_bytes(key_bytes), _value_offset(value_offset) {
  }
  constexpr inline void writeKmer(Dna4SequenceConstIter in, uint8_t kEff) const
    requires(!std::is_const_v<T>)
  {
    writeMerTo(in, this->_data, kEff);
    return;
  }
  constexpr inline void writeKmer(const Kmer &kmer) const
    requires(!std::is_const_v<T>)
  {
    std::memcpy(this->_data, kmer.data(), _key_bytes);
    return;
  }
  inline uint8_t readEdge() const {
    return *(this->_data + _value_offset) & 0b111;
  }
  inline uint8_t keyMSB() const { return *(this->_data + _key_bytes - 1); }
};

template <typename T>
class KmerBufferRandomAccessIterator
    : public BaseKmerRandomIterator<T, KmerBufferRandomAccessIterator<T>> {
  using BaseIterType =
      BaseKmerRandomIterator<T, KmerBufferRandomAccessIterator<T>>;

protected:
  uint32_t _value_bytes;

public:
  KmerBufferRandomAccessIterator(uint32_t key_bytes, uint32_t value_bytes,
                                 uint32_t value_offset)
      : BaseIterType(key_bytes, value_offset), _value_bytes(value_bytes) {}
  KmerBufferRandomAccessIterator(T *p, size_t w, uint32_t key_bytes,
                                 uint32_t value_bytes, uint32_t value_offset)
      : BaseIterType(p, w, key_bytes, value_offset), _value_bytes(value_bytes) {
  }
  constexpr inline void writeValue(BufferValue val) const
    requires(!std::is_const_v<T>)
  {
    val.flush(this->_data + this->_value_offset, _value_bytes);
    return;
  }
  inline BufferValue readValue() const {
    return BufferValue(this->_data + this->_value_offset, _value_bytes);
  }
};

class KmerBuffer : public BaseKmerBuffer<KmerBufferRandomAccessIterator> {
  using BaseBufferType = BaseKmerBuffer<KmerBufferRandomAccessIterator>;

public:
  KmerBuffer(uint64_t num_records, uint32_t value_bytes, uint8_t _k,
             uint8_t _k_eff)
      : BaseBufferType(num_records, value_bytes, _k, _k_eff),
        _value_bytes(value_bytes) {}
  KmerBuffer(uint64_t num_records, uint32_t value_bytes, uint8_t k_eff)
      : BaseBufferType(num_records, value_bytes, k_eff, k_eff),
        _value_bytes(value_bytes) {}
  KmerBuffer(uint32_t value_bytes, uint8_t k, uint8_t k_eff,
             std::initializer_list<std::initializer_list<uint8_t>> data)
      : BaseBufferType(value_bytes, k, k_eff, data), _value_bytes(value_bytes) {
  }

protected:
  uint32_t _value_bytes;

public:
  inline iterator at(size_type idx) {
    return _data.at(idx, _key_bytes, _value_bytes, _value_offset);
  }
  inline const_iterator constAt(size_type idx) const {
    return _data.at(idx, _key_bytes, _value_bytes, _value_offset);
  }
  inline const_iterator at(size_type idx) const {
    return _data.at(idx, _key_bytes, _value_bytes, _value_offset);
  }
  inline uint32_t valueBytes() const noexcept { return _value_bytes; }
  inline iterator begin() {
    return _data.begin(_key_bytes, _value_bytes, _value_offset);
  }
  inline sentinel_type end() {
    return _num_records == _record_capacity
               ? _data.end(_key_bytes, _value_bytes, _value_offset)
               : _data.begin(_key_bytes, _value_bytes, _value_offset) +
                     _num_records;
  }

public:
  // --------------
  // Insert methods
  // --------------

  iterator insertRegion(iterator it, const random_dna4_range auto &sequence,
                        uint64_t colour, bool endIsTerminal) {
    auto seq_iter = std::ranges::begin(sequence);
    auto seq_end = std::ranges::end(sequence);

    Kmer kmer(k, seq_iter);
    seq_iter += k;

    uint8_t edge;

    while (seq_iter != seq_end) {
      edge = parsing::dna4ToShort(*seq_iter++);

      it.writeKmer(kmer);
      it.writeValue(BufferValue(colour, edge));
      ++it;

      kmer.roll(edge);
    }

    if (endIsTerminal) {
      it.writeKmer(kmer);
      it.writeValue(BufferValue((uint64_t)0u, terminalEdge));
      ++it;
    }

    return it;
  }

  iterator insertRegion(iterator it, const random_dna4_range auto &sequence,
                        uint64_t colour, bool endIsTerminal, ShortSuffix key) {
    auto seq_iter = std::ranges::begin(sequence);
    auto seq_end = std::ranges::end(sequence);

    uint8_t s = key.size();
    seq_iter += k - s;

    // initialise
    ShortSuffix suffix(s, seq_iter);
    seq_iter += s;

    uint64_t edge;

    // fill
    while (seq_iter != seq_end) {
      edge = parsing::dna4ToLong(*seq_iter);

      if (suffix == key) {
        it.writeKmer(seq_iter - k, k_eff);
        it.writeValue(BufferValue(colour, edge));
        ++it;
      }

      suffix.roll(edge);
      ++seq_iter;
    }

    // terminal edge
    if ((suffix == key) && endIsTerminal) {
      it.writeKmer(seq_iter - k, k_eff);
      it.writeValue(BufferValue((uint64_t)0u, terminalEdge));
      ++it;
    }

    return it;
  }

  template <random_dna4_range R>
  iterator insertFragment(iterator it, SequenceFragment<R> &&fmt) {
    return insertRegion(it, fmt.data(), fmt.id(), fmt.endIsTerminal());
  }

  template <random_dna4_range R>
  iterator insertFragment(iterator it, SequenceFragment<R> &&fmt,
                          ShortSuffix sfx_) {
    return insertRegion(it, fmt.data(), fmt.id(), fmt.endIsTerminal(), sfx_);
  }

  // --------------------
  // Container interfaces
  // --------------------

  iterator insertKmers(iterator it, const sequence_like auto &data) {
    return insertFragment(it, data.fragments(k));
  }

  iterator insertKmers(iterator it, const sequence_container_like auto &data) {
    for (auto &&fmt : data.fragments(k))
      it = insertFragment(it, std::forward<decltype(fmt)>(fmt));

    return it;
  }

  iterator insertKmers(iterator it,
                       const stranded_sequence_container_like auto &data) {
    for (auto &&fmt : data.forwardFragments(k))
      it = insertFragment(it, std::forward<decltype(fmt)>(fmt));

    for (auto &&fmt : data.reverseFragments(k))
      it = insertFragment(it, std::forward<decltype(fmt)>(fmt));

    return it;
  }

  iterator insertKmers(iterator it, const sequence_like auto &data,
                       ShortSuffix sfx_) {
    return insertFragment(it, data.fragments(k), sfx_);
  }

  iterator insertKmers(iterator it, const sequence_container_like auto &data,
                       ShortSuffix sfx_) {
    for (auto &&fmt : data.fragments(k))
      it = insertFragment(it, std::forward<decltype(fmt)>(fmt), sfx_);

    return it;
  }

  iterator insertKmers(iterator it,
                       const stranded_sequence_container_like auto &data,
                       ShortSuffix sfx_) {
    for (auto &&fmt : data.forwardFragments(k))
      it = insertFragment(it, std::forward<decltype(fmt)>(fmt), sfx_);

    for (auto &&fmt : data.reverseFragments(k))
      it = insertFragment(it, std::forward<decltype(fmt)>(fmt), sfx_);

    return it;
  }

  // ------------------------
  // Generic parallel wrapper
  // ------------------------
  void fillKmers(const container_span auto &data_,
                 std::span<const std::size_t> blocks_) {
    oneapi::tbb::parallel_for(
        (std::size_t)0, data_.size(), (std::size_t)1,
        [&](std::size_t i) { insertKmers(at(blocks_[i]), data_[i]); });
  }

  void fillKmers(const container_span auto &data_, std::span<const SuffixTable> blocks_,
                 ShortSuffix sfx_) {
    oneapi::tbb::parallel_for(
        (std::size_t)0, data_.size(), (std::size_t)1, [&](std::size_t i) {
          insertKmers(at(i == 0 ? 0 : blocks_[i - 1][sfx_]), data_[i], sfx_);
        });
  }

  void sort(KmerBuffer *temp);
  void sort();
};

struct KmerDiff {
  KmerDiff(uint32_t key_bytes) : _seqWidth(key_bytes) {}

  // `_lhs` and `_rhs` point to the least significant bytes of each k-mer
  KmerDiffClass operator()(const uint8_t *_lhs, const uint8_t *_rhs) const;

private:
  uint32_t _seqWidth;
};

struct MaskedBytesDiff {
  MaskedBytesDiff(uint8_t k_)
      : _width(key_size(k_)),
        _msb_mask(k_ % 4 == 0 ? 0xFFu : (0xFFu >> (2 * (4 - (k_ % 4))))) {}

protected:
  size_t _width;
  uint8_t _msb_mask;

public:
  /**
   * Compare bytes of two arrays, where `l_` and `r_` point to the least
   * significant bytes. A mask is applied to the most significant byte of `l_`.
   */
  KmerDiffClass operator()(const uint8_t *l_, const uint8_t *r_) const;
};

void adjacentDifference(KmerBuffer &, sdsl::int_vector<2> &);
sdsl::int_vector<2> adjacentDifference(KmerBuffer &);
