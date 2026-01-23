
#pragma once
#include "base_buffer.hpp"
#include "maki/build/graph/build_colours.hpp"
#include "maki/build/utils/bits.hpp"
#include "maki/core/seq/seq_io.hpp"
#include "nt_encoding.hpp"
#include <cstdint>

template <typename T, typename Derived>
class BaseKmerRandomIterator : public ndim::StridedIteratorBase<T, Derived> {
  using BaseIterType = ndim::StridedIteratorBase<T, Derived>;
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
  /**
   * Fill with k-mer sequences that have a specific suffix.
   */
  void fill(const std::vector<const SequenceContainer *> &data_,
            const std::vector<SuffixTable> &blocks_, ShortSuffix sfx_);
  iterator insert(iterator it, Dna4SequenceConstIter begin,
                  Dna4SequenceConstIter end, uint64_t colour,
                  bool endIsTerminal);
  iterator insert(iterator it, Dna4SequenceConstIter begin,
                  Dna4SequenceConstIter end, uint64_t colour,
                  bool endIsTerminal, ShortSuffix key);
  void sort(KmerBuffer *temp);
  void sort();
};

enum KmerDiffClass : uint8_t {
  IS_0 = 0b00u,  // corresponds to pattern `+  +  + ... +`
  IS_K = 0b01u,  // corresponds to pattern `-  +  + ... +`
  BW_0_K = 0b10u // corresponds to pattern `* ... - ... *`
};

struct KmerDiff {
  static constexpr uint8_t result_width = 2;

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
