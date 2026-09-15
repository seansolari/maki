
#pragma once
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ranges>
#include <sdsl/int_vector.hpp>
#include <type_traits>

#include <oneapi/tbb/parallel_for.h>

#include "base_buffer.hpp"
#include "kmers.hpp"
#include "maki/build/utils/bits.hpp"
#include "maki/core/seq/concepts.hpp"
#include "nt_encoding.hpp"

template <typename T>
class TerminalBufferRandomAccessIterator
    : public ndim::StridedIteratorBase<T,
                                       TerminalBufferRandomAccessIterator<T>> {
  using BaseIterType =
      ndim::StridedIteratorBase<T, TerminalBufferRandomAccessIterator<T>>;

public:
  TerminalBufferRandomAccessIterator() = default;
  TerminalBufferRandomAccessIterator(uint32_t lengthBytes, uint32_t seqBytes,
                                     uint32_t valueOffset, uint8_t kEff)
      : BaseIterType(), _lengthBytes(lengthBytes), _seqBytes(seqBytes),
        _valueOffset(valueOffset), _k_eff(kEff) {}

  TerminalBufferRandomAccessIterator(T *p, size_t w, uint32_t lengthBytes,
                                     uint32_t seqBytes, uint32_t valueOffset,
                                     uint8_t kEff)
      : BaseIterType(p, w), _lengthBytes(lengthBytes), _seqBytes(seqBytes),
        _valueOffset(valueOffset), _k_eff(kEff) {}

  TerminalBufferRandomAccessIterator(
      const TerminalBufferRandomAccessIterator &) = default;
  TerminalBufferRandomAccessIterator &
  operator=(const TerminalBufferRandomAccessIterator &) = default;

  using difference_type = std::ptrdiff_t;

  constexpr inline void writeTerminal(const LongSuffix &terminal) const
    requires(!std::is_const_v<T>)
  {
    std::memcpy(this->_data + _lengthBytes, terminal.data(), _seqBytes);
  }
  constexpr inline void writeSize(size_t v_) const
    requires(!std::is_const_v<T>)
  {
    ser(v_, this->_data, _lengthBytes);
  }
  constexpr inline void writeKey(LongSuffix const &tl) const
    requires(!std::is_const_v<T>)
  {
    writeSize(tl.terminalLength());
    writeTerminal(tl);
  }
  constexpr inline void writeKey(random_dna4_iter auto &&_it,
                                 std::size_t size_) const
    requires(!std::is_const_v<T>)
  {
    writeSize(size_);
    reverseWriteMerTo(_it, this->_data + _lengthBytes, _k_eff, size_);
  }
  constexpr inline void writeEdge(T edge_) const
    requires(!std::is_const_v<T>)
  {
    *(this->_data + _valueOffset) = edge_;
  }
  inline uint32_t lengthBytes() const noexcept { return _lengthBytes; }
  inline uint32_t keyBytes() const noexcept { return _seqBytes + _lengthBytes; }
  inline size_t readSize() const { return deser(this->_data, _lengthBytes); }
  inline T readEdge() const { return *(this->_data + _valueOffset); }
  inline T readValue() const { return readEdge(); }
  inline T *data() { return this->_data; }
  inline const T *cdata() const { return this->_data; }
  inline T keyMSB() const { return *(this->_data + keyBytes() - 1); }

protected:
  uint32_t _lengthBytes, _seqBytes, _valueOffset;
  uint8_t _k_eff;
};

class TerminalRange;

/**
 * Values are serialised in orientation (LSB --> MSB):
 *
 *      | --- length --- | --- sequence --- | 0 ... 0 | edge |
 */
class TerminalBuffer
    : public BaseKmerBuffer<TerminalBufferRandomAccessIterator> {
  friend struct TerminalsLessThan;

protected:
  using BaseBufferType = BaseKmerBuffer<TerminalBufferRandomAccessIterator>;

public:
  static constexpr struct autofit_tag_t {
  } autofit_tag{};

  explicit TerminalBuffer(uint64_t numRecords, uint8_t k_, uint8_t k_eff_,
                          uint32_t lengthBytes)
      : BaseBufferType(numRecords, _valueBytes, k_, k_eff_, lengthBytes),
        _lengthBytes(lengthBytes), _seqBytes(_key_bytes - _lengthBytes) {}

  explicit TerminalBuffer(uint64_t numRecords, uint8_t k_, uint8_t k_eff_,
                          autofit_tag_t)
      : TerminalBuffer(numRecords, k_, k_eff_, required_bytes(k_eff_)) {}

  explicit TerminalBuffer(uint64_t numRecords, uint8_t k_, uint32_t lengthBytes)
      : BaseBufferType(numRecords, _valueBytes, k_, k_, lengthBytes),
        _lengthBytes(lengthBytes), _seqBytes(_key_bytes - _lengthBytes) {}

  explicit TerminalBuffer(uint64_t numRecords, uint8_t k_, autofit_tag_t)
      : TerminalBuffer(numRecords, k_, required_bytes(k_)) {}

  explicit TerminalBuffer(
      uint8_t k_, uint8_t k_eff_, uint32_t lengthBytes,
      std::initializer_list<std::initializer_list<uint8_t>> data)
      : BaseBufferType(_valueBytes, k_, k_eff_, lengthBytes, data),
        _lengthBytes(lengthBytes), _seqBytes(_key_bytes - _lengthBytes) {}

  explicit TerminalBuffer(
      uint8_t k_, uint8_t k_eff_,
      std::initializer_list<std::initializer_list<uint8_t>> data, autofit_tag_t)
      : TerminalBuffer(k_, k_eff_, required_bytes(k_eff_), data) {}

  explicit TerminalBuffer(
      uint8_t k_, uint32_t lengthBytes,
      std::initializer_list<std::initializer_list<uint8_t>> data)
      : BaseBufferType(_valueBytes, k_, k_, lengthBytes, data),
        _lengthBytes(lengthBytes), _seqBytes(_key_bytes - _lengthBytes) {}

  explicit TerminalBuffer(
      uint8_t k_, std::initializer_list<std::initializer_list<uint8_t>> data,
      autofit_tag_t)
      : TerminalBuffer(k_, required_bytes(k_), data) {}

  template <std::input_iterator InputIt>
  TerminalBuffer(const TerminalBuffer &src_, InputIt it_, InputIt end_)
      : TerminalBuffer(std::distance(it_, end_), src_.getK(), src_.getEffK(),
                       src_.lengthBytes()) {
    tbb::parallel_for((size_t)0, (size_t)_num_records, [&](size_t i) {
      auto __srcspan = *src_.at(*(it_ + i));
      std::copy(__srcspan.begin(), __srcspan.end(), data(i));
    });
  }

private:
  uint32_t _lengthBytes, _seqBytes;
  static constexpr uint32_t _valueBytes = 1u;

public:
  inline uint32_t lengthBytes() const noexcept { return _lengthBytes; }
  inline static constexpr uint32_t valueBytes() noexcept { return _valueBytes; }

  /**
   * Will not reallocate - can only reduce size. Does not modify
   * values (even those now out of reach).
   */
  void shrink(size_t new_size_);

public:
  inline iterator begin() {
    return _data.begin(_lengthBytes, _seqBytes, valueOffset(), k_eff);
  }
  inline iterator end() {
    return _num_records == _record_capacity
               ? _data.end(_lengthBytes, _seqBytes, valueOffset(), k_eff)
               : _data.begin(_lengthBytes, _seqBytes, valueOffset(), k_eff) +
                     _num_records;
  }
  inline const_iterator constBegin() const noexcept {
    return _data.cbegin(_lengthBytes, _seqBytes, valueOffset(), k_eff);
  }
  inline const_iterator constEnd() const noexcept {
    return _num_records == _record_capacity
               ? _data.cend(_lengthBytes, _seqBytes, valueOffset(), k_eff)
               : _data.cbegin(_lengthBytes, _seqBytes, valueOffset(), k_eff) +
                     _num_records;
  }
  inline const_iterator begin() const noexcept { return constBegin(); }
  inline const_iterator end() const noexcept { return constEnd(); }
  inline iterator at(size_type idx) {
    return _data.at(idx, _lengthBytes, _seqBytes, valueOffset(), k_eff);
  }
  inline const_iterator constAt(size_type idx) const {
    return _data.at(idx, _lengthBytes, _seqBytes, valueOffset(), k_eff);
  }
  inline const_iterator at(size_type idx) const {
    return _data.at(idx, _lengthBytes, _seqBytes, valueOffset(), k_eff);
  }

public:
  // --------------
  // Insert methods
  // --------------

  // insert a sequence range
  iterator insertRegion(iterator it, const random_dna4_range auto &sequence,
                        bool endIsTerminal) {
    LongSuffix tl(k_eff /* here, will be same as k */);

    for (uint8_t edge : sequence | std::views::transform([](seqan3::dna4 &&nt) {
                          return parsing::dna4ToShort(std::move(nt));
                        })) {
      it.writeKey(tl);
      it.writeEdge(edge);
      ++it;

      tl.push(edge);
    }

    if (endIsTerminal) {
      it.writeKey(tl);
      it.writeEdge(terminalEdge);
      ++it;
    }

    return it;
  }

  // insert terminals from a sequence range with a given suffix
  iterator insertRegion(iterator it, const random_dna4_range auto &sequence,
                        bool endIsTerminal, ShortSuffix key) {
    auto seq_iter = std::ranges::begin(sequence);
    auto seq_end = std::ranges::end(sequence);

    uint8_t s = key.size();
    size_t distanceFromStart = 0, effTerminalSize;

    // initialise

    ShortSuffix suffix(s, seq_iter);
    seq_iter += s;

    // fill

    uint64_t edge;

    while (seq_iter != seq_end) {
      edge = parsing::dna4ToLong(*seq_iter);

      if (suffix == key) {
        effTerminalSize = std::min(distanceFromStart, (size_t)k_eff);
        it.writeKey(seq_iter - (s + 1), effTerminalSize);
        it.writeEdge(edge);
        ++it;
      }

      suffix.roll(edge);
      ++seq_iter;
      ++distanceFromStart;
    }

    // terminal edge

    if ((suffix == key) && endIsTerminal) {
      effTerminalSize = std::min(distanceFromStart, (size_t)k_eff);
      it.writeKey(seq_iter - (s + 1), effTerminalSize);
      it.writeEdge(terminalEdge);
      ++it;
    }

    return it;
  }

  // competitively insert a suffix if it has not been seen before (`lock` must
  // be shared between threads)
  void insertRegion(LongSuffixGate &lock, const random_dna4_range auto &sequence,
                    bool endIsTerminal) {
    LongSuffix tl(k_eff /* here, will be same as k */);

    std::optional<std::size_t> p = std::nullopt;

    for (uint8_t edge : sequence | std::views::transform([](seqan3::dna4 &&nt) {
                          return parsing::dna4ToShort(std::move(nt));
                        })) {
      if (p = lock.trySet(tl, edge); p) {
        auto it = at(p.value());
        it.writeKey(tl);
        it.writeEdge(edge);
      }

      tl.push(edge);
    }

    if (endIsTerminal) {
      if (p = lock.trySet(tl, terminalEdge); p) {
        auto it = at(p.value());
        it.writeKey(tl);
        it.writeEdge(terminalEdge);
      }
    }
  }

  template <random_dna4_range R>
  iterator insertFragmentTerminals(iterator it, SequenceFragment<R> &&fmt) {
    return insertRegion(it, fmt.data() | std::views::take(k), false);
  }

  template <random_dna4_range R>
  void insertFragmentTerminals(LongSuffixGate &lock, SequenceFragment<R> &&fmt) {
    insertRegion(lock, fmt.data() | std::views::take(k), false);
  }

  template <random_dna4_range R>
  iterator insertFragmentKmers(iterator it, SequenceFragment<R> &&fmt, ShortSuffix sfx_) {
    return insertRegion(it, fmt.data(), fmt.endIsTerminal(), sfx_);
  }

  // --------------------
  // Container interfaces
  // --------------------

  iterator insertTerminals(iterator it, const sequence_like auto &data) {
    return insertFragmentTerminals(it, data.terminals());
  }

  iterator insertTerminals(iterator it, const sequence_container_like auto &data) {
    for (auto &&fmt : data.terminals())
      it = insertFragmentTerminals(it, std::forward<decltype(fmt)>(fmt));

    return it;
  }

  iterator insertTerminals(iterator it, const stranded_sequence_container_like auto &data) {
    for (auto &&fmt : data.forwardTerminals())
      it = insertFragmentTerminals(it, std::forward<decltype(fmt)>(fmt));
    
    for (auto &&fmt : data.reverseTerminals())
      it = insertFragmentTerminals(it, std::forward<decltype(fmt)>(fmt));

    return it;
  }

  void insertTerminals(LongSuffixGate &lock, const sequence_like auto &data) {
    insertFragmentTerminals(lock, data.terminals());
  }

  void insertTerminals(LongSuffixGate &lock, const sequence_container_like auto &data) {
    for (auto &&fmt : data.terminals())
      insertFragmentTerminals(lock, std::forward<decltype(fmt)>(fmt));
  }

  void insertTerminals(LongSuffixGate &lock, const stranded_sequence_container_like auto &data) {
    for (auto &&fmt : data.forwardTerminals())
      insertFragmentTerminals(lock, std::forward<decltype(fmt)>(fmt));

    for (auto &&fmt : data.reverseTerminals())
      insertFragmentTerminals(lock, std::forward<decltype(fmt)>(fmt));
  }

  iterator insertKmers(iterator it, const sequence_like auto &data, ShortSuffix sfx_) {
    return insertFragmentKmers(it, data.fragments(k), sfx_);
  }

  iterator insertKmers(iterator it, const sequence_container_like auto &data, ShortSuffix sfx_) {
    for (auto &&fmt : data.fragments(k))
      it = insertFragmentKmers(it, std::forward<decltype(fmt)>(fmt), sfx_);

    return it;
  }

  iterator insertKmers(iterator it, const stranded_sequence_container_like auto &data, ShortSuffix sfx_) {
    for (auto &&fmt : data.forwardFragments(k))
      it = insertFragmentKmers(it, std::forward<decltype(fmt)>(fmt), sfx_);
    
    for (auto &&fmt : data.reverseFragments(k))
      it = insertFragmentKmers(it, std::forward<decltype(fmt)>(fmt), sfx_);

    return it;
  }

  // ------------------------
  // Generic parallel wrapper
  // ------------------------

  template <sequence_fragment_container T>
  void fillTerminals(std::span<const T> data_, std::span<const std::size_t> blocks_) {
    oneapi::tbb::parallel_for(
        (std::size_t)0, data_.size(), (std::size_t)1, [&](std::size_t i) {
          insertTerminals(at(i == 0 ? 0 : blocks_[i - 1]), data_[i]);
        });
  }

  template <sequence_fragment_container T>
  void fillTerminals(std::span<const T> data_, LongSuffixGate &lock) {
    oneapi::tbb::parallel_for(
        (std::size_t)0, data_.size(), (std::size_t)1, [&](std::size_t i) {
          insertTerminals(lock, data_[i]);
        });
  }

  template <sequence_fragment_container T>
  void fillKmers(std::span<const T> data_, std::span<const std::size_t> blocks_, ShortSuffix sfx_) {
    oneapi::tbb::parallel_for(
        (std::size_t)0, data_.size(), (std::size_t)1, [&](std::size_t i) {
          insertKmers(at(i == 0 ? 0 : blocks_[i - 1]), data_[i], sfx_);
        });
  }

  // -------------------
  // Downstream analysis
  // -------------------

  void sort(TerminalBuffer *temp);
  void sort();
  void unique();
  TerminalBuffer OOPsort(bool makeUnique = true) const;
  TerminalRange asRange() const;
};

using TerminalConstIter = typename TerminalBuffer::const_iterator;

class TerminalRange {
public:
  using iterator = typename TerminalBuffer::const_iterator;
  using const_iterator = typename TerminalBuffer::const_iterator;

public:
  TerminalRange() = default;
  TerminalRange(uint8_t k_eff_, long lengthBytes_, TerminalConstIter begin_,
                TerminalConstIter end_)
      : _k_eff(k_eff_), _lengthBytes(lengthBytes_), _begin(begin_), _end(end_) {
  }
  TerminalRange(TerminalRange &rhs, TerminalConstIter pivot)
      : _k_eff(rhs._k_eff), _lengthBytes(rhs._lengthBytes), _begin(rhs._begin),
        _end(pivot) {
    rhs._begin = pivot;
  }

protected:
  uint8_t _k_eff;
  long _lengthBytes;
  TerminalConstIter _begin, _end;

public:
  inline TerminalConstIter constBegin() const noexcept { return _begin; }
  inline TerminalConstIter constEnd() const noexcept { return _end; }
  inline TerminalConstIter begin() const noexcept { return constBegin(); }
  inline TerminalConstIter end() const noexcept { return constEnd(); }
  inline uint8_t getEffK() const noexcept { return _k_eff; }
  inline uint32_t lengthBytes() const noexcept { return _lengthBytes; }
  inline size_t size() const noexcept { return _end - _begin; }
  inline bool empty() const noexcept { return _begin == _end; }
  inline void stepForward(size_t di_) { _begin += di_; }

public:
  /**
   * Find lower bound `target` in the range `[begin, end)`.
   */
  TerminalConstIter lowerBound(ndim::Span<uint8_t *>, uint8_t) const;
  TerminalConstIter lowerBound(ndim::Span<uint8_t *>) const;

  /**
   * Find upper bound `target` in the range `[begin, end)`.
   */
  TerminalConstIter upperBound(ndim::Span<uint8_t *>, uint8_t) const;
  TerminalConstIter upperBound(ndim::Span<uint8_t *>) const;

public:
  TerminalRange endsWith(ShortSuffix) const;
  TerminalRange retrieve(ShortSuffix) const;
};

// Operators

struct TerminalsLessThan {
public:
  TerminalsLessThan(const TerminalBuffer &buffer_)
      : _buffer(buffer_), _keyWidth(buffer_.keyBytes()),
        _valueOffset(buffer_.valueOffset()) {}

  bool operator()(size_t i_, size_t j_) const;

protected:
  static int _rng_cmp(const uint8_t *l_, const uint8_t *r_, uint32_t width);

  TerminalBuffer const &_buffer;
  uint32_t _keyWidth, _valueOffset;
};

struct TerminalIsSuffix {
  TerminalIsSuffix(uint8_t k_, uint8_t s_, long offset_ = 0)
      : _trg_width(key_size(k_) + offset_), _sfx_shift((k_ - s_) % 4),
        _sfx_width((s_ + _sfx_shift + 3) / 4), _offset(offset_) {
    assert(s_ <= k_);
  }

  bool operator()(uint8_t const *trg_, uint64_t sfx_) const;
  bool operator()(TerminalConstIter it, ShortSuffix const &sfx) const;

protected:
  long _trg_width;   // number of bytes per search target (k-mer)
  size_t _sfx_shift; // shift each input suffix to align with k-mer byte
                     // orientation
  long _sfx_width;   // number of bytes per (shifted) search query (terminal)
  long _offset;
};

/**
 * Compare two sequences whose bytes first need to be aligned (i.e. compare
 * suffix with a terminal)
 */
template <bool mask_lsb> struct UnalignedBytesLessThan {
  /**
   * `k_` defines the template orientation, and we need to modify `s_`
   * to fit the template.
   */
  UnalignedBytesLessThan(uint8_t k_, uint8_t s_, uint8_t offset_ = 0)
      : _trg_width(key_size(k_) + offset_), _sfx_shift((k_ - s_) % 4),
        _sfx_width((s_ + _sfx_shift + 3) / 4),
        _trg_lsb_msk((0xFFu << (2 * _sfx_shift)) & 0xFFu) {
    assert(s_ <= k_);
  }

  template <typename Op>
  bool operator()(ndim::Span<uint8_t const *> trg_, uint64_t sfx_,
                  Op cmp) const {
    uint8_t const *trg_aln = trg_.data() + _trg_width - _sfx_width;
    sfx_ <<= 2 * _sfx_shift;

    for (long i = _sfx_width - 1; i > -1; --i) {
      uint8_t trg_v = *(trg_aln + i), sfx_v = (sfx_ >> (8 * i)) & 0xFFu;

      if constexpr (mask_lsb) {
        if (i == 0)
          trg_v &= _trg_lsb_msk;
      }

      if (cmp(trg_v, sfx_v))
        return true;
      else if (cmp(sfx_v, trg_v))
        return false;
    }

    return false;
  }
  inline bool operator()(ndim::Span<const uint8_t *> trg_,
                         uint64_t sfx_) const {
    return operator()(trg_, sfx_, std::less{});
  }
  inline bool operator()(uint64_t sfx_,
                         ndim::Span<const uint8_t *> trg_) const {
    return operator()(trg_, sfx_, std::greater{});
  }

protected:
  long _trg_width;      // number of bytes per search target (k-mer)
  size_t _sfx_shift;    // shift each input suffix to align with k-mer byte
                        // orientation
  long _sfx_width;      // number of bytes per (shifted) search query (terminal)
  uint8_t _trg_lsb_msk; // apply mask to least-significant-byte of search target
};

enum mask_parity { LHS, RHS };

/**
 * Compare two sequences whose bytes are aligned, according to a given k-mer
 * size. Used for comparing k-mers and terminals, where the terminal size may be
 * greater than the k-mer size (and so its most significant bytes should be
 * masked). The constructor parameter `k_` should be the smaller.
 */
template <mask_parity P> struct MaskedAlignedBytesLessThan {
  MaskedAlignedBytesLessThan(uint8_t k_, long offset_ = 0)
      : _width(key_size(k_)),
        _msb_mask(k_ % 4 == 0 ? 0xFFu : (0xFFu >> (2 * (4 - (k_ % 4))))),
        _offset(offset_) {}

  /**
   * Compare bytes of two arrays, where `l_` and `r_` point to the least
   * significant bytes. A mask is applied to the most significant byte, the
   * target of which is defined by the template parameter `P`.
   */
  bool operator()(const uint8_t *l_, const uint8_t *r_) const {
    uint8_t msk = _msb_mask;

    for (const uint8_t *lhs = l_ + _width - 1, *rhs = r_ + _width - 1;
         lhs != l_ - 1; (void)--lhs, (void)--rhs) {
      uint8_t l_val = *lhs, r_val = *rhs;

      if constexpr (P == LHS)
        l_val &= msk;
      else
        r_val &= msk;

      if (l_val < r_val)
        return true;
      else if (l_val > r_val)
        return false;

      msk = 0xFFu;
    }
    return false;
  }

  template <typename X, typename Y>
  inline bool operator()(ndim::Span<X *> lhs_, ndim::Span<Y *> rhs_) const {
    if constexpr (P == LHS)
      return operator()(lhs_.data() + _offset, rhs_.data());
    else
      return operator()(lhs_.data(), rhs_.data() + _offset);
  }

protected:
  size_t _width;
  uint8_t _msb_mask;
  long _offset;
};

struct TerminalDiff : public KmerDiff {
  TerminalDiff(uint32_t lengthBytes, uint8_t k_eff)
      : KmerDiff(key_size(k_eff)), _lengthBytes(lengthBytes), _k_eff(k_eff) {}

  // `_lhs` and `_rhs` point to the least significant bytes of each k-mer
  KmerDiffClass operator()(const uint8_t *_lhs, const uint8_t *_rhs) const;

protected:
  uint32_t _lengthBytes;
  size_t _k_eff;
};

void adjacentDifference(TerminalBuffer &, sdsl::int_vector<2> &);
void adjacentDifference(TerminalRange, sdsl::int_vector<2> &);
sdsl::int_vector<2> adjacentDifference(TerminalBuffer &);
sdsl::int_vector<2> adjacentDifference(TerminalRange);
