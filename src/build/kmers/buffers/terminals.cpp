
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/build/kmers/buffers/base_buffer.hpp"
#include "maki/build/kmers/buffers/sort.hpp"
#include "maki/maki.h"
#include <execution>
#include <optional>

void TerminalBuffer::shrink(size_t new_size_) {
  assert(new_size_ <= _num_records);
  _num_records = new_size_;
}

TerminalBuffer::iterator TerminalBuffer::insert(iterator it,
                                                Dna4SequenceConstIter seqIter,
                                                Dna4SequenceConstIter seqEnd,
                                                bool endIsTerminal) {
  LongSuffix tl(k_eff /* here, will be same as k */);
  uint8_t edge;

  while (seqIter != seqEnd) {
    edge = parsing::dna4ToShort(*seqIter++);

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

TerminalBuffer::iterator TerminalBuffer::insert(iterator it,
                                                Dna4SequenceConstIter seqIter,
                                                Dna4SequenceConstIter seqEnd,
                                                bool endIsTerminal,
                                                ShortSuffix key) {
  uint8_t s = key.size();
  size_t distanceFromStart = 0, effTerminalSize;

  // initialise

  ShortSuffix suffix(s, seqIter);
  seqIter += s;

  // fill

  uint64_t edge;

  while (seqIter != seqEnd) {
    edge = parsing::dna4ToLong(*seqIter);

    if (suffix == key) {
      effTerminalSize = std::min(distanceFromStart, (size_t)k_eff);
      it.writeKey(seqIter - (s + 1), effTerminalSize);
      it.writeEdge(edge);
      ++it;
    }

    suffix.roll(edge);
    ++seqIter;
    ++distanceFromStart;
  }

  // terminal edge

  if ((suffix == key) && endIsTerminal) {
    effTerminalSize = std::min(distanceFromStart, (size_t)k_eff);
    it.writeKey(seqIter - (s + 1), effTerminalSize);
    it.writeEdge(terminalEdge);
    ++it;
  }

  return it;
}

void TerminalBuffer::insert(LongSuffixGate &lock, Dna4SequenceConstIter seqIter,
                            Dna4SequenceConstIter seqEnd, bool endIsTerminal) {
  LongSuffix tl(k_eff /* here, will be same as k */);

  uint8_t edge;
  std::optional<std::size_t> p = std::nullopt;

  while (seqIter != seqEnd) {
    edge = parsing::dna4ToShort(*seqIter++);

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

void TerminalBuffer::fill(const std::vector<const SequenceContainer *> &data_,
                          std::vector<size_t> const &blocks_) {
  oneapi::tbb::parallel_for(
      (std::size_t)0, data_.size(), (std::size_t)1, [&](std::size_t i) {
        auto it = at(i == 0 ? 0 : blocks_[i - 1]);
        for (auto seq : data_[i]->terminals()) {
          it = insert(it, seq.begin(), seq.begin() + k, false);
        }
      });
}

void TerminalBuffer::fill(const std::vector<const SequenceContainer *> &data_,
                          LongSuffixGate &lock) {
  oneapi::tbb::parallel_for(
      (std::size_t)0, data_.size(), (std::size_t)1, [&](std::size_t i) {
        for (auto seq : data_[i]->terminals()) {
          insert(lock, seq.begin(), seq.begin() + k, false);
        }
      });
}

void TerminalBuffer::fill(const std::vector<const SequenceContainer *> &data_,
                          const std::vector<SuffixTable> &blocks_,
                          ShortSuffix sfx_) {
  oneapi::tbb::parallel_for(
      (std::size_t)0, data_.size(), (std::size_t)1, [&](std::size_t i) {
        auto it = at(i == 0 ? 0 : blocks_[i - 1][sfx_]);
        for (auto seq : data_[i]->fragments(k)) {
          it = insert(it, seq.begin(), seq.end(), seq.endIsTerminal(), sfx_);
        }
      });
}

TerminalRange TerminalBuffer::asRange() const {
  return TerminalRange(k_eff, _lengthBytes, constBegin(), constEnd());
}

TerminalConstIter TerminalRange::lowerBound(ndim::Span<uint8_t *> trg_,
                                            uint8_t k_) const {
  assert(k_ <= _k_eff);
  return std::lower_bound(constBegin(), constEnd(), trg_,
                          MaskedAlignedBytesLessThan<LHS>(k_, _lengthBytes));
}

TerminalConstIter TerminalRange::lowerBound(ndim::Span<uint8_t *> trg_) const {
  return lowerBound(trg_, _k_eff);
}

TerminalConstIter TerminalRange::upperBound(ndim::Span<uint8_t *> trg_,
                                            uint8_t k_) const {
  assert(k_ <= _k_eff);
  /**
   * Note that in `std::upper_bound`, `trg_` will be supplied as the LHS
   * argument, and `*Iter` (i.e. the terminal) will be supplied as the RHS
   * (which needs to be masked).
   */
  return std::upper_bound(constBegin(), constEnd(), trg_,
                          MaskedAlignedBytesLessThan<RHS>(k_, _lengthBytes));
}

TerminalConstIter TerminalRange::upperBound(ndim::Span<uint8_t *> trg_) const {
  return upperBound(trg_, _k_eff);
}

TerminalRange TerminalRange::endsWith(ShortSuffix sfx_) const {
  size_t s_ = sfx_.size();
  UnalignedBytesLessThan<true> LessThan(_k_eff, (uint8_t)s_, _lengthBytes);

  auto _begin =
           std::lower_bound(constBegin(), constEnd(), sfx_.data(), LessThan),
       _end = std::upper_bound(constBegin(), constEnd(), sfx_.data(), LessThan);

  while ((_begin < _end) && (_begin.readSize() < s_))
    ++_begin;

  return TerminalRange(_k_eff, _lengthBytes, _begin, _end);
}

TerminalRange TerminalRange::retrieve(ShortSuffix sfx_) const {
  size_t s_ = sfx_.size();
  TerminalIsSuffix Eq(_k_eff, (uint8_t)s_, _lengthBytes);

  auto _end = constEnd(),
       _lhs = std::lower_bound(
           constBegin(), _end, sfx_.data(),
           UnalignedBytesLessThan<false>(_k_eff, (uint8_t)s_, _lengthBytes));

  while ((_lhs < _end) && Eq(_lhs, sfx_) && (_lhs.readSize() < s_))
    ++_lhs;

  auto _rhs = _lhs;

  while ((_rhs < _end) && Eq(_rhs, sfx_) && (_rhs.readSize() == s_))
    ++_rhs;

  return TerminalRange(_k_eff, _lengthBytes, _lhs, _rhs);
}

bool TerminalsLessThan::operator()(size_t i_, size_t j_) const {
  const uint8_t *ip_ = _buffer.cdata(i_), *jp_ = _buffer.cdata(j_);

  // compare sequence
  switch (_rng_cmp(ip_, jp_, _keyWidth)) {
  case -1:
    return true;
  case 1:
    return false;
  default:
    // compare edge
    return *(ip_ + _valueOffset) < *(jp_ + _valueOffset);
  }
}

int TerminalsLessThan::_rng_cmp(const uint8_t *l_, const uint8_t *r_,
                                uint32_t width) {
  for (const uint8_t *lhs = l_ + width - 1, *rhs = r_ + width - 1;
       lhs != l_ - 1; (void)--lhs, (void)--rhs) {
    uint8_t lhs_v = *lhs, rhs_v = *rhs;
    if (lhs_v < rhs_v) {
      return -1;
    } else if (lhs_v > rhs_v) {
      return 1;
    }
  }
  return 0;
}

bool TerminalIsSuffix::operator()(TerminalConstIter it,
                                  ShortSuffix const &sfx) const {
  return operator()(it.cdata(), sfx.data());
}

bool TerminalIsSuffix::operator()(uint8_t const *trg_, uint64_t sfx_) const {
  uint8_t const *trg_aln_ = trg_ + _trg_width - _sfx_width;
  sfx_ <<= 2 * _sfx_shift;

  for (long i = 0; i < _sfx_width; ++i) {
    uint8_t trg_v = *(trg_aln_ + i), sfx_v = (sfx_ >> (8 * i)) & 0xFFu;
    if (trg_v != sfx_v)
      return false;
  }

  for (long i = _trg_width - _sfx_width - _offset; i > 0; --i) {
    if (*(trg_aln_ - i) != 0u)
      return false;
  }

  return true;
}

void TerminalBuffer::sort(TerminalBuffer *temp) {
  TerminalBuffer *result = radixSort(this, temp, KMER_RADIX_CHUNKSIZE);
  if (result != this)
    _data.swap(result->_data);
}

void TerminalBuffer::sort() {
  TerminalBuffer temp(_record_capacity, k, k_eff, _lengthBytes);
  temp.shrink(_num_records);
  sort(&temp);
}

TerminalBuffer TerminalBuffer::OOPsort(bool makeUnique) const {
  std::vector<size_t> indices(_num_records); // capacity constructor
  std::iota(indices.begin(), indices.end(), 0ul);

  // sort indices
  std::sort(std::execution::par_unseq, indices.begin(), indices.end(),
            TerminalsLessThan{*this});
  auto indices_end = indices.end();

  if (makeUnique) {
    // get rid of duplicate values
    indices_end = std::unique(
        std::execution::par_unseq, indices.begin(), indices_end,
        [&](size_t i, size_t j) -> bool { return *at(i) == *at(j); });
  }

  // take values in order
  return TerminalBuffer(*this, indices.begin(), indices_end);
}

void TerminalBuffer::unique() {
  auto __first = begin(), __last = end();

  __first = std::adjacent_find(__first, __last);
  if (__first == __last)
    return;

  auto __out = __first;
  ++__first;
  while (++__first != __last)
    if (*__out != *__first)
      (*++__out).write(*__first);

  shrink(++__out - begin());
}

KmerDiffClass TerminalDiff::operator()(const uint8_t *_lhs,
                                       const uint8_t *_rhs) const {
  KmerDiffClass dseq =
      KmerDiff::operator()(_lhs + _lengthBytes, _rhs + _lengthBytes);
  if (dseq == BW_0_K) {
    return dseq;
  } else {
    // Here, `dseq` is either `IS_0` or `IS_K`. We need the size of each k-mer.
    size_t small = deser(_lhs, _lengthBytes), big = deser(_rhs, _lengthBytes);
    if (small > big)
      std::swap(small, big);

    // If it's `IS_0`, we need to compare the sizes.
    //    - If they're the same size, it's fine.
    //    - If they differ by 1, and the larger sequence is length `k_eff`, then
    //    it's actually `IS_K`.
    //    - Otherwise, its `BW_0_K`.
    // If its `IS_K`, we need to compare the sizes.
    //    - If they're the same size, or they differ by 1, it's fine.
    //    - If they differ by more than 1, it's actually `BW_0_K`.
    //
    if (dseq == IS_0) {
      if (small == big)
        return dseq;
      else if ((big - small == 1u) && (big == _k_eff))
        return IS_K;
      else
        return BW_0_K;
    } else { // dseq == IS_K
      if (big - small > 1u)
        return BW_0_K;
      else
        return IS_K;
    }
  }
}

void adjacentDifference(TerminalBuffer &buffer, sdsl::int_vector<2> &arr) {
  adjacentDifference(buffer.constBegin(), buffer.constEnd(),
                     TerminalDiff(buffer.lengthBytes(), buffer.getEffK()), arr);
}

void adjacentDifference(TerminalRange buffer, sdsl::int_vector<2> &arr) {
  adjacentDifference(buffer.constBegin(), buffer.constEnd(),
                     TerminalDiff(buffer.lengthBytes(), buffer.getEffK()), arr);
}

sdsl::int_vector<2> adjacentDifference(TerminalBuffer &buffer) {
  return adjacentDifference(
      buffer.constBegin(), buffer.constEnd(),
      TerminalDiff(buffer.lengthBytes(), buffer.getEffK()));
}

sdsl::int_vector<2> adjacentDifference(TerminalRange buffer) {
  return adjacentDifference(
      buffer.constBegin(), buffer.constEnd(),
      TerminalDiff(buffer.lengthBytes(), buffer.getEffK()));
}
