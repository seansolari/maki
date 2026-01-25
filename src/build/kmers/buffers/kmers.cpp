#include "maki/build/kmers/buffers/kmers.hpp"
#include "maki/build/kmers/buffers/base_buffer.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/build/kmers/buffers/sort.hpp"
#include "maki/maki.h"

KmerBuffer::iterator KmerBuffer::insert(iterator it,
                                        Dna4SequenceConstIter seqIter,
                                        Dna4SequenceConstIter endIter,
                                        uint64_t colour, bool endIsTerminal) {
  Kmer kmer(k, seqIter);
  seqIter += k;

  uint8_t edge;

  while (seqIter != endIter) {
    edge = parsing::dna4ToShort(*seqIter++);

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

KmerBuffer::iterator KmerBuffer::insert(iterator it,
                                        Dna4SequenceConstIter seqIter,
                                        Dna4SequenceConstIter endIter,
                                        uint64_t colour, bool endIsTerminal,
                                        ShortSuffix key) {
  uint8_t s = key.size();
  seqIter += k - s;

  // initialise
  ShortSuffix suffix(s, seqIter);
  seqIter += s;

  uint64_t edge;

  // fill
  while (seqIter != endIter) {
    edge = parsing::dna4ToLong(*seqIter);

    if (suffix == key) {
      it.writeKmer(seqIter - k, k_eff);
      it.writeValue(BufferValue(colour, edge));
      ++it;
    }

    suffix.roll(edge);
    ++seqIter;
  }

  // terminal edge
  if ((suffix == key) && endIsTerminal) {
    it.writeKmer(seqIter - k, k_eff);
    it.writeValue(BufferValue((uint64_t)0u, terminalEdge));
    ++it;
  }

  return it;
}

void KmerBuffer::fill(const std::vector<const SequenceContainer *> &data_,
                      const std::vector<SuffixTable> &blocks_,
                      ShortSuffix sfx_) {
  oneapi::tbb::parallel_for(
      (std::size_t)0, data_.size(), (std::size_t)1, [&](std::size_t i) {
        auto it = at(i == 0 ? 0 : blocks_[i - 1][sfx_]);
        for (auto seq : data_[i]->terminals()) {
          it = insert(it, seq.begin(), seq.begin() + k, seq.id(),
                      seq.endIsTerminal(), sfx_);
        }
      });
}

void KmerBuffer::sort(KmerBuffer *temp) {
  KmerBuffer *result = radixSort(this, temp, KMER_RADIX_CHUNKSIZE);
  if (result != this)
    _data.swap(result->_data);
}

void KmerBuffer::sort() {
  KmerBuffer temp(_num_records, _value_bytes, k, k_eff);
  sort(&temp);
}

KmerDiffClass KmerDiff::operator()(const uint8_t *_lhs,
                                   const uint8_t *_rhs) const {
  // - check most significant bytes
  if ((_seqWidth > 1) && (!std::equal(_lhs + 1, _lhs + _seqWidth, _rhs + 1))) {
    return BW_0_K;
  } else {
    // - check least significant byte
    uint8_t _lhs_v = *_lhs, _rhs_v = *_rhs;

    if (_lhs_v == _rhs_v)
      return IS_0;
    else if ((_lhs_v >> 2) == (_rhs_v >> 2))
      return IS_K;
    else
      return BW_0_K;
  }
}

KmerDiffClass MaskedBytesDiff::operator()(const uint8_t *l_,
                                          const uint8_t *r_) const {
  uint8_t msk = _msb_mask;

  for (int i = _width - 1; i >= 0; --i) {
    uint8_t l_val = *(l_ + i), r_val = *(r_ + i);

    l_val &= msk;

    if (l_val != r_val) {
      if (i == 0) {
        if ((l_val >> 2) == (r_val >> 2)) {
          return IS_K;
        }
      }
      return BW_0_K;
    }

    msk = 0xFFu;
  }

  return IS_0;
}

void adjacentDifference(KmerBuffer &buffer, sdsl::int_vector<2> &arr) {
  adjacentDifference(buffer.begin(), buffer.end(), KmerDiff{buffer.keyBytes()},
                     arr);
}

sdsl::int_vector<2> adjacentDifference(KmerBuffer &buffer) {
  return adjacentDifference(buffer.begin(), buffer.end(),
                            KmerDiff{buffer.keyBytes()});
}
