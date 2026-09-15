#include "maki/build/kmers/buffers/kmers.hpp"
#include "maki/build/kmers/buffers/base_buffer.hpp"
#include "maki/build/kmers/buffers/sort.hpp"
#include "maki/maki.h"

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
