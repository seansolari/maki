
#pragma once
#include <cstdint>

constexpr uint32_t ALIGNMENT = 256U;

void lsdRadixSort(uint8_t *input, uint8_t *tmp, uint64_t n_recs,
                  uint32_t rec_size, uint32_t key_size, std::size_t chunksize);

template <typename Buffer_>
Buffer_ *radixSort(Buffer_ *full, Buffer_ *temp, std::size_t grainsize) {
  lsdRadixSort(full->data(), temp->data(), full->size(), full->recordBytes(),
               full->keyBytes(), grainsize);
  return full->keyBytes() % 2 ? temp : full;
}
