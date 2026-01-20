#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <numeric>
#include <vector>
#include <oneapi/tbb.h>
#include <maki/maki.h>

namespace sort
{
    constexpr uint32_t ALIGNMENT = 256U;

    void touchArray([[maybe_unused]] uint8_t *tmp,
                    [[maybe_unused]] uint64_t n_recs,
                    [[maybe_unused]] uint32_t rec_size,
                    [[maybe_unused]] uint32_t n_threads);

    void radixSort(uint8_t *input,
                   uint8_t *tmp,
                   uint64_t n_recs,
                   uint32_t rec_size,
                   uint32_t key_size,
                   uint32_t n_threads);

    void lsdRadixSort(uint8_t *input,
                      uint8_t *tmp,
                      uint64_t n_recs,
                      uint32_t rec_size,
                      uint32_t key_size,
                      uint32_t n_threads);

    template <typename Buffer_>
    Buffer_* sort(Buffer_ *full, Buffer_ *temp, uint32_t n_threads)
    {
        radixSort(full->data(), temp->data(), full->size(), full->recordBytes(), full->keyBytes(), n_threads);
        return full->keyBytes() % 2
            ? temp
            : full;
    }

}
