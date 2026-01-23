#include "maki/build/kmers/buffers/sort.hpp"
#include <algorithm>
#include <array>
#include <numeric>
#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <vector>

void lsdRadixSort(uint8_t *input, uint8_t *tmp, uint64_t n_recs,
                  uint32_t rec_size, uint32_t key_size, std::size_t chunksize) {
  if (n_recs <= 1)
    return;

  for (std::size_t j = 0; j < key_size; ++j) {
    // count
    std::size_t num_chunk = (n_recs + chunksize - 1) / chunksize;
    std::array<std::size_t, 256> initial_block{};
    std::vector<std::array<std::size_t, 256>> blocks(
        num_chunk, std::array<std::size_t, 256>{});

    tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, num_chunk, 1),
        [&](const tbb::blocked_range<std::size_t> &r) {
          for (std::size_t chk = r.begin(); chk != r.end(); ++chk) {
            auto &table = blocks[chk];
            const uint8_t *p = input + (chk * chunksize * rec_size),
                          *end = input +
                                 rec_size * (chk == (num_chunk - 1)
                                                 ? n_recs
                                                 : ((chk + 1) * chunksize));
            while (p != end) {
              ++table[*(p + j)];
              p += rec_size;
            }
          }
        });

    // partial sum
    for (std::size_t chk = 1; chk < num_chunk; ++chk) {
      for (std::size_t i = 0; i < 256; ++i) {
        blocks[chk][i] += blocks[chk - 1][i];
      }
    }

    std::array<std::size_t, 256> char_h1_blocks{};
    std::partial_sum(blocks.back().cbegin(), blocks.back().cend() - 1,
                     char_h1_blocks.begin() + 1);

    // unwind
    tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, num_chunk, 1),
        [&](const tbb::blocked_range<std::size_t> &r) {
          for (std::size_t chk = r.begin(); chk < r.end(); ++chk) {
            auto &block = chk == 0 ? initial_block : blocks[chk - 1];
            const uint8_t *p = input + (chk * chunksize * rec_size),
                          *end = input +
                                 rec_size * (chk == (num_chunk - 1)
                                                 ? n_recs
                                                 : ((chk + 1) * chunksize));
            while (p != end) {
              uint8_t key = *(p + j);
              std::copy(p, p + rec_size,
                        tmp +
                            ((char_h1_blocks[key] + block[key]++) * rec_size));
              p += rec_size;
            }
          }
        });

    std::swap(input, tmp);
  }
  return;
}
