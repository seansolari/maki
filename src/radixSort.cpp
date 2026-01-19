#include <maki/radixSort.hpp>

namespace sort
{

    void touchArray([[maybe_unused]] uint8_t *tmp,
                    [[maybe_unused]] uint64_t n_recs,
                    [[maybe_unused]] uint32_t rec_size,
                    [[maybe_unused]] uint32_t n_threads)
    {
        return;
    }

    void radixSort(uint8_t *input,
                   uint8_t *tmp,
                   uint64_t n_recs,
                   uint32_t rec_size,
                   uint32_t key_size,
                   uint32_t n_threads)
    {
        lsdRadixSort(input, tmp, n_recs, rec_size, key_size, n_threads);
        return;
    }

    void lsdRadixSort(uint8_t *input,
                      uint8_t *tmp,
                      uint64_t n_recs,
                      uint32_t rec_size,
                      uint32_t key_size,
                      [[maybe_unused]] uint32_t n_threads)
    {
        if (n_recs <= 1)
            return;

        for (size_t j = 0; j < key_size; ++j)
        {
            // count
            size_t num_chunk = (n_recs + KMER_RADIX_CHUNKSIZE - 1) / KMER_RADIX_CHUNKSIZE;
            auto initial_block = std::array<size_t, 256>{};
            auto blocks = std::vector<std::array<size_t, 256>>(num_chunk, std::array<size_t, 256>{});

            tbb::parallel_for(tbb::blocked_range<size_t>(0, num_chunk, 1),
                              [&](const tbb::blocked_range<size_t> &r)
                              {
                                  for (size_t chk = r.begin(); chk != r.end(); ++chk)
                                  {
                                      auto &table = blocks[chk];
                                      const uint8_t
                                          *p = input + (chk * KMER_RADIX_CHUNKSIZE * rec_size),
                                          *end = input + rec_size * (chk == (num_chunk - 1) ? n_recs : ((chk + 1) * KMER_RADIX_CHUNKSIZE));
                                      while (p != end)
                                      {
                                          ++table[*(p + j)];
                                          p += rec_size;
                                      }
                                  }
                              });

            // partial sum
            for (size_t chk = 1; chk < num_chunk; ++chk)
            {
                for (size_t i = 0; i < 256; ++i)
                {
                    blocks[chk][i] += blocks[chk - 1][i];
                }
            }

            auto char_h1_blocks = std::array<size_t, 256>{};
            std::partial_sum(blocks.back().cbegin(),
                             blocks.back().cend() - 1,
                             char_h1_blocks.begin() + 1);

            // unwind
            tbb::parallel_for(tbb::blocked_range<size_t>(0, num_chunk, 1),
                              [&](const tbb::blocked_range<size_t> &r)
                              {
                                  for (size_t chk = r.begin(); chk < r.end(); ++chk)
                                  {
                                      auto &block = chk == 0 ? initial_block : blocks[chk - 1];
                                      const uint8_t
                                          *p = input + (chk * KMER_RADIX_CHUNKSIZE * rec_size),
                                          *end = input + rec_size * (chk == (num_chunk - 1) ? n_recs : ((chk + 1) * KMER_RADIX_CHUNKSIZE));
                                      while (p != end)
                                      {
                                          uint8_t key = *(p + j);
                                          std::copy(p, p + rec_size, tmp + ((char_h1_blocks[key] + block[key]++) * rec_size));
                                          p += rec_size;
                                      }
                                  }
                              });

            std::swap(input, tmp);
        }
        return;
    }

}
