
#include "maki/build/kmers/buffers/kmers.hpp"
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using ::testing::ContainerEq;

using seqan3::operator""_dna4;

TEST(KmerBufferSort, SortNoChange) {
  constexpr uint8_t k = 3;

  KmerBuffer buffer(5, sizeof(uint32_t), k, k), temp(5, sizeof(uint32_t), k, k);

  buffer.insertFragment(buffer.begin(),
                        SequenceFragment{"AAAACCCC"_dna4, 1, false});
  buffer.sort(&temp);

  std::vector<std::string> expected{"AAA", "AAA", "AAC", "ACC", "CCC"},
      observed;

  for (auto it = buffer.begin(), end = buffer.end(); it != end; ++it) {
    observed.push_back(toString(it.readKmer(k)));
  }

  EXPECT_THAT(observed, ContainerEq(expected));
}

TEST(KmerBufferSort, SmallSort) {
  constexpr uint8_t k = 3;

  KmerBuffer buffer(5, sizeof(uint32_t), k), temp(5, sizeof(uint32_t), k);

  buffer.insertFragment(buffer.begin(),
                        SequenceFragment{"ATAACGT"_dna4, 5, true});
  buffer.sort(&temp);

  std::vector<std::string> expected{"TAA", "ATA", "AAC", "ACG", "CGT"},
      observed;

  for (auto it = buffer.begin(), end = buffer.end(); it != end; ++it) {
    observed.push_back(toString(it.readKmer(k)));
  }

  EXPECT_THAT(observed, ContainerEq(expected));
}

TEST(KmerBufferSort, LongSortWithColour) {
  constexpr uint8_t k = 6;

  KmerBuffer buffer(6, sizeof(uint32_t), k), temp(6, sizeof(uint32_t), k);

  buffer.insertFragment(buffer.begin(),
                        SequenceFragment{"ACGTACGTACGT"_dna4, 42u, false});
  buffer.sort(&temp);

  std::vector<std::string> expected{"TACGTA", "ACGTAC", "ACGTAC",
                                    "CGTACG", "CGTACG", "GTACGT"},
      observed;

  for (auto it = buffer.begin(), end = buffer.end(); it != end; ++it) {
    observed.push_back(toString(it.readKmer(k)));

    EXPECT_EQ(it.readValue().colour(), 42u);
  }

  EXPECT_THAT(observed, ContainerEq(expected));
}
