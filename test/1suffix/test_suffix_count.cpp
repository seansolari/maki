#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/core/seq/concepts.hpp"
#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

using seqan3::operator""_dna4;

TEST(SuffixTableCountKmerSuffixes, CountsSingleTerminalSuffix) {
  SuffixTable table;
  table.resize(2);

  auto seq = Dna4Sequence("AC"_dna4);
  SequenceFragment<Dna4Sequence> fragment(std::move(seq), 0, true);

  table.countKmerSuffixes(fragment, 2, 0);

  EXPECT_EQ(table[ShortSuffix("AC"_dna4)], 1u);
  EXPECT_EQ(table.size(), 16u);
}

TEST(SuffixTableCountKmerSuffixes,
     NonTerminalFragmentOfLengthKProducesNoCounts) {
  SuffixTable table;
  table.resize(2);

  auto seq = Dna4Sequence("AC"_dna4);
  SequenceFragment<Dna4Sequence> fragment(std::move(seq), 0, false);

  table.countKmerSuffixes(fragment, 2, 0);

  for (auto v : table)
    EXPECT_EQ(v, 0u);
}

TEST(SuffixTableCountKmerSuffixes, CountsAllTerminalWindows) {
  SuffixTable table;
  table.resize(2);

  auto seq = Dna4Sequence("ACGT"_dna4);
  SequenceFragment<Dna4Sequence> fragment(std::move(seq), 0, true);

  table.countKmerSuffixes(fragment, 2, 0);

  EXPECT_EQ(table[ShortSuffix("AC"_dna4)], 1u);
  EXPECT_EQ(table[ShortSuffix("CG"_dna4)], 1u);
  EXPECT_EQ(table[ShortSuffix("GT"_dna4)], 1u);
}

TEST(SuffixTableCountKmerSuffixes, ExcludesLastWindowForInteriorFragment) {
  SuffixTable table;
  table.resize(2);

  auto seq = Dna4Sequence("ACGT"_dna4);
  SequenceFragment<Dna4Sequence> fragment(std::move(seq), 0, false);

  table.countKmerSuffixes(fragment, 2, 0);

  EXPECT_EQ(table[ShortSuffix("AC"_dna4)], 1u);
  EXPECT_EQ(table[ShortSuffix("CG"_dna4)], 1u);

  EXPECT_EQ(table[ShortSuffix("GT"_dna4)], 0u);
}

TEST(SuffixTableCountKmerSuffixes, OffsetSkipsLeadingWindows) {
  SuffixTable table;
  table.resize(2);

  auto seq = Dna4Sequence("ACGTAC"_dna4);
  SequenceFragment<Dna4Sequence> fragment(std::move(seq), 0, true);

  table.countKmerSuffixes(fragment, 2, 2);

  EXPECT_EQ(table[ShortSuffix("GT"_dna4)], 1u);
  EXPECT_EQ(table[ShortSuffix("TA"_dna4)], 1u);
  EXPECT_EQ(table[ShortSuffix("AC"_dna4)], 1u);

  EXPECT_EQ(table[ShortSuffix("CG"_dna4)], 0u);
}

TEST(SuffixTableCountKmerSuffixes, RepeatedSequenceAccumulatesCounts) {
  SuffixTable table;
  table.resize(2);

  auto seq = Dna4Sequence("AAAAA"_dna4);
  SequenceFragment<Dna4Sequence> fragment(std::move(seq), 0, true);

  table.countKmerSuffixes(fragment, 2, 0);

  EXPECT_EQ(table[ShortSuffix("AA"_dna4)], 4u);
}

TEST(SuffixTableCountKmerSuffixes, LengthBelowKIsIgnored) {
  SuffixTable table;
  table.resize(2);

  auto seq = Dna4Sequence("ACG"_dna4);
  SequenceFragment<Dna4Sequence> fragment(std::move(seq), 0, true);

  table.countKmerSuffixes(fragment, 4, 0);

  for (auto v : table)
    EXPECT_EQ(v, 0u);
}

TEST(SuffixTableCountKmerSuffixes, InteriorFragmentMayFailKFilter) {
  SuffixTable table;
  table.resize(2);

  auto seq = Dna4Sequence("ACGT"_dna4);
  SequenceFragment<Dna4Sequence> fragment(std::move(seq), 0, false);

  table.countKmerSuffixes(fragment, 4, 0);

  for (auto v : table)
    EXPECT_EQ(v, 0u);
}

TEST(SuffixTableCountKmerSuffixes, OffsetCanReduceToSingleWindow) {
  SuffixTable table;
  table.resize(2);

  auto seq = Dna4Sequence("ACGT"_dna4);
  SequenceFragment<Dna4Sequence> fragment(std::move(seq), 0, true);

  table.countKmerSuffixes(fragment, 2, 2);

  EXPECT_EQ(table[ShortSuffix("GT"_dna4)], 1u);

  EXPECT_EQ(table[ShortSuffix("AC"_dna4)], 0u);
  EXPECT_EQ(table[ShortSuffix("CG"_dna4)], 0u);
}
