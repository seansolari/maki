#include <gtest/gtest.h>
#include <vector>

#include "maki/build/io/filter.hpp"
#include "maki/build/io/gff.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/core/utils/tempfile.hpp"

class CreateSuffixPlanTest : public testing::Test {
protected:
  Colours c;

  Dna4Genome makeGenome(const std::string &sequence) {
    tempio::TempFastaFile fnaFile(sequence);
    return parseGenome(fnaFile.path(), c, 0);
  }

  std::vector<Dna4Genome>
  makeGenomes(const std::vector<std::string> &sequences) {
    std::vector<Dna4Genome> result;
    result.reserve(sequences.size());

    for (const std::string &sequence : sequences) {
      result.emplace_back(makeGenome(sequence));
    }

    return result;
  }
};

TEST_F(CreateSuffixPlanTest, SingleGenomeProducesExpectedCounts) {
  auto genomes = makeGenomes({">g1\nACGT\n"});

  auto tables = createSuffixPlan(genomes, 2, 2);

  ASSERT_EQ(tables.size(), 1u);

  EXPECT_EQ(tables[0][ShortSuffix("AC"_dna4)], 2u);
  EXPECT_EQ(tables[0][ShortSuffix("CG"_dna4)], 2u);
  EXPECT_EQ(tables[0][ShortSuffix("GT"_dna4)], 2u);
}

TEST_F(CreateSuffixPlanTest, IndependentCountsWithoutAccumulation) {
  auto genomes = makeGenomes({">g1\nACGT\n", ">g2\nAAAA\n"});

  auto tables = createSuffixPlan(genomes, 2, 2, 0, false);

  ASSERT_EQ(tables.size(), 2u);

  EXPECT_EQ(tables[0][ShortSuffix("AC"_dna4)], 2u);
  EXPECT_EQ(tables[0][ShortSuffix("CG"_dna4)], 2u);
  EXPECT_EQ(tables[0][ShortSuffix("GT"_dna4)], 2u);

  EXPECT_EQ(tables[1][ShortSuffix("AA"_dna4)], 3u);

  EXPECT_EQ(tables[1][ShortSuffix("AC"_dna4)], 0u);
}

TEST_F(CreateSuffixPlanTest, CountsAreAccumulatedAcrossGenomes) {
  auto genomes = makeGenomes({">g1\nACGT\n", ">g2\nAAAA\n"});

  auto tables = createSuffixPlan(genomes, 2, 2, 0, true);

  ASSERT_EQ(tables.size(), 2u);

  EXPECT_EQ(tables[0][ShortSuffix("AC"_dna4)], 2u);

  EXPECT_EQ(tables[1][ShortSuffix("AC"_dna4)], 2u);
  EXPECT_EQ(tables[1][ShortSuffix("CG"_dna4)], 2u);
  EXPECT_EQ(tables[1][ShortSuffix("GT"_dna4)], 2u);
  EXPECT_EQ(tables[1][ShortSuffix("AA"_dna4)], 3u);
}

TEST_F(CreateSuffixPlanTest, AccumulationAddsMatchingSuffixCounts) {
  auto genomes = makeGenomes({">g1\nAAAA\n", ">g2\nAAAA\n"});

  auto tables = createSuffixPlan(genomes, 2, 2, 0, true);

  EXPECT_EQ(tables[0][ShortSuffix("AA"_dna4)], 3u);
  EXPECT_EQ(tables[1][ShortSuffix("AA"_dna4)], 6u);
}

TEST_F(CreateSuffixPlanTest, EmptyGenomeProducesZeroCounts) {
  auto genomes = makeGenomes({">g1\n\n"});

  auto tables = createSuffixPlan(genomes, 2, 2);

  ASSERT_EQ(tables.size(), 1u);

  for (auto v : tables[0])
    EXPECT_EQ(v, 0u);
}

TEST_F(CreateSuffixPlanTest, EmptyGenomeDoesNotAffectAccumulation) {
  auto genomes = makeGenomes({">g1\nAAAA\n", ">g2\n\n", ">g3\nAAAA\n"});

  auto tables = createSuffixPlan(genomes, 2, 2, 0, true);

  EXPECT_EQ(tables[0][ShortSuffix("AA"_dna4)], 3u);

  EXPECT_EQ(tables[1][ShortSuffix("AA"_dna4)], 3u);

  EXPECT_EQ(tables[2][ShortSuffix("AA"_dna4)], 6u);
}

TEST_F(CreateSuffixPlanTest, AllEmptyGenomesProduceZeroTables) {
  auto genomes = makeGenomes({">g1\n\n", ">g2\n\n", ">g3\n\n"});

  auto tables = createSuffixPlan(genomes, 2, 2);

  ASSERT_EQ(tables.size(), 3u);

  for (auto const &table : tables) {
    for (auto v : table)
      EXPECT_EQ(v, 0u);
  }
}

TEST_F(CreateSuffixPlanTest, OffsetIsAppliedToAllGenomes) {
  auto genomes = makeGenomes({">g1\nACGT\n"});

  auto tables = createSuffixPlan(genomes, 2, 2, 2, false);

  EXPECT_EQ(tables[0][ShortSuffix("GT"_dna4)], 2u);

  EXPECT_EQ(tables[0][ShortSuffix("CG"_dna4)], 0u);
}

TEST_F(CreateSuffixPlanTest, GenomeShorterThanKProducesZeroCounts) {
  auto genomes = makeGenomes({">g1\nAC\n"});

  auto tables = createSuffixPlan(genomes, 4, 2);

  for (auto v : tables[0])
    EXPECT_EQ(v, 0u);
}
