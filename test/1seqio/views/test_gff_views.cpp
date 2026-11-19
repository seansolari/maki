#include "maki/build/io/filter.hpp"
#include "maki/build/io/gff.hpp"
#include <gtest/gtest.h>

#include "maki/core/utils/tempfile.hpp"
#include "view_common.hpp"

namespace {

std::vector<std::string> collectSequences(sequence_fragment_input_range auto &&range) {
  std::vector<std::string> result;

  for (const auto &fragment : range)
    result.push_back(toString(Dna4Sequence(fragment.data())));

  return result;
}

std::vector<uint64_t> collectIds(sequence_fragment_input_range auto &&range) {
  std::vector<uint64_t> ids;

  for (auto const &fragment : range)
    ids.push_back(fragment.id());

  return ids;
}

} // namespace

class Dna4GenomeTests : public testing::Test {
protected:
  Colours c;

  Dna4Genome makeGenome(const std::string &sequence) {
    tempio::TempFastaFile fnaFile(sequence);
    return parseGenome(fnaFile.path(), c, 0);
  }
};

TEST_F(Dna4GenomeTests, SingleContigProducesForwardAndReverseComplement) {
  auto genome = makeGenome(">chr1\nACGT\n");

  auto seqs = collectSequences(genome.fragments(3));
  auto ids = collectIds(genome.fragments(3));

  ASSERT_EQ(seqs.size(), 2u);

  EXPECT_EQ(seqs[0], "ACGT");
  EXPECT_EQ(seqs[1], revcomp("ACGT"));

  EXPECT_EQ(ids[0], 1u);
  EXPECT_EQ(ids[1], 1u);
}

TEST_F(Dna4GenomeTests, SplitsContigsOnN) {
  auto genome = makeGenome(">c1\nAAAANCCCCNGGG\n");

  auto seqs = collectSequences(genome.fragments(3));
  auto ids = collectIds(genome.fragments(3));

  ASSERT_EQ(seqs.size(), 6u);

  EXPECT_EQ(seqs[0], "AAAA");
  EXPECT_EQ(seqs[1], "CCCC");
  EXPECT_EQ(seqs[2], "GGG");
  
  EXPECT_EQ(seqs[3], revcomp("AAAA"));
  EXPECT_EQ(seqs[4], revcomp("CCCC"));
  EXPECT_EQ(seqs[5], revcomp("GGG"));

  EXPECT_EQ(ids[0], 1u);
  EXPECT_EQ(ids[1], 1u);
  EXPECT_EQ(ids[2], 1u);

  EXPECT_EQ(ids[3], 1u);
  EXPECT_EQ(ids[4], 1u);
  EXPECT_EQ(ids[5], 1u);
}

TEST_F(Dna4GenomeTests, LeadingAndTrailingNsIgnored) {
  auto genome = makeGenome(">c1\nNNNAAACCCNN\n");

  auto seqs = collectSequences(genome.fragments(3));

  ASSERT_EQ(seqs.size(), 2u);

  EXPECT_EQ(seqs[0], "AAACCC");
  EXPECT_EQ(seqs[1], revcomp("AAACCC"));
}

TEST_F(Dna4GenomeTests, ConsecutiveNsDoNotCreateEmptyFragments) {
  auto genome = makeGenome(">c1\nAAAANNNNNCCCC\n");

  auto seqs = collectSequences(genome.fragments(3));

  ASSERT_EQ(seqs.size(), 4u);

  EXPECT_EQ(seqs[0], "AAAA");
  EXPECT_EQ(seqs[1], "CCCC");
  
  EXPECT_EQ(seqs[2], revcomp("AAAA"));
  EXPECT_EQ(seqs[3], revcomp("CCCC"));
}

TEST_F(Dna4GenomeTests, TerminalsMatchFragmentsForSmallGenome) {
  auto genome = makeGenome(">c1\nAAAANCCCC\n");

  auto termSeqs = collectSequences(genome.terminals());
  auto fragSeqs = collectSequences(genome.fragments(3));

  EXPECT_EQ(termSeqs, fragSeqs);
}

TEST_F(Dna4GenomeTests, EmptyGenomeProducesNoFragments) {
  auto genome = makeGenome(">c1\n\n");

  EXPECT_TRUE(collectSequences(genome.fragments(3)).empty());
  EXPECT_TRUE(collectSequences(genome.terminals()).empty());
}