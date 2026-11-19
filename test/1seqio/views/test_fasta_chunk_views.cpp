#include "maki/build/io/filter.hpp"
#include "maki/core/seq/concepts.hpp"
#include "maki/core/utils/tempfile.hpp"
#include "view_common.hpp"
#include <gtest/gtest.h>

std::string toString(is_sequence_fragment auto &&rng) {
  return toString(Dna4Sequence(rng.data()));
}

std::string toString(is_sequence_fragment auto &rng) {
  return toString(Dna4Sequence(rng.data()));
}

class ChunkedDna4GenomeTest : public testing::Test {
protected:
  Colours c;

  ChunkedDna4Genome makeChunkedGenome(const std::string &sequence,
                                      std::size_t chunkSize,
                                      std::size_t overlap) {
    tempio::TempFastaFile fnaFile(sequence);
    return ChunkedDna4Genome(parseGenome(fnaFile.path(), c, 1), chunkSize,
                             overlap);
  }
};

TEST_F(ChunkedDna4GenomeTest, FirstChunkOfContigContainsTerminalSequence) {
  auto genome = makeChunkedGenome(">c1\nAAAACCCCGGGG\n", 4, 1);

  auto const &chunk = genome[0];

  EXPECT_FALSE(toString(chunk.terminals()).empty());
}

TEST_F(ChunkedDna4GenomeTest, InternalChunksReturnEmptyTerminalSequence) {
  auto genome = makeChunkedGenome(">c1\nAAAACCCCGGGGTTTT\n", 4, 1);

  ASSERT_EQ(genome.size(), 10u);

  for (size_t i = 0; i < genome.size(); ++i) {
    if ((i == 0) || (i == 5))
      EXPECT_FALSE(toString(genome[i].terminals()).empty());
    else
      EXPECT_TRUE(toString(genome[i].terminals()).empty());
  }
}

TEST_F(ChunkedDna4GenomeTest, EndChunkReportsTerminal) {
  auto genome = makeChunkedGenome(">c1\nAAAACCCCGGGGTTTT\n", 4, 1);

  EXPECT_TRUE(genome[genome.size() - 1].endIsTerminal());
}

TEST_F(ChunkedDna4GenomeTest, InternalChunksAreNotTerminal) {
  auto genome = makeChunkedGenome(">c1\nAAAACCCCGGGGTTTT\n", 4, 1);

  ASSERT_EQ(genome.size(), 10u);

  for (size_t i = 0; i < genome.size(); ++i) {
    if ((i == 4) || (i == 9))
      EXPECT_TRUE(genome[i].endIsTerminal());
    else
      EXPECT_FALSE(genome[i].endIsTerminal());
  }
}

TEST_F(ChunkedDna4GenomeTest, NCreatesSeparateTerminalRegions) {
  auto genome = makeChunkedGenome(">c1\nAAAAACCCCNNNNGGGGTTTT\n", 4, 1);

  bool foundTerminalAfterN = false;
  std::size_t found = 0;

  for (size_t i = 0; i < genome.size(); ++i) {
    auto seq = toString(genome[i].terminals());

    if (!seq.empty()) {
      ++found;

      if (seq.find("GGGG") != std::string::npos) {
        foundTerminalAfterN = true;
      }
    }
  }

  EXPECT_TRUE(foundTerminalAfterN);
  EXPECT_EQ(found, 4);
}

TEST_F(ChunkedDna4GenomeTest, ReverseComplementFragmentsExist) {
  auto genome = makeChunkedGenome(">c1\nAAAANCCCC\n", 4, 1);

  bool foundForward = false;
  bool foundReverse = false;

  for (size_t i = 0; i < genome.size(); ++i) {
    auto seq = toString(genome[i].fragments(3));

    if (seq == "AAAA")
      foundForward = true;

    if (seq == revcomp("AAAA"))
      foundReverse = true;
  }

  EXPECT_TRUE(foundForward);
  EXPECT_TRUE(foundReverse);
}

TEST_F(ChunkedDna4GenomeTest, ColourIdsPreservedAcrossReverseComplement) {
  auto genome = makeChunkedGenome(">c1\nAAAAANCCCCC\n", 4, 1);

  std::unordered_map<std::string, uint64_t> forwardIds;
  std::unordered_map<std::string, uint64_t> reverseIds;

  for (size_t i = 0; i < genome.size(); ++i) {
    auto fragment = genome[i].fragments(3);

    auto seq = toString(fragment);

    if (seq == "AAAAA")
      forwardIds["A"] = fragment.id();

    if (seq == revcomp("AAAAA"))
      reverseIds["A"] = fragment.id();
  }

  EXPECT_EQ(forwardIds["A"], reverseIds["A"]);
}