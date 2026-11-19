#include "maki/build/io/filter.hpp"
#include "maki/core/utils/tempfile.hpp"
#include <gtest/gtest.h>

#include <string>
#include <vector>

std::vector<std::string> chunkSequences(ChunkedDna4Genome const &genome) {
  std::vector<std::string> result;

  for (size_t i = 0; i < genome.size(); ++i)
    result.push_back(toString(genome[i].sequence()));

  return result;
}

class GenomeChunking : public testing::Test {
protected:
  Colours c;

  Dna4Genome makeGenome(const std::string &sequence) {
    tempio::TempFastaFile fnaFile(sequence);
    return parseGenome(fnaFile.path(), c, 0);
  }
};

TEST_F(GenomeChunking, SingleContigSmallerThanChunkSize) {
  auto genome = makeGenome(">chr1\nACGTACGT\n");

  ChunkedDna4Genome chunks(std::move(genome), 100, 0);

  ASSERT_EQ(chunks.size(), 2u);

  EXPECT_EQ(toString(chunks[0].sequence()), "ACGTACGT");
  EXPECT_EQ(toString(chunks[1].sequence()), "ACGTACGT");

  EXPECT_EQ(chunks[0].id(), 1u);
  EXPECT_EQ(chunks[1].id(), 1u);
}

TEST_F(GenomeChunking, ExactMultipleWithoutOverlap) {
  auto genome = makeGenome(">chr1\nAAAACCCCGGGGTTTT\n");

  ChunkedDna4Genome chunks(std::move(genome), 4, 0);

  ASSERT_EQ(chunks.size(), 8u);

  EXPECT_EQ(toString(chunks[0].sequence()), "AAAA");
  EXPECT_EQ(toString(chunks[1].sequence()), "CCCC");
  EXPECT_EQ(toString(chunks[2].sequence()), "GGGG");
  EXPECT_EQ(toString(chunks[3].sequence()), "TTTT");

  for (int i = 0; i < 4; ++i)
    EXPECT_EQ(chunks[i].id(), 1u);
}

TEST_F(GenomeChunking, FinalChunkContainsRemainder) {
  auto genome = makeGenome(">chr1\nACGTACGTAC\n");

  ChunkedDna4Genome chunks(std::move(genome), 4, 0);

  ASSERT_EQ(chunks.size(), 6u);

  EXPECT_EQ(toString(chunks[0].sequence()), "ACGT");
  EXPECT_EQ(toString(chunks[1].sequence()), "ACGT");
  EXPECT_EQ(toString(chunks[2].sequence()), "AC");
}

TEST_F(GenomeChunking, FixedOverlapApplied) {
  auto genome = makeGenome(">chr1\nACGTACGTACGT\n");

  ChunkedDna4Genome chunks(std::move(genome), 5, 2);

  EXPECT_EQ(toString(chunks[0].sequence()), "ACGTA");

  EXPECT_EQ(toString(chunks[1].sequence()), "TACGT");

  EXPECT_EQ(toString(chunks[2].sequence()), "GTACG");

  EXPECT_EQ(toString(chunks[3].sequence()), "CGT");
}

TEST_F(GenomeChunking, OverlapOneLessThanChunkSize) {
  auto genome = makeGenome(">chr1\nACGTACGT\n");

  ChunkedDna4Genome chunks(std::move(genome), 4, 3);

  EXPECT_GT(chunks.size(), 4u);

  EXPECT_EQ(toString(chunks[0].sequence()), "ACGT");
  EXPECT_EQ(toString(chunks[1].sequence()), "CGTA");
}

TEST_F(GenomeChunking, ContigIdsIncreasePerContig) {
  auto genome = makeGenome(">c1\nAAAAAA\n"
                           ">c2\nCCCCCC\n");

  ChunkedDna4Genome chunks(std::move(genome), 3, 0);

  EXPECT_EQ(chunks.size(), 8);

  EXPECT_EQ(chunks[0].id(), 1u);
  EXPECT_EQ(chunks[1].id(), 1u);
  EXPECT_EQ(chunks[2].id(), 1u);
  EXPECT_EQ(chunks[3].id(), 1u);

  EXPECT_EQ(chunks[4].id(), 2u);
  EXPECT_EQ(chunks[5].id(), 2u);
  EXPECT_EQ(chunks[6].id(), 2u);
  EXPECT_EQ(chunks[7].id(), 2u);
}

TEST_F(GenomeChunking, NsCreateIndependentFragments) {
  auto genome = makeGenome(">chr1\n"
                           "AAAAANNNCCCC\n");

  ChunkedDna4Genome chunks(std::move(genome), 100, 0);

  ASSERT_EQ(chunks.size(), 4u);

  // Check sequences
  std::vector<std::string> sequences = chunkSequences(chunks),
    expectedSequences = {"AAAAA", "CCCC", "GGGG", "TTTTT"};

  std::sort(sequences.begin(), sequences.end());

  EXPECT_EQ(sequences, expectedSequences);

  // Check IDs
  for (std::size_t i = 0; i < chunks.size(); ++i)
    EXPECT_EQ(chunks[i].id(), 1u);
}

TEST_F(GenomeChunking, ContigFragmentsShareId) {
  auto genome = makeGenome(">chr1\n"
                           "AAAANNNNCCCCNNNNGGGG\n");

  ChunkedDna4Genome chunks(std::move(genome), 2, 0);

  uint64_t expected = chunks[0].id();

  for (size_t i = 0; i < chunks.size(); ++i) {
    EXPECT_EQ(chunks[i].id(), expected);
  }
}

TEST_F(GenomeChunking, MultipleContigsReceiveDifferentIds) {
  auto genome = makeGenome(">chr1\nAAAA\n"
                           ">chr2\nCCCC\n"
                           ">chr3\nGGGG\n");

  ChunkedDna4Genome chunks(std::move(genome), 100, 0);

  std::set<uint64_t> ids;

  for (size_t i = 0; i < chunks.size(); ++i)
    ids.insert(chunks[i].id());

  EXPECT_EQ(ids, (std::set<uint64_t>{1, 2, 3}));
}

TEST_F(GenomeChunking, ChunksRemainGroupedByContig) {
  auto genome = makeGenome(">c1\nAAAAAAA\n"
                           ">c2\nCCCCCCC\n");

  ChunkedDna4Genome chunks(std::move(genome), 3, 0);

  bool seenContig2 = false;

  for (size_t i = 0; i < chunks.size(); ++i) {
    if (chunks[i].id() == 2) {
      seenContig2 = true;
    }

    if (seenContig2) {
      EXPECT_EQ(chunks[i].id(), 2u);
    }
  }

  EXPECT_TRUE(seenContig2);
}

TEST_F(GenomeChunking, ChunksReconstructOriginalSequence) {
  auto genome = makeGenome(">chr1\n"
                           "ACGTACGTACGTACGT\n");

  constexpr size_t chunkSize = 5;
  constexpr size_t overlap = 2;

  ChunkedDna4Genome chunks(std::move(genome), chunkSize, overlap);

  std::string reconstructed;

  reconstructed += toString(chunks[0].sequence());

  for (size_t i = 1; i < chunks.size()/2 /*forward chunks only*/; ++i) {
    auto seq = toString(chunks[i].sequence());

    reconstructed += seq.substr(overlap);
  }

  EXPECT_EQ(reconstructed, "ACGTACGTACGTACGT");
}
