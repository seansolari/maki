#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "maki/build/io/fastq.hpp"
#include "maki/core/seq/io.hpp"
#include "maki/core/utils/tempfile.hpp"

using tempio::TempFastqFile;

namespace {

std::vector<std::string> readSequences(reads::Read const &read) {
  std::vector<std::string> result;

  for (auto const &fragment : read.fragments())
    result.push_back(toString(fragment.sequence));

  return result;
}

std::vector<std::string> chunkReadIds(reads::FastqDatasetRange const &range) {
  std::vector<std::string> ids;

  for (auto const &read : range.reads)
    ids.push_back(read.id());

  return ids;
}

} // namespace

TEST(FastqDatasetChunkView, SingleChunkContainsEntireDataset) {
  TempFastqFile fq("@r1\nAAAA\n+\nIIII\n"
                   "@r2\nCCCC\n+\nIIII\n"
                   "@r3\nGGGG\n+\nIIII\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  reads::FastqDatasetChunkView chunks(ds, 1);

  ASSERT_EQ(chunks.size(), 1u);

  auto range = chunks.range(0);

  ASSERT_EQ(range.reads.size(), 3u);

  EXPECT_EQ(range.reads[0].id(), "r1");
  EXPECT_EQ(range.reads[1].id(), "r2");
  EXPECT_EQ(range.reads[2].id(), "r3");
}

TEST(FastqDatasetChunkView, EvenDivision) {
  TempFastqFile fq("@r1\nAAAA\n+\nIIII\n"
                   "@r2\nCCCC\n+\nIIII\n"
                   "@r3\nGGGG\n+\nIIII\n"
                   "@r4\nTTTT\n+\nIIII\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  reads::FastqDatasetChunkView chunks(ds, 2);

  auto c0 = chunks.range(0);
  auto c1 = chunks.range(1);

  ASSERT_EQ(c0.reads.size(), 2u);
  ASSERT_EQ(c1.reads.size(), 2u);

  EXPECT_EQ(c0.reads[0].id(), "r1");
  EXPECT_EQ(c0.reads[1].id(), "r2");

  EXPECT_EQ(c1.reads[0].id(), "r3");
  EXPECT_EQ(c1.reads[1].id(), "r4");
}

TEST(FastqDatasetChunkView, UnevenDivision) {
  TempFastqFile fq("@r1\nAAAA\n+\nIIII\n"
                   "@r2\nCCCC\n+\nIIII\n"
                   "@r3\nGGGG\n+\nIIII\n"
                   "@r4\nTTTT\n+\nIIII\n"
                   "@r5\nACGT\n+\nIIII\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  reads::FastqDatasetChunkView chunks(ds, 2);

  auto c0 = chunks.range(0);
  auto c1 = chunks.range(1);

  EXPECT_EQ(c0.reads.size(), 3u);
  EXPECT_EQ(c1.reads.size(), 2u);
}

TEST(FastqDatasetChunkView, MoreChunksThanReads) {
  TempFastqFile fq("@r1\nAAAA\n+\nIIII\n"
                   "@r2\nCCCC\n+\nIIII\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  reads::FastqDatasetChunkView chunks(ds, 5);

  ASSERT_EQ(chunks.size(), 5u);

  std::size_t totalReads = 0;

  for (std::size_t i = 0; i < chunks.size(); ++i)
    totalReads += chunks.range(i).reads.size();

  EXPECT_EQ(totalReads, ds.size());
}

TEST(FastqDatasetChunkView, EveryReadAppearsExactlyOnce) {
  TempFastqFile fq("@r1\nAAAA\n+\nIIII\n"
                   "@r2\nCCCC\n+\nIIII\n"
                   "@r3\nGGGG\n+\nIIII\n"
                   "@r4\nTTTT\n+\nIIII\n"
                   "@r5\nACGT\n+\nIIII\n"
                   "@r6\nTGCA\n+\nIIII\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  reads::FastqDatasetChunkView chunks(ds, 4);

  std::vector<std::string> observedIds;

  for (std::size_t c = 0; c < chunks.size(); ++c) {
    auto range = chunks.range(c);

    for (auto const &read : range.reads)
      observedIds.push_back(read.id());
  }

  ASSERT_EQ(observedIds.size(), ds.size());

  std::vector<std::string> expected{"r1", "r2", "r3", "r4", "r5", "r6"};

  EXPECT_EQ(observedIds, expected);
}

TEST(FastqDatasetChunkView, FragmentSequencesPreserved) {
  TempFastqFile fq("@r1\nAAAA\n+\nIIII\n"
                   "@r2\nCCCC\n+\nIIII\n"
                   "@r3\nGGGG\n+\nIIII\n"
                   "@r4\nTTTT\n+\nIIII\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  reads::FastqDatasetChunkView chunks(ds, 2);

  auto c0 = chunks.range(0);
  auto c1 = chunks.range(1);

  EXPECT_EQ(readSequences(c0.reads[0]), std::vector<std::string>{"AAAA"});

  EXPECT_EQ(readSequences(c0.reads[1]), std::vector<std::string>{"CCCC"});

  EXPECT_EQ(readSequences(c1.reads[0]), std::vector<std::string>{"GGGG"});

  EXPECT_EQ(readSequences(c1.reads[1]), std::vector<std::string>{"TTTT"});
}

TEST(FastqDatasetChunkView, PreservesFragmentedReads) {
  TempFastqFile fq("@r1\nACGTNNNNAAAA\n+\nIIIIIIIIIIII\n"
                   "@r2\nCCCC\n+\nIIII\n"
                   "@r3\nGGGGNNNNTTTT\n+\nJJJJJJJJJJJJ\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  reads::FastqDatasetChunkView chunks(ds, 2);

  auto c0 = chunks.range(0);
  auto c1 = chunks.range(1);

  ASSERT_EQ(c0.reads[0].numFragments(), 2u);

  EXPECT_EQ(readSequences(c0.reads[0]),
            (std::vector<std::string>{"ACGT", "AAAA"}));

  EXPECT_EQ(readSequences(c1.reads[0]),
            (std::vector<std::string>{"GGGG", "TTTT"}));
}

TEST(FastqDatasetChunkView, ExpectedBoundariesForTenReadsThreeChunks) {
  std::ostringstream fastq;

  for (int i = 0; i < 10; ++i) {
    fastq << "@r" << i << "\n"
          << "ACGT\n+\nIIII\n";
  }

  TempFastqFile fq(fastq.str());

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  reads::FastqDatasetChunkView chunks(ds, 3);

  auto c0 = chunkReadIds(chunks.range(0));
  auto c1 = chunkReadIds(chunks.range(1));
  auto c2 = chunkReadIds(chunks.range(2));

  EXPECT_EQ(c0, (std::vector<std::string>{"r0", "r1", "r2", "r3"}));

  EXPECT_EQ(c1, (std::vector<std::string>{"r4", "r5", "r6"}));

  EXPECT_EQ(c2, (std::vector<std::string>{"r7", "r8", "r9"}));
}
