#include <gtest/gtest.h>

#include <string>

#include "maki/build/io/fastq.hpp"
#include "maki/core/seq/io.hpp"
#include "maki/core/utils/tempfile.hpp"

using tempio::TempFastqFile;

void expectFragment(reads::ReadFragment const &fragment,
                    std::string const &sequence, std::string const &qualities) {
  EXPECT_EQ(toString(fragment.sequence), sequence);
  EXPECT_EQ(toString(fragment.qualities), qualities);
}

TEST(FastqDatasetUnpaired, ParsesSingleRead) {
  TempFastqFile fq("@read1\n"
                   "ACGTACGT\n"
                   "+\n"
                   "IIIIIIII\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  ASSERT_EQ(ds.size(), 1u);

  auto const &read = ds[0];

  EXPECT_EQ(read.id(), "read1");
  ASSERT_EQ(read.numFragments(), 1u);

  expectFragment(read.fragments()[0], "ACGTACGT", "IIIIIIII");
}

TEST(FastqDatasetUnpaired, ParsesMultipleReads) {
  TempFastqFile fq("@r1\nACGT\n+\nIIII\n"
                   "@r2\nTTTT\n+\nJJJJ\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  ASSERT_EQ(ds.size(), 2u);

  EXPECT_EQ(ds[0].id(), "r1");
  EXPECT_EQ(ds[1].id(), "r2");
}

TEST(FastqDatasetUnpaired, SplitsReadOnSingleN) {
  TempFastqFile fq("@read1\n"
                   "ACGTNACGT\n"
                   "+\n"
                   "IIIIIIIII\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  auto const &read = ds[0];

  ASSERT_EQ(read.numFragments(), 2u);

  expectFragment(read.fragments()[0], "ACGT", "IIII");

  expectFragment(read.fragments()[1], "ACGT", "IIII");
}

TEST(FastqDatasetUnpaired, ConsecutiveNsCreateSingleGap) {
  TempFastqFile fq("@r1\n"
                   "ACGTNNNNAC\n"
                   "+\n"
                   "JJJJJJJJJJ\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  auto const &read = ds[0];

  ASSERT_EQ(read.numFragments(), 2u);

  expectFragment(read.fragments()[0], "ACGT", "JJJJ");
  expectFragment(read.fragments()[1], "AC", "JJ");
}

TEST(FastqDatasetUnpaired, ReadContainingOnlyNsProducesNoFragments) {
  TempFastqFile fq("@r1\n"
                   "NNNNNN\n"
                   "+\n"
                   "IIIIII\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  ASSERT_EQ(ds.size(), 1u);

  EXPECT_TRUE(ds[0].empty());
  EXPECT_EQ(ds[0].numFragments(), 0u);
}

TEST(FastqDatasetUnpaired, SplitsOnArbitraryNonAcgtSymbols) {
  TempFastqFile fq("@r1\n"
                   "ACGTRACGTYTT\n"
                   "+\n"
                   "IIIIIIJJIIII\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  auto const &read = ds[0];

  ASSERT_EQ(read.numFragments(), 3u);

  expectFragment(read.fragments()[0], "ACGT", "IIII");
  expectFragment(read.fragments()[1], "ACGT", "IJJI");
  expectFragment(read.fragments()[2], "TT", "II");
}

TEST(FastqDatasetUnpaired, FragmentQualitiesRemainAligned) {
  TempFastqFile fq("@r1\n"
                   "ACGTNAC\n"
                   "+\n"
                   "ABCDABC\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  auto const &read = ds[0];

  ASSERT_EQ(read.numFragments(), 2u);

  expectFragment(read.fragments()[0], "ACGT", "ABCD");

  expectFragment(read.fragments()[1], "AC", "BC");
}

TEST(FastqDatasetUnpaired, RejectsSequenceQualityLengthMismatch) {
  TempFastqFile fq("@r1\n"
                   "ACGT\n"
                   "+\n"
                   "III\n");

  EXPECT_THROW(reads::FastqDataset::fromUnpaired(fq.path()),
               reads::FastqDatasetError);
}

TEST(FastqDatasetUnpaired, RejectsInvalidQualityCharacter) {
  TempFastqFile fq("@r1\n"
                   "ACGT\n"
                   "+\n"
                   "\x01\x01\x01\x01");

  EXPECT_THROW(reads::FastqDataset::fromUnpaired(fq.path()),
               reads::FastqDatasetError);
}

TEST(FastqDatasetUnpaired, RejectsIncompleteRecord) {
  TempFastqFile fq("@r1\n"
                   "ACGT\n"
                   "+\n");

  EXPECT_THROW(reads::FastqDataset::fromUnpaired(fq.path()),
               reads::FastqDatasetError);
}

TEST(FastqDatasetPaired, CombinesMatchingPairs) {
  TempFastqFile r1("@read1\nACGT\n+\nIIII\n");

  TempFastqFile r2("@read1\nTTTT\n+\nJJJJ\n");

  auto ds = reads::FastqDataset::fromPairedFiles(r1.path(), r2.path());

  ASSERT_EQ(ds.size(), 1u);

  auto const &read = ds[0];

  ASSERT_EQ(read.numFragments(), 2u);

  expectFragment(read.fragments()[0], "ACGT", "IIII");
  expectFragment(read.fragments()[1], "TTTT", "JJJJ");
}

TEST(FastqDatasetPaired, RejectsIdMismatch) {
  TempFastqFile r1("@readA\nACGT\n+\nIIII\n");

  TempFastqFile r2("@readB\nTTTT\n+\nJJJJ\n");

  EXPECT_THROW(reads::FastqDataset::fromPairedFiles(r1.path(), r2.path(), true),
               reads::FastqDatasetError);
}

TEST(FastqDatasetPaired, RejectsDifferentRecordCounts) {
  TempFastqFile r1("@r1\nACGT\n+\nIIII\n"
                   "@r2\nACGT\n+\nIIII\n");

  TempFastqFile r2("@r1\nTTTT\n+\nJJJJ\n");

  EXPECT_THROW(reads::FastqDataset::fromPairedFiles(r1.path(), r2.path()),
               reads::FastqDatasetError);
}

TEST(FastqDatasetInterleaved, ParsesReadPairs) {
  TempFastqFile fq("@read1/1\nACGT\n+\nIIII\n"
                   "@read1/2\nTTTT\n+\nJJJJ\n");

  auto ds = reads::FastqDataset::fromInterleaved(fq.path());

  ASSERT_EQ(ds.size(), 1u);

  auto const &read = ds[0];

  ASSERT_EQ(read.numFragments(), 2u);

  expectFragment(read.fragments()[0], "ACGT", "IIII");
  expectFragment(read.fragments()[1], "TTTT", "JJJJ");
}

TEST(FastqDatasetInterleaved, RejectsOddRecordCount) {
  TempFastqFile fq("@r1\nACGT\n+\nIIII\n");

  EXPECT_THROW(reads::FastqDataset::fromInterleaved(fq.path()),
               reads::FastqDatasetError);
}

TEST(FastqDatasetInterleaved, RejectsMismatchedInterleavedIds) {
  TempFastqFile fq("@read1/1\nACGT\n+\nIIII\n"
                   "@read2/2\nTTTT\n+\nJJJJ\n");

  EXPECT_THROW(reads::FastqDataset::fromInterleaved(fq.path(), true),
               reads::FastqDatasetError);
}

TEST(FastqDatasetMetadata, CountsDiscardedBases) {
  TempFastqFile fq("@r1\n"
                   "ACGTNNNACGT\n"
                   "+\n"
                   "IIIIIIIIIII\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  EXPECT_EQ(ds.metadata().discarded_base_count, 3u);
}

TEST(FastqDatasetMetadata, RecordsInputAndLogicalReadCounts) {
  TempFastqFile fq("@r1\nACGT\n+\nIIII\n"
                   "@r2\nACGT\n+\nIIII\n");

  auto ds = reads::FastqDataset::fromUnpaired(fq.path());

  EXPECT_EQ(ds.metadata().input_record_count, 2u);
  EXPECT_EQ(ds.metadata().logical_read_count, 2u);
}
