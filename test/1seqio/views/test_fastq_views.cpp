#include "maki/build/io/fastq.hpp"
#include "maki/core/utils/tempfile.hpp"
#include "view_common.hpp"
#include <gtest/gtest.h>

#include <sstream>
#include <string>
#include <vector>

inline std::string makeFastq(std::vector<std::string> const &sequences) {
  std::ostringstream out;

  for (std::size_t i = 0; i < sequences.size(); ++i) {
    auto const &seq = sequences[i];

    out << "@read_" << i << '\n';
    out << seq << '\n';
    out << "+\n";
    out << std::string(seq.size(), 'I') << '\n';
  }

  return out.str();
}

class FastqDatasetChunkViewTests : public testing::Test {
protected:
  inline reads::FastqDatasetChunkView
  makeFastqView(std::vector<std::string> const &sequences, std::size_t chunk) {
    tempio::TempFastqFile fq(makeFastq(sequences));
    ds = reads::FastqDataset::fromUnpaired(fq.path());
    return reads::FastqDatasetChunkView(ds, chunk);
  }

  reads::FastqDataset ds;
};

TEST_F(FastqDatasetChunkViewTests, ReadsRemainForwardOnlyInternally) {
  auto view = makeFastqView({"ACGTACGT"}, 1);

  auto chunk = view.range(0);

  auto seqs = collectSequences(chunk.forwardFragments(3));

  ASSERT_EQ(seqs.size(), 1u);
  EXPECT_EQ(seqs[0], "ACGTACGT");
}

TEST_F(FastqDatasetChunkViewTests, ReverseFragmentsGeneratedOnDemand) {
  auto view = makeFastqView({"ACGTACGTA"}, 1);

  auto chunk = view.range(0);

  auto forward = collectSequences(chunk.forwardFragments(3));
  auto reverse = collectSequences(chunk.reverseFragments(3));

  ASSERT_EQ(forward.size(), 1u);
  ASSERT_EQ(reverse.size(), 1u);

  EXPECT_EQ(reverse[0], revcomp(forward[0]));
}

TEST_F(FastqDatasetChunkViewTests, AllReadsHaveColourZero) {
  auto view = makeFastqView({"AAAA", "CCCC", "GGGG"}, 1);

  auto chunk = view.range(0);

  for (auto const &fragment : chunk.forwardFragments(3)) {
    EXPECT_EQ(fragment.id(), 0u);
  }

  for (auto const &fragment : chunk.reverseFragments(3)) {
    EXPECT_EQ(fragment.id(), 0u);
  }
}

TEST_F(FastqDatasetChunkViewTests, SplitsReadsOnN) {
  auto view = makeFastqView({"AAAANCCCC"}, 1);

  auto chunk = view.range(0);

  auto seqs = collectSequences(chunk.forwardFragments(3));

  ASSERT_EQ(seqs.size(), 2u);

  EXPECT_EQ(seqs[0], "AAAA");
  EXPECT_EQ(seqs[1], "CCCC");
}

TEST_F(FastqDatasetChunkViewTests, ConsecutiveNsDoNotCreateEmptyReads) {
  auto view = makeFastqView({"AAAANNNNCCCC"}, 1);

  auto chunk = view.range(0);

  auto seqs = collectSequences(chunk.forwardFragments(3));

  ASSERT_EQ(seqs.size(), 2u);
}

TEST_F(FastqDatasetChunkViewTests,
       ForwardAndReverseContainSameNumberOfFragments) {
  auto view = makeFastqView({"AAAA", "CCCC", "GGGG"}, 1);

  auto chunk = view.range(0);

  auto forward = collectSequences(chunk.forwardFragments(3));

  auto reverse = collectSequences(chunk.reverseFragments(3));

  EXPECT_EQ(forward.size(), reverse.size());
}

TEST_F(FastqDatasetChunkViewTests, ForwardAndReverseTerminalsMatch) {
  auto view = makeFastqView({"AAAANCCCC"}, 1);

  auto chunk = view.range(0);

  auto forward = collectSequences(chunk.forwardTerminals());

  auto reverse = collectSequences(chunk.reverseTerminals());

  ASSERT_EQ(forward.size(), reverse.size());

  for (size_t i = 0; i < forward.size(); ++i) {
    EXPECT_EQ(reverse[i], revcomp(forward[i]));
  }
}

TEST_F(FastqDatasetChunkViewTests, EmptyInputProducesNoFragments) {
  auto view = makeFastqView({}, 1);

  auto chunk = view.range(0);

  EXPECT_TRUE(collectSequences(chunk.forwardFragments(3)).empty());
}