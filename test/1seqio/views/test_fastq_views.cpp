#include <gtest/gtest.h>

namespace {

std::vector<std::string> collectSequences(auto const &range) {
  std::vector<std::string> result;

  for (auto const &fragment : range)
    result.push_back(to_string(fragment));

  return result;
}

} // namespace

TEST(FastqDatasetChunkView, ReadsRemainForwardOnlyInternally) {
  auto view = makeFastqView({"ACGTACGT"}, 1);

  auto chunk = view.range(0);

  auto seqs = collectSequences(chunk.forwardFragments(3));

  ASSERT_EQ(seqs.size(), 1u);
  EXPECT_EQ(seqs[0], "ACGTACGT");
}

TEST(FastqDatasetChunkView, ReverseFragmentsGeneratedOnDemand) {
  auto view = makeFastqView({"ACGTACGT"}, 1);

  auto chunk = view.range(0);

  auto forward = collectSequences(chunk.forwardFragments(3));

  auto reverse = collectSequences(chunk.reverseFragments(3));

  ASSERT_EQ(forward.size(), 1u);
  ASSERT_EQ(reverse.size(), 1u);

  EXPECT_EQ(reverse[0], revcomp(forward[0]));
}

TEST(FastqDatasetChunkView, AllReadsHaveColourZero) {
  auto view = makeFastqView({"AAAA", "CCCC", "GGGG"}, 1);

  auto chunk = view.range(0);

  for (auto const &fragment : chunk.forwardFragments(3)) {
    EXPECT_EQ(fragment.id(), 0u);
  }

  for (auto const &fragment : chunk.reverseFragments(3)) {
    EXPECT_EQ(fragment.id(), 0u);
  }
}

TEST(FastqDatasetChunkView, SplitsReadsOnN) {
  auto view = makeFastqView({"AAAANCCCC"}, 1);

  auto chunk = view.range(0);

  auto seqs = collectSequences(chunk.forwardFragments(3));

  ASSERT_EQ(seqs.size(), 2u);

  EXPECT_EQ(seqs[0], "AAAA");
  EXPECT_EQ(seqs[1], "CCCC");
}

TEST(FastqDatasetChunkView, ConsecutiveNsDoNotCreateEmptyReads) {
  auto view = makeFastqView({"AAAANNNNCCCC"}, 1);

  auto chunk = view.range(0);

  auto seqs = collectSequences(chunk.forwardFragments(3));

  ASSERT_EQ(seqs.size(), 2u);
}

TEST(FastqDatasetChunkView, ForwardAndReverseContainSameNumberOfFragments) {
  auto view = makeFastqView({"AAAA", "CCCC", "GGGG"}, 1);

  auto chunk = view.range(0);

  auto forward = collectSequences(chunk.forwardFragments(3));

  auto reverse = collectSequences(chunk.reverseFragments(3));

  EXPECT_EQ(forward.size(), reverse.size());
}

TEST(FastqDatasetChunkView, ForwardAndReverseTerminalsMatch) {
  auto view = makeFastqView({"AAAANCCCC"}, 1);

  auto chunk = view.range(0);

  auto forward = collectSequences(chunk.forwardTerminals());

  auto reverse = collectSequences(chunk.reverseTerminals());

  ASSERT_EQ(forward.size(), reverse.size());

  for (size_t i = 0; i < forward.size(); ++i) {
    EXPECT_EQ(reverse[i], revcomp(forward[i]));
  }
}

TEST(FastqDatasetChunkView, EmptyInputProducesNoFragments) {
  auto view = makeFastqView({}, 1);

  auto chunk = view.range(0);

  EXPECT_TRUE(collectSequences(chunk.forwardFragments(3)).empty());
}