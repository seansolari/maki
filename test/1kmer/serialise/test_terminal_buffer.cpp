
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/core/seq/concepts.hpp"
#include "maki/core/seq/io.hpp"
#include <gtest/gtest.h>

using seqan3::operator""_dna4;

TEST(TerminalAutoFitTests, TestMin) {
  TerminalBuffer buffer(1, 1, TerminalBuffer::autofit_tag);
  EXPECT_EQ(buffer.lengthBytes(), 1);
}

TEST(TerminalAutoFitTests, TestMax) {
  TerminalBuffer buffer(1, (uint8_t)-1, TerminalBuffer::autofit_tag);
  EXPECT_EQ(buffer.lengthBytes(), 1);
}

TEST(TerminalBufferIterator, ReadWriteFullKmer) {
  constexpr uint8_t k = 17;

  TerminalBuffer buf(1, k, TerminalBuffer::autofit_tag);
  auto it = buf.begin();

  Dna4Sequence seq{"ACGTACGTACGTACGTA"_dna4};

  it.writeKey(seq.end() - 1, seq.size());

  auto read = it.readKmer();

  EXPECT_EQ(read.size(), seq.size());
  EXPECT_EQ(read, seq);
}

TEST(TerminalBufferIterator, ReadWritePartialTerminal) {
  constexpr uint8_t k = 13;

  TerminalBuffer buf(1, k, TerminalBuffer::autofit_tag);
  auto it = buf.begin();

  Dna4Sequence seq{"ACGT"_dna4};

  it.writeKey(seq.end() - 1, seq.size());

  auto read = it.readKmer();

  EXPECT_EQ(read.size(), 4);
  EXPECT_EQ(read, seq);
}

TEST(TerminalBufferIterator, ReadWriteLongSuffix) {
  constexpr uint8_t k = 11;

  LongSuffix suffix(k);

  for (auto nt : "ACGTA"_dna4)
    suffix.push(seqan3::to_rank(nt));

  TerminalBuffer buf(1, k, TerminalBuffer::autofit_tag);

  auto it = buf.begin();
  it.writeKey(suffix);

  auto decoded = it.readTerminal();

  EXPECT_EQ(decoded.terminalLength(), suffix.terminalLength());
  EXPECT_EQ(decoded.view(), suffix.view());
}

TEST(TerminalBufferIterator, LengthDistinguishesEncodedPadding) {
  constexpr uint8_t k = 5;

  TerminalBuffer buf(2, k, TerminalBuffer::autofit_tag);

  auto it = buf.begin();

  LongSuffix a(k);
  a.push(0); // A

  it.writeKey(a);
  ++it;

  LongSuffix aa(k);
  aa.push(0);
  aa.push(0);

  it.writeKey(aa);

  auto first = buf.begin().readTerminal();
  auto second = buf.at(1).readTerminal();

  EXPECT_EQ(first.terminalLength(), 1);
  EXPECT_EQ(second.terminalLength(), 2);

  EXPECT_NE(first.terminalLength(), second.terminalLength());
}

TEST(TerminalBuffer, InsertRegionTerminalSequence) {
  constexpr uint8_t k = 5;

  auto seq = "ACGTAC"_dna4;

  TerminalBuffer buf(6, k, TerminalBuffer::autofit_tag);

  auto end = buf.insertRegion(buf.begin(), seq, false);

  EXPECT_EQ(std::distance(buf.begin(), end), 6);

  EXPECT_EQ(buf.at(0).readKmer(), Dna4Sequence{""_dna4});
  EXPECT_EQ(buf.at(1).readKmer(), Dna4Sequence{"A"_dna4});
  EXPECT_EQ(buf.at(2).readKmer(), Dna4Sequence{"AC"_dna4});
  EXPECT_EQ(buf.at(3).readKmer(), Dna4Sequence{"ACG"_dna4});
  EXPECT_EQ(buf.at(4).readKmer(), Dna4Sequence{"ACGT"_dna4});
  EXPECT_EQ(buf.at(5).readKmer(), Dna4Sequence{"ACGTA"_dna4});
}

TEST(TerminalBuffer, InsertRegionIncludesFullLengthKmer) {
  constexpr uint8_t k = 5;

  auto seq = "ACGTA"_dna4;

  TerminalBuffer buf(6, k, TerminalBuffer::autofit_tag);

  auto end = buf.insertRegion(buf.begin(), seq, true);
  --end;

  EXPECT_EQ(end.readSize(), 5);
  EXPECT_EQ(end.readKmer(), Dna4Sequence{"ACGTA"_dna4});
  EXPECT_EQ(end.readEdge(), TerminalBuffer::terminalEdge);
}

TEST(TerminalBuffer, InsertRegionWithSuffixMatch) {
  constexpr uint8_t k = 7;

  auto seq = "ACGTACGAT"_dna4;
  ShortSuffix suffix{"TA"_dna4};

  TerminalBuffer buf(8, k, k - 2, TerminalBuffer::autofit_tag);

  auto end = buf.insertRegion(buf.begin(), seq, false, suffix);

  EXPECT_EQ(std::distance(buf.begin(), end), 1);

  --end;
  EXPECT_EQ(end.readKmer(), Dna4Sequence{"ACG"_dna4});
}

TEST(TerminalBuffer, InsertRegionWithSuffixMiss) {
  constexpr uint8_t k = 7;

  auto seq = "ACGTACGT"_dna4;

  ShortSuffix suffix{"GG"_dna4};

  TerminalBuffer buf(8, k, TerminalBuffer::autofit_tag);

  auto end = buf.insertRegion(buf.begin(), seq, false, suffix);

  EXPECT_EQ(std::distance(buf.begin(), end), 0);
}

TEST(TerminalBuffer, SuffixComparisonMustRespectTerminalLength) {
  constexpr uint8_t k = 5;

  auto seq = "ACCCC"_dna4;

  ShortSuffix suffix{"AA"_dna4};

  TerminalBuffer buf(4, k, TerminalBuffer::autofit_tag);

  auto end = buf.insertRegion(buf.begin(), seq, true, suffix);

  EXPECT_EQ(std::distance(buf.begin(), end), 0);
}

TEST(TerminalBuffer, LockInsertDistinctTerminals) {
  constexpr uint8_t k = 5;

  LongSuffixGate lock(k);

  auto seq = "ACGTA"_dna4;

  TerminalBuffer buf(16, k, TerminalBuffer::autofit_tag);

  buf.insertRegion(lock, seq, true);

  EXPECT_GT(lock.count(), 0);
}

TEST(TerminalBuffer, LockInsertDeduplicates) {
  constexpr uint8_t k = 5;

  LongSuffixGate lock(k);

  auto seq = "ACGTA"_dna4;

  TerminalBuffer buf(32, k, TerminalBuffer::autofit_tag);

  buf.insertRegion(lock, seq, true);

  auto before = lock.count();

  buf.insertRegion(lock, seq, true);

  EXPECT_EQ(lock.count(), before);
}

TEST(TerminalBuffer, LockDistinguishesPaddingFromRealBases) {
  constexpr uint8_t k = 5;

  LongSuffixGate lock(k);

  TerminalBuffer buf(10, k, TerminalBuffer::autofit_tag);

  buf.insertRegion(lock, "AAAAA"_dna4, true);

  EXPECT_EQ(lock.count(), 6);
}

TEST(TerminalBuffer, InsertFragmentTerminalsExactK) {
  constexpr uint8_t k = 5;

  TerminalBuffer buf(16, k, TerminalBuffer::autofit_tag);

  auto end = buf.insertFragmentTerminals(
      buf.begin(), SequenceFragment("ACGTA"_dna4, 0, true));

  ASSERT_EQ(std::distance(buf.begin(), end), 5);

  EXPECT_EQ(buf.at(0).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{""_dna4});

  EXPECT_EQ(buf.at(1).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"A"_dna4});

  EXPECT_EQ(buf.at(2).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"AC"_dna4});

  EXPECT_EQ(buf.at(3).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"ACG"_dna4});

  EXPECT_EQ(buf.at(4).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"ACGT"_dna4});
}

TEST(TerminalBuffer, InsertFragmentTerminalsIgnoreCompleteKmers) {
  constexpr uint8_t k = 5;

  TerminalBuffer buf(16, k, TerminalBuffer::autofit_tag);

  auto end = buf.insertFragmentTerminals(
      buf.begin(), SequenceFragment("ACGTACGT"_dna4, 0, true));

  ASSERT_EQ(std::distance(buf.begin(), end), 5);

  EXPECT_EQ(buf.at(0).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{""_dna4});

  EXPECT_EQ(buf.at(1).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"A"_dna4});

  EXPECT_EQ(buf.at(2).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"AC"_dna4});

  EXPECT_EQ(buf.at(3).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"ACG"_dna4});

  EXPECT_EQ(buf.at(4).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"ACGT"_dna4});
}

TEST(TerminalBuffer, InsertFragmentKmersExactK) {
  constexpr uint8_t k = 5;

  TerminalBuffer buf(16, k, k - 2, TerminalBuffer::autofit_tag);

  auto end = buf.insertFragmentKmers(buf.begin(),
                                     SequenceFragment("ACGTA"_dna4, 0, true),
                                     ShortSuffix{"AC"_dna4});

  ASSERT_EQ(std::distance(buf.begin(), end), 1);

  EXPECT_EQ(buf.at(0).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{""_dna4});

  EXPECT_EQ(buf.at(0).readEdge(), 0b10);
}

TEST(TerminalBuffer, InsertFragmentKmersCompleteSequence) {
  constexpr uint8_t k = 5;

  TerminalBuffer buf(32, k, k - 2, TerminalBuffer::autofit_tag);

  auto end = buf.insertFragmentKmers(buf.begin(),
                                     SequenceFragment("ACGTACGT"_dna4, 0, true),
                                     ShortSuffix{"AC"_dna4});

  ASSERT_EQ(std::distance(buf.begin(), end), 2);

  EXPECT_EQ(buf.at(0).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{""_dna4});

  EXPECT_EQ(buf.at(0).readEdge(), 0b10);

  EXPECT_EQ(buf.at(1).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"CGT"_dna4});

  EXPECT_EQ(buf.at(0).readEdge(), 0b10);
}

TEST(TerminalBuffer, InsertFragmentKmersPreservesTerminalLengths) {
  constexpr uint8_t k = 5;

  TerminalBuffer buf(16, k, k-1, TerminalBuffer::autofit_tag);

  auto end = buf.insertFragmentKmers(
      buf.begin(), SequenceFragment("AAAAAAAA"_dna4, 0, false),
      ShortSuffix{"A"_dna4});

  ASSERT_EQ(std::distance(buf.begin(), end), 7);

  EXPECT_EQ(buf.at(0).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{""_dna4});

  EXPECT_EQ(buf.at(1).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"A"_dna4});

  EXPECT_EQ(buf.at(2).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"AA"_dna4});

  EXPECT_EQ(buf.at(3).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"AAA"_dna4});

  EXPECT_EQ(buf.at(4).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"AAAA"_dna4});

  EXPECT_EQ(buf.at(5).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"AAAA"_dna4});

  EXPECT_EQ(buf.at(6).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"AAAA"_dna4});
}
