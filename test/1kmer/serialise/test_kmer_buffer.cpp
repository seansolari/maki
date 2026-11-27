
#include "maki/build/kmers/buffers/kmers.hpp"
#include "maki/core/seq/io.hpp"
#include "seqan3/alphabet/nucleotide/dna4.hpp"
#include <cstdint>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using seqan3::operator""_dna4;

static std::string kmerString(auto const &it, uint8_t k) {
  return toString(it.readKmer(k));
}

TEST(KmerEncodingTests, WriteKmerBits) {
  Dna4Sequence original{"ACGTACGTAC"_dna4};

  std::array<uint8_t, 3> buffer{};
  writeMerTo(original.begin(), buffer.data(), 10);

  auto decoded = readMerFrom(buffer.data(), 10);

  EXPECT_EQ(original, decoded);
}

TEST(KmerBufferInsertFragment, NonTerminalFragment) {
  constexpr uint8_t k = 3;

  KmerBuffer buffer(10, sizeof(uint32_t), k, k);

  std::vector<seqan3::dna4> seq = "ATACGT"_dna4;

  SequenceFragment fragment(std::move(seq), 17, false);

  auto it = buffer.begin();

  it = buffer.insertFragment(it, fragment);

  auto r = buffer.begin();

  EXPECT_EQ(kmerString(r, k), "ATA");
  EXPECT_EQ(r.readEdge(), seqan3::dna4{}.assign_char('C').to_rank());
  EXPECT_EQ(r.readValue().colour(), 17u);

  ++r;

  EXPECT_EQ(kmerString(r, k), "TAC");
  EXPECT_EQ(r.readEdge(), seqan3::dna4{}.assign_char('G').to_rank());
  EXPECT_EQ(r.readValue().colour(), 17u);

  ++r;

  EXPECT_EQ(kmerString(r, k), "ACG");
  EXPECT_EQ(r.readEdge(), seqan3::dna4{}.assign_char('T').to_rank());
  EXPECT_EQ(r.readValue().colour(), 17u);
}

TEST(KmerBufferInsertFragment, TerminalFragmentAddsTerminalEdge) {
  constexpr uint8_t k = 3;

  KmerBuffer buffer(10, sizeof(uint32_t), k, k);

  std::vector<seqan3::dna4> seq = "ATACGT"_dna4;

  SequenceFragment fragment(std::move(seq), 5, true);

  buffer.insertFragment(buffer.begin(), fragment);

  auto r = buffer.begin();

  EXPECT_EQ(kmerString(r, k), "ATA");
  ++r;

  EXPECT_EQ(kmerString(r, k), "TAC");
  ++r;

  EXPECT_EQ(kmerString(r, k), "ACG");
  ++r;

  EXPECT_EQ(kmerString(r, k), "CGT");

  EXPECT_EQ(r.readEdge(), KmerBuffer::terminalEdge);
  EXPECT_EQ(r.readValue().colour(), 0u);
}

TEST(KmerBufferInsertFragment, PreservesInputOrder) {
  constexpr uint8_t k = 3;

  KmerBuffer buffer(16, sizeof(uint32_t), k, k);

  std::vector<seqan3::dna4> seq = "AAAACCCC"_dna4;

  SequenceFragment fragment(std::move(seq), 1, false);

  buffer.insertFragment(buffer.begin(), fragment);

  std::vector<std::string> expected{"AAA", "AAA", "AAC", "ACC", "CCC"};

  auto r = buffer.begin();

  for (auto const &kmer : expected) {
    EXPECT_EQ(kmerString(r, k), kmer);
    ++r;
  }
}

TEST(KmerBufferInsertFragment, SuffixFilterSingleMatch) {
  constexpr uint8_t k = 3;

  KmerBuffer buffer(10, sizeof(uint32_t), k, k - 1);

  std::vector<seqan3::dna4> seq = "ATACGT"_dna4;

  SequenceFragment fragment(std::move(seq), 11, false);

  ShortSuffix suffix("G"_dna4);

  buffer.insertFragment(buffer.begin(), fragment, suffix);

  auto r = buffer.begin();

  EXPECT_EQ(kmerString(r, k - 1), "AC");

  EXPECT_EQ(r.readEdge(), seqan3::dna4{}.assign_char('T').to_rank());

  EXPECT_EQ(r.readValue().colour(), 11u);
}

TEST(KmerBufferInsertFragment, SuffixFilterMultipleMatches) {
  constexpr uint8_t k = 3;

  KmerBuffer buffer(10, sizeof(uint32_t), k, 1);

  std::vector<seqan3::dna4> seq = "ACGACGAC"_dna4;

  SequenceFragment fragment(std::move(seq), 9, false);

  ShortSuffix suffix("CG"_dna4);

  buffer.insertFragment(buffer.begin(), fragment, suffix);

  auto r = buffer.begin();

  EXPECT_EQ(kmerString(r, 1), "A");
  EXPECT_EQ(r.readEdge(), seqan3::dna4{}.assign_char('A').to_rank());

  ++r;

  EXPECT_EQ(kmerString(r, 1), "A");
  EXPECT_EQ(r.readEdge(), seqan3::dna4{}.assign_char('A').to_rank());
}

TEST(KmerBufferInsertFragment, SuffixFilterPreservesTerminalEdge) {
  constexpr uint8_t k = 3;

  KmerBuffer buffer(10, sizeof(uint32_t), k, 1);

  std::vector<seqan3::dna4> seq = "ACGT"_dna4;

  SequenceFragment fragment(std::move(seq), 3, true);

  ShortSuffix suffix("GT"_dna4);

  buffer.insertFragment(buffer.begin(), fragment, suffix);

  auto r = buffer.begin();

  EXPECT_EQ(kmerString(r, 1), "C");
  EXPECT_EQ(r.readEdge(), TerminalEdge);
  EXPECT_EQ(r.readValue().colour(), 0u);
}

TEST(KmerBufferInsertFragment, SuffixFilterAcrossByteBoundary) {
  constexpr uint8_t k = 9;

  KmerBuffer buffer(10, sizeof(uint32_t), k, 6);

  std::vector<seqan3::dna4> seq = "ACGTACGTAC"_dna4;

  SequenceFragment fragment(std::move(seq), 42, false);

  ShortSuffix suffix("GTA"_dna4);

  buffer.insertFragment(buffer.begin(), fragment, suffix);

  auto r = buffer.begin();

  EXPECT_EQ(kmerString(r, 6), "ACGTAC");

  EXPECT_EQ(r.readEdge(), seqan3::dna4{}.assign_char('C').to_rank());

  EXPECT_EQ(r.readValue().colour(), 42u);
}
