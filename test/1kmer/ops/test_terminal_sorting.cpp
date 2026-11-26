
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "seqan3/alphabet/nucleotide/dna4.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using ::testing::ContainerEq;

using seqan3::operator""_dna4;

namespace {

LongSuffix makeTerminal(uint8_t k, std::string_view seq) {
  LongSuffix t(k);

  for (char c : seq) {
    t.push(seqan3::to_rank(seqan3::dna4{}.assign_char(c)));
  }

  return t;
}

void writeRecord(TerminalBuffer::iterator it, LongSuffix const &t,
                 unsigned char edge_) {
  it.writeKey(t);
  it.writeEdge(edge_);
}

std::vector<std::string> exportKmers(const TerminalBuffer &buffer) {
  std::vector<std::string> expected;

  for (auto it = buffer.begin(), end = buffer.end(); it != end; ++it) {
    expected.push_back(toString(it.readKmer()));
  }

  return expected;
}

std::vector<unsigned char> exportEdges(const TerminalBuffer &buffer) {
  std::vector<unsigned char> expected;

  for (auto it = buffer.begin(), end = buffer.end(); it != end; ++it) {
    expected.push_back(it.readEdge());
  }

  return expected;
}

} // namespace

TEST(TerminalBufferSort, SortTerminalsOnly) {
  constexpr uint8_t k = 5;

  TerminalBuffer buf(4, k, TerminalBuffer::autofit_tag);

  writeRecord(buf.at(0), makeTerminal(k, "AC"), 0);
  writeRecord(buf.at(1), makeTerminal(k, "A"), 0);
  writeRecord(buf.at(2), makeTerminal(k, "G"), 0);
  writeRecord(buf.at(3), makeTerminal(k, "AA"), 0);

  buf.sort();

  EXPECT_EQ(buf.at(0).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"A"_dna4});

  EXPECT_EQ(buf.at(1).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"AA"_dna4});

  EXPECT_EQ(buf.at(2).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"AC"_dna4});

  EXPECT_EQ(buf.at(3).readKmer(),
            seqan3::bitpacked_sequence<seqan3::dna4>{"G"_dna4});
}

TEST(TerminalBufferSort, DistinguishTerminalLengths) {
  constexpr uint8_t k = 5;

  TerminalBuffer buf(2, k, TerminalBuffer::autofit_tag);

  writeRecord(buf.at(0), makeTerminal(k, "AA"), 0);
  writeRecord(buf.at(1), makeTerminal(k, "A"), 0);

  buf.sort();

  auto a = buf.at(0).readTerminal();
  auto aa = buf.at(1).readTerminal();

  EXPECT_EQ(a.terminalLength(), 1);
  EXPECT_EQ(aa.terminalLength(), 2);
}

TEST(TerminalBufferSort, SortMixedTerminalAndKmerRecords) {
  constexpr uint8_t k = 5;

  TerminalBuffer buf(6, k, TerminalBuffer::autofit_tag);

  writeRecord(buf.at(0), makeTerminal(k, "CGTAC"), 0);
  writeRecord(buf.at(1), makeTerminal(k, "AAAAA"), 0);
  writeRecord(buf.at(2), makeTerminal(k, "A"), 0);
  writeRecord(buf.at(3), makeTerminal(k, "ACGT"), 0);
  writeRecord(buf.at(4), makeTerminal(k, "AC"), 0);
  writeRecord(buf.at(5), makeTerminal(k, "ACGTA"), 0);

  buf.sort();

  std::vector<std::string> observed = exportKmers(buf),
                           expected = {"A",  "AAAAA", "ACGTA",
                                       "AC", "CGTAC", "ACGT"};

  EXPECT_THAT(observed, ContainerEq(expected));
}

TEST(TerminalBufferSort, TemporaryBufferSortMatchesInPlaceSort) {
  constexpr uint8_t k = 5;

  TerminalBuffer a(8, k, TerminalBuffer::autofit_tag),
      b(8, k, TerminalBuffer::autofit_tag),
      tmp(8, k, TerminalBuffer::autofit_tag);

  for (size_t i = 0; i < 8; ++i) {
    auto t = makeTerminal(k, std::string(i + 1, 'A'));

    writeRecord(a.at(i), t, 0);
    writeRecord(b.at(i), t, 0);
  }

  a.sort();
  b.sort(&tmp);

  for (size_t i = 0; i < 8; ++i) {
    EXPECT_EQ(a.at(i).readSize(), b.at(i).readSize());

    EXPECT_EQ(a.at(i).readKmer(), b.at(i).readKmer());
  }
}

TEST(TerminalBufferSort, DifferentEdgesNotIdentical) {
  constexpr uint8_t k = 5;

  TerminalBuffer buf(2, k, TerminalBuffer::autofit_tag);

  auto t = makeTerminal(k, "ACGTA");

  writeRecord(buf.at(0), t, 1);
  writeRecord(buf.at(1), t, 2);

  auto sorted = buf.OOPsort(true);

  EXPECT_EQ(sorted.size(), 2);

  EXPECT_NE(sorted.at(0).readEdge(), sorted.at(1).readEdge());
}

TEST(TerminalBufferSort, UniqueDoesNotMergeDifferentEdges) {
  constexpr uint8_t k = 5;

  TerminalBuffer buf(3, k, TerminalBuffer::autofit_tag);

  auto t = makeTerminal(k, "ACGTA");

  writeRecord(buf.at(0), t, 1);
  writeRecord(buf.at(1), t, 1);
  writeRecord(buf.at(2), t, 2);

  buf.unique();

  EXPECT_EQ(buf.size(), 2);

  EXPECT_EQ(buf.at(0).readEdge(), 1);
  EXPECT_EQ(buf.at(1).readEdge(), 2);
}

TEST(TerminalBufferSort, OOPsortRemovesTrueDuplicates) {
  constexpr uint8_t k = 5;

  TerminalBuffer buf(3, k, TerminalBuffer::autofit_tag);

  auto t = makeTerminal(k, "ACGTA");

  writeRecord(buf.at(0), t, 2);
  writeRecord(buf.at(1), t, 2);
  writeRecord(buf.at(2), t, 2);

  auto sorted = buf.OOPsort(true);

  EXPECT_EQ(sorted.size(), 1);
}

TEST(TerminalBufferSort, OOPsortDistinguishesPaddingFromBases) {
  constexpr uint8_t k = 5;

  TerminalBuffer buf(2, k, TerminalBuffer::autofit_tag);

  writeRecord(buf.at(0), makeTerminal(k, "A"), 0);
  writeRecord(buf.at(1), makeTerminal(k, "AA"), 0);

  auto sorted = buf.OOPsort(true);

  EXPECT_EQ(sorted.size(), 2);

  EXPECT_EQ(sorted.at(0).readTerminal().terminalLength(), 1);
  EXPECT_EQ(sorted.at(1).readTerminal().terminalLength(), 2);
}

TEST(TerminalBufferSort, MixedTerminalsKmersColoursAndDuplicates) {
  constexpr uint8_t k = 5;

  TerminalBuffer buf(8, k, TerminalBuffer::autofit_tag);

  writeRecord(buf.at(0), makeTerminal(k, "A"), 0);
  writeRecord(buf.at(1), makeTerminal(k, "AA"), 0);
  writeRecord(buf.at(2), makeTerminal(k, "ACGTA"), 1);
  writeRecord(buf.at(3), makeTerminal(k, "ACGTA"), 2);
  writeRecord(buf.at(4), makeTerminal(k, "CGTAC"), 1);
  writeRecord(buf.at(5), makeTerminal(k, "CGTAC"), 1);
  writeRecord(buf.at(6), makeTerminal(k, "AC"), 0);
  writeRecord(buf.at(7), makeTerminal(k, "ACGT"), 0);

  auto sorted = buf.OOPsort(true);

  EXPECT_EQ(sorted.size(), 7);

  std::vector<std::string> observedKmers = exportKmers(sorted),
                           expectedKmers = {"A",  "AA",    "ACGTA", "ACGTA",
                                            "AC", "CGTAC", "ACGT"};

  EXPECT_THAT(observedKmers, ContainerEq(expectedKmers));

  std::vector<unsigned char> observedEdges = exportEdges(sorted),
                             expectedEdges = {0, 0, 1, 2, 0, 1, 0};

  EXPECT_THAT(observedEdges, ContainerEq(expectedEdges));
}
