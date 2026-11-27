
#include "maki/build/kmers/buffers/base_buffer.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "test_common.hpp"

class SmallTerminalsLessThan : public testing::Test {
protected:
  SmallTerminalsLessThan()
      : k(3),
        buffer(6, k, TerminalBuffer::autofit_tag)
  {
    writeRecord(buffer.at(0), makeTerminal(k, ""), 0b10);
    writeRecord(buffer.at(1), makeTerminal(k, "A"), 0b10);
    writeRecord(buffer.at(2), makeTerminal(k, "C"), 0b10);
    writeRecord(buffer.at(3), makeTerminal(k, "TT"), 0b11);
    writeRecord(buffer.at(4), makeTerminal(k, "TT"), 0b11);
    writeRecord(buffer.at(5), makeTerminal(k, "T"), TerminalEdge);
  }

  uint8_t k;
  TerminalBuffer buffer;
};

TEST_F(SmallTerminalsLessThan, TerminalsLessThan_Small_lt_size) {
  TerminalsLessThan Cmp{buffer};
  ASSERT_TRUE(Cmp(0, 1));
}

TEST_F(SmallTerminalsLessThan, TerminalsLessThan_Small_lt_kmer) {
  TerminalsLessThan Cmp{buffer};
  ASSERT_TRUE(Cmp(1, 2));
}

TEST_F(SmallTerminalsLessThan, TerminalsLessThan_Small_lt_size_and_kmer) {
  TerminalsLessThan Cmp{buffer};
  ASSERT_TRUE(Cmp(2, 3));
}

TEST_F(SmallTerminalsLessThan, TerminalsLessThan_Small_eq) {
  TerminalsLessThan Cmp{buffer};
  ASSERT_FALSE(Cmp(3, 4));
}

TEST_F(SmallTerminalsLessThan, TerminalsLessThan_Small_gt) {
  TerminalsLessThan Cmp{buffer};
  ASSERT_FALSE(Cmp(4, 5));
}

class BigTerminalsLessThan : public testing::Test {
protected:
  BigTerminalsLessThan()
      : k(9),
        buffer(10, k, TerminalBuffer::autofit_tag)
  {
    writeRecord(buffer.at(0), makeTerminal(k, "AAAA"), 0b010);
    writeRecord(buffer.at(1), makeTerminal(k, "AAAAAAAA"), 0b010);
    writeRecord(buffer.at(2), makeTerminal(k, "AAAAAAAA"), 0b011);
    writeRecord(buffer.at(3), makeTerminal(k, "TACGAGGT"), 0b010);
    writeRecord(buffer.at(4), makeTerminal(k, "GGGT"), 0b010);
    writeRecord(buffer.at(5), makeTerminal(k, "GGGT"), 0b010);
    writeRecord(buffer.at(6), makeTerminal(k, "GGTT"), 0b010);
    writeRecord(buffer.at(7), makeTerminal(k, "AAAAGGTT"), 0b010);
    writeRecord(buffer.at(8), makeTerminal(k, "AAAAGGTT"), 0b010);
    writeRecord(buffer.at(9), makeTerminal(k, "AAAAGGTT"), 0b110);
  }

  uint8_t k;
  TerminalBuffer buffer;
};

TEST_F(BigTerminalsLessThan, TerminalsLessThan_Large_lt_size) {
  TerminalsLessThan Cmp{buffer};
  ASSERT_TRUE(Cmp(0, 1));
}

TEST_F(BigTerminalsLessThan, TerminalsLessThan_Large_lt_edge) {
  TerminalsLessThan Cmp{buffer};
  ASSERT_TRUE(Cmp(1, 2));
}

TEST_F(BigTerminalsLessThan, TerminalsLessThan_Large_lt_kmer) {
  TerminalsLessThan Cmp{buffer};
  ASSERT_TRUE(Cmp(3, 4));
}

TEST_F(BigTerminalsLessThan, TerminalsLessThan_Large_eq) {
  TerminalsLessThan Cmp{buffer};
  ASSERT_FALSE(Cmp(4, 5));
}

TEST_F(BigTerminalsLessThan, TerminalsLessThan_Large_gt_kmer) {
  TerminalsLessThan Cmp{buffer};
  ASSERT_FALSE(Cmp(6, 5));
}

TEST_F(BigTerminalsLessThan, TerminalsLessThan_Large_gt_size) {
  TerminalsLessThan Cmp{buffer};
  ASSERT_FALSE(Cmp(7, 6));
}
