
#include "maki/build/io/fasta.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/build/kmers/construct_terminals.hpp"
#include "test_common.hpp"
#include <cstdint>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>

using ::testing::ElementsAreArray;

using namespace seqan3::literals;

// custom matchers
// -----------------------------------------------------------------------------

void CheckMatrixEq(TerminalBuffer const &res, ndim::Matrix<> const &exp) {
  auto __lfirst = res.constBegin(), __llast = res.constEnd();
  auto __rfirst = exp.cbegin(), __rlast = exp.cend();

  for (; (__lfirst != __llast) && (__rfirst != __rlast);
       (void)++__lfirst, (void)++__rfirst)
    ASSERT_TRUE(*__lfirst == *__rfirst);

  ASSERT_EQ(__lfirst, __llast);
  ASSERT_EQ(__rfirst, __rlast);
}

// rolling suffix
// ------------------------------------------------------------------

class RollingTerminalTests : public testing::Test {
protected:
  RollingTerminalTests() : _sfx(7) {}

  LongSuffix _sfx;
};

TEST_F(RollingTerminalTests, Initial) {
  EXPECT_EQ(_sfx.terminalLength(), 0);
  ASSERT_THAT(_sfx.view(), ElementsAreArray({0, 0}));
}

TEST_F(RollingTerminalTests, RollBackNotFull) {
  _sfx.push(0b11);
  _sfx.push(0b10);
  _sfx.push(0b01);
  _sfx.push(0b00);
  _sfx.push(0b11);
  _sfx.push(0b10);
  EXPECT_EQ(_sfx.terminalLength(), 6);
  ASSERT_THAT(_sfx.view(), ElementsAreArray({0b01101100, 0b00101100}));
}

TEST_F(RollingTerminalTests, RollBackFull) {
  for (std::size_t i = 0; i < 7; ++i)
    _sfx.push(0b10);
  EXPECT_EQ(_sfx.terminalLength(), 7);
  ASSERT_THAT(_sfx.view(), ElementsAreArray({0b10101010, 0b00101010}));
}

TEST_F(RollingTerminalTests, RollBackExceed) {
  for (std::size_t i = 0; i < 14; ++i)
    _sfx.push(0b10);
  EXPECT_EQ(_sfx.terminalLength(), 7);
  ASSERT_THAT(_sfx.view(), ElementsAreArray({0b10101010, 0b00101010}));
}

// initialisation tests
// --------------------------------------------------------------------

// serialise

class TerminalAutoFitTests : public testing::Test {};

TEST_F(TerminalAutoFitTests, TestMin) {
  TerminalBuffer buffer(1, 1, TerminalBuffer::autofit_tag);
  EXPECT_EQ(buffer.lengthBytes(), 1);
}

TEST_F(TerminalAutoFitTests, TestMax) {
  TerminalBuffer buffer(1, (uint8_t)-1, TerminalBuffer::autofit_tag);
  EXPECT_EQ(buffer.lengthBytes(), 1);
}

class TerminalBufferSerialiseTests : public testing::Test {
protected:
  TerminalBufferSerialiseTests()
      : k(8), keff(7), buffer(1,    // num entries
                              k,    // k-mer size
                              keff, // effective k-mer size
                              3)    // length bytes
  {}

  uint8_t k, keff;
  TerminalBuffer buffer;
};

TEST_F(TerminalBufferSerialiseTests, MemCpySizeT) {
  std::size_t index = 3519u;

  auto it = buffer.at(0);
  it.writeSize(index);

  std::size_t recovered_index = it.readSize();
  ASSERT_EQ(index, recovered_index);
}

TEST_F(TerminalBufferSerialiseTests, MemCpyValue) {
  uint8_t e_ = 0b01u;

  auto it = buffer.at(0);
  it.writeEdge(e_);

  uint8_t recoveredValue = it.readEdge();
  ASSERT_EQ(e_, recoveredValue);
}

TEST_F(TerminalBufferSerialiseTests, WriteEmptyKey) {
  // create key

  LongSuffix _sfx(keff);

  // write key

  buffer.at(0).writeKey(_sfx);

  // check key value

  uint8_t expData[] = {0, 0, 0, 0, 0};
  ndim::Span<uint8_t *> exp{&expData[0], 5};

  EXPECT_EQ(buffer[0].first(buffer.keyBytes()), exp);
  EXPECT_EQ(buffer.at(0).readSize(), 0);
  EXPECT_EQ(buffer.at(0).keyMSB(), 0);
}

TEST_F(TerminalBufferSerialiseTests, WritePartialKey) {
  // create key

  LongSuffix _sfx(keff);
  _sfx.push(0b11);
  _sfx.push(0b10);
  _sfx.push(0b01);
  _sfx.push(0b00);
  _sfx.push(0b11);
  _sfx.push(0b10);

  // write key

  buffer.at(0).writeKey(_sfx);

  // check key value

  uint8_t expData[] = {6, 0, 0, 0b01101100, 0b00101100};
  ndim::Span<uint8_t *> exp{&expData[0], 5};

  EXPECT_EQ(buffer[0].first(buffer.keyBytes()), exp);
  EXPECT_EQ(buffer.at(0).readSize(), 6);
  EXPECT_EQ(buffer.at(0).keyMSB(), 0b00101100);
}

TEST_F(TerminalBufferSerialiseTests, WriteFullKey) {
  // create key

  LongSuffix _sfx(keff);
  for (std::size_t i = 0; i < 7; ++i)
    _sfx.push(0b10);

  // write key

  buffer.at(0).writeKey(_sfx);

  // check key value

  uint8_t expData[] = {7, 0, 0, 0b10101010, 0b00101010};
  ndim::Span<uint8_t *> exp{&expData[0], 5};

  EXPECT_EQ(buffer[0].first(buffer.keyBytes()), exp);
  EXPECT_EQ(buffer.at(0).readSize(), 7);
  EXPECT_EQ(buffer.at(0).keyMSB(), 0b00101010);
}

TEST_F(TerminalBufferSerialiseTests, WriteOverflowKey) {
  // create key

  LongSuffix _sfx(keff);
  for (std::size_t i = 0; i < 14; ++i)
    _sfx.push(0b10);

  // write key

  buffer.at(0).writeKey(_sfx);

  // check key value

  uint8_t expData[] = {7, 0, 0, 0b10101010, 0b00101010};
  ndim::Span<uint8_t *> exp{&expData[0], 5};

  EXPECT_EQ(buffer[0].first(buffer.keyBytes()), exp);
  EXPECT_EQ(buffer.at(0).readSize(), 7);
  EXPECT_EQ(buffer.at(0).keyMSB(), 0b00101010);
}

TEST_F(TerminalBufferSerialiseTests, ReversePartialKey) {
  // create key

  Dna4Sequence data{"ACGT"_dna4};

  // write key

  buffer.at(0).writeKey(data.cbegin() + (data.size() - 1), data.size());

  // check key value

  uint8_t expData[] = {4, 0, 0, 0b00000000, 0b00111001};
  ndim::Span<uint8_t *> exp{&expData[0], 5};

  EXPECT_EQ(buffer[0].first(buffer.keyBytes()), exp);
  EXPECT_EQ(buffer.at(0).readSize(), 4);
  EXPECT_EQ(buffer.at(0).keyMSB(), 0b00111001);
}

TEST_F(TerminalBufferSerialiseTests, ReverseFullKey) {
  // create key

  Dna4Sequence data{"ACGTACG"_dna4};

  // write key

  buffer.at(0).writeKey(data.cbegin() + (data.size() - 1), data.size());

  // check key value

  uint8_t expData[] = {7, 0, 0, 0b11100100, 0b00100100};
  ndim::Span<uint8_t *> exp{&expData[0], 5};

  EXPECT_EQ(buffer[0].first(buffer.keyBytes()), exp);
  EXPECT_EQ(buffer.at(0).readSize(), 7);
  EXPECT_EQ(buffer.at(0).keyMSB(), 0b00100100);
}

// comparison operators

class TerminalsLessThanTests : public testing::Test {
protected:
  TerminalsLessThanTests()
      : ksmall(3), kbig(9), num_genomes(2),
        small_buffer(ksmall,
                     {{/* size */ 0b00000000, /* kmer */ 0b00000000,
                       /* edge */ 0b00000010},
                      {/* size */ 0b00000001, /* kmer */ 0b00000000,
                       /* edge */ 0b00000010},
                      {/* size */ 0b00000001, /* kmer */ 0b00010000,
                       /* edge */ 0b00000010},
                      {/* size */ 0b00000010, /* kmer */ 0b00111100,
                       /* edge */ 0b00000011},
                      {/* size */ 0b00000010, /* kmer */ 0b00111100,
                       /* edge */ 0b00000011},
                      {/* size */ 0b00000001, /* kmer */ 0b00110000,
                       /* edge */ 0b00000111}},
                     TerminalBuffer::autofit_tag),
        big_buffer(kbig,
                   {{/* size */ 0b00000100, /* kmer */ 0b00000000, 0b00000000,
                     0b00000000, /* edge */ 0b00000010},
                    {/* size */ 0b00001000, /* kmer */ 0b00000000, 0b00000000,
                     0b00000000, /* edge */ 0b00000010},
                    {/* size */ 0b00001000, /* kmer */ 0b00000000, 0b00000000,
                     0b00000000, /* edge */ 0b00000011},
                    {/* size */ 0b00001000, /* kmer */ 0b01001100, 0b10100010,
                     0b00000011, /* edge */ 0b00000010},
                    {/* size */ 0b00000100, /* kmer */ 0b00000000, 0b10101000,
                     0b00000011, /* edge */ 0b00000010},
                    {/* size */ 0b00000100, /* kmer */ 0b00000000, 0b10101000,
                     0b00000011, /* edge */ 0b00000010},
                    {/* size */ 0b00000100, /* kmer */ 0b00000000, 0b11101000,
                     0b00000011, /* edge */ 0b00000010},
                    {/* size */ 0b00001000, /* kmer */ 0b00000000, 0b11101000,
                     0b00000011, /* edge */ 0b00000010},
                    {/* size */ 0b00001000, /* kmer */ 0b00000000, 0b11101000,
                     0b00000011, /* edge */ 0b00000010},
                    {/* size */ 0b00001000, /* kmer */ 0b00000000, 0b11101000,
                     0b00000011, /* edge */ 0b00000110}},
                   TerminalBuffer::autofit_tag) {}

  uint8_t ksmall, kbig;
  std::size_t num_genomes;
  TerminalBuffer small_buffer, big_buffer;
};

TEST_F(TerminalsLessThanTests, TerminalsLessThan_Small_lt_size) {
  TerminalsLessThan Cmp{small_buffer};
  ASSERT_TRUE(Cmp(0, 1));
}

TEST_F(TerminalsLessThanTests, TerminalsLessThan_Small_lt_kmer) {
  TerminalsLessThan Cmp{small_buffer};
  ASSERT_TRUE(Cmp(1, 2));
}

TEST_F(TerminalsLessThanTests, TerminalsLessThan_Small_lt_size_and_kmer) {
  TerminalsLessThan Cmp{small_buffer};
  ASSERT_TRUE(Cmp(2, 3));
}

TEST_F(TerminalsLessThanTests, TerminalsLessThan_Small_eq) {
  TerminalsLessThan Cmp{small_buffer};
  ASSERT_FALSE(Cmp(3, 4));
}

TEST_F(TerminalsLessThanTests, TerminalsLessThan_Small_gt) {
  TerminalsLessThan Cmp{small_buffer};
  ASSERT_FALSE(Cmp(4, 5));
}

TEST_F(TerminalsLessThanTests, TerminalsLessThan_Large_lt_size) {
  TerminalsLessThan Cmp{big_buffer};
  ASSERT_TRUE(Cmp(0, 1));
}

TEST_F(TerminalsLessThanTests, TerminalsLessThan_Large_lt_edge) {
  TerminalsLessThan Cmp{big_buffer};
  ASSERT_TRUE(Cmp(1, 2));
}

TEST_F(TerminalsLessThanTests, TerminalsLessThan_Large_lt_kmer) {
  TerminalsLessThan Cmp{big_buffer};
  ASSERT_TRUE(Cmp(3, 4));
}

TEST_F(TerminalsLessThanTests, TerminalsLessThan_Large_eq) {
  TerminalsLessThan Cmp{big_buffer};
  ASSERT_FALSE(Cmp(4, 5));
}

TEST_F(TerminalsLessThanTests, TerminalsLessThan_Large_gt_kmer) {
  TerminalsLessThan Cmp{big_buffer};
  ASSERT_FALSE(Cmp(6, 5));
}

TEST_F(TerminalsLessThanTests, TerminalsLessThan_Large_gt_size) {
  TerminalsLessThan Cmp{big_buffer};
  ASSERT_FALSE(Cmp(7, 6));
}

class TerminalSuffixEqTests : public testing::Test {
protected:
  TerminalSuffixEqTests() {}
};

TEST_F(TerminalSuffixEqTests, SmallKSmallS) {
  TerminalIsSuffix Eq(3, 1);

  uint64_t sfx = 0b11;

  uint8_t seq1 = 0b00110000, seq2 = 0b00111010, seq3 = 0b00100000;
  ASSERT_TRUE(Eq(&seq1, sfx));
  ASSERT_FALSE(Eq(&seq2, sfx));
  ASSERT_FALSE(Eq(&seq3, sfx));
}

TEST_F(TerminalSuffixEqTests, SmallKBigS) {
  TerminalIsSuffix Eq(3, 2);

  uint64_t sfx = 0b1110;

  uint8_t seq1 = 0b00111000, seq2 = 0b00111010, seq3 = 0b00100000;
  ASSERT_TRUE(Eq(&seq1, sfx));
  ASSERT_FALSE(Eq(&seq2, sfx));
  ASSERT_FALSE(Eq(&seq3, sfx));
}

TEST_F(TerminalSuffixEqTests, BigKSmallS) {
  TerminalIsSuffix Eq(9, 1);

  uint64_t sfx = 0b11;
  uint8_t seq1[3] = {0b00000000, 0b00000000, 0b00000011},
          seq2[3] = {0b10000000, 0b00000000, 0b00000011};

  ASSERT_TRUE(Eq(seq1, sfx));
  ASSERT_FALSE(Eq(seq2, sfx));
}

TEST_F(TerminalSuffixEqTests, BigKBigS) {
  TerminalIsSuffix Eq(9, 3);

  uint64_t sfx = 0b110111;
  uint8_t seq1[3] = {0b00000000, 0b01110000, 0b00000011},
          seq2[3] = {0b00000000, 0b00111000, 0b00000011};

  ASSERT_TRUE(Eq(seq1, sfx));
  ASSERT_FALSE(Eq(seq2, sfx));
}

class UnalignedBytesCmpTests : public testing::Test {
protected:
  UnalignedBytesCmpTests()
      : ksmall(3), kbig(9), ssmall(1), sbig(2), kwsmall(key_size(ksmall)),
        kwbig(key_size(kbig)) {}

  uint8_t ksmall, kbig, ssmall, sbig, kwsmall, kwbig;
};

TEST_F(UnalignedBytesCmpTests, MaskLSB_kSmallsSmall_lt) {
  UnalignedBytesLessThan<true> Cmp{ksmall, ssmall};
  std::vector<uint8_t> v1 = {0b00010000};
  ASSERT_TRUE(Cmp(ndim::Span<const uint8_t *>(v1.data(), kwsmall), 0b0011));
  ASSERT_FALSE(Cmp(0b0011, ndim::Span<const uint8_t *>(v1.data(), kwsmall)));
}

TEST_F(UnalignedBytesCmpTests, MaskLSB_kSmallsSmall_eq) {
  UnalignedBytesLessThan<true> Cmp{ksmall, ssmall};
  std::vector<uint8_t> v1 = {0b00110000};
  ASSERT_FALSE(Cmp(ndim::Span<const uint8_t *>(v1.data(), kwsmall), 0b0011));
  ASSERT_FALSE(Cmp(0b0011, ndim::Span<const uint8_t *>(v1.data(), kwsmall)));
}

TEST_F(UnalignedBytesCmpTests, MaskLSB_kSmallsSmall_gt) {
  UnalignedBytesLessThan<true> Cmp{ksmall, ssmall};
  std::vector<uint8_t> v1 = {0b00110000};
  ASSERT_FALSE(Cmp(ndim::Span<const uint8_t *>(v1.data(), kwsmall), 0b0001));
  ASSERT_TRUE(Cmp(0b0001, ndim::Span<const uint8_t *>(v1.data(), kwsmall)));
}

TEST_F(UnalignedBytesCmpTests, MaskLSB_kBigsSmall_lt) {
  UnalignedBytesLessThan<true> Cmp{kbig, ssmall};
  std::vector<uint8_t> v1 = {0b00110000, 0b00111111, 0b00000001};
  ASSERT_TRUE(Cmp(ndim::Span<const uint8_t *>(v1.data(), kbig), 0b0011));
  ASSERT_FALSE(Cmp(0b0011, ndim::Span<const uint8_t *>(v1.data(), kbig)));
}

TEST_F(UnalignedBytesCmpTests, MaskLSB_kBigsSmall_eq) {
  UnalignedBytesLessThan<true> Cmp{kbig, ssmall};
  std::vector<uint8_t> v1 = {0b00110000, 0b00111111, 0b00000001};
  ASSERT_FALSE(Cmp(ndim::Span<const uint8_t *>(v1.data(), kbig), 0b0001));
  ASSERT_FALSE(Cmp(0b0001, ndim::Span<const uint8_t *>(v1.data(), kbig)));
}

TEST_F(UnalignedBytesCmpTests, MaskLSB_kBigsSmall_gt) {
  UnalignedBytesLessThan<true> Cmp{kbig, ssmall};
  std::vector<uint8_t> v1 = {0b00110000, 0b00111111, 0b00000001};
  ASSERT_FALSE(Cmp(ndim::Span<const uint8_t *>(v1.data(), kbig), 0b0000));
  ASSERT_TRUE(Cmp(0b0000, ndim::Span<const uint8_t *>(v1.data(), kbig)));
}

TEST_F(UnalignedBytesCmpTests, MaskLSB_kBigsBig_lt) {
  UnalignedBytesLessThan<true> Cmp{kbig, sbig};
  std::vector<uint8_t> v1 = {0b00110000, 0b00111111, 0b00000001};
  ASSERT_TRUE(Cmp(ndim::Span<const uint8_t *>(v1.data(), kbig), 0b0101));
  ASSERT_FALSE(Cmp(0b0101, ndim::Span<const uint8_t *>(v1.data(), kbig)));
}

TEST_F(UnalignedBytesCmpTests, MaskLSB_kBigsBig_eq) {
  UnalignedBytesLessThan<true> Cmp{kbig, sbig};
  std::vector<uint8_t> v1 = {0b00110000, 0b00111111, 0b00000001};
  ASSERT_FALSE(Cmp(ndim::Span<const uint8_t *>(v1.data(), kbig), 0b0100));
  ASSERT_FALSE(Cmp(0b0100, ndim::Span<const uint8_t *>(v1.data(), kbig)));
}

TEST_F(UnalignedBytesCmpTests, MaskLSB_kBigsBig_gt) {
  UnalignedBytesLessThan<true> Cmp{kbig, sbig};
  std::vector<uint8_t> v1 = {0b00110000, 0b00111111, 0b00000001};
  ASSERT_FALSE(Cmp(ndim::Span<const uint8_t *>(v1.data(), kbig), 0b0001));
  ASSERT_TRUE(Cmp(0b0001, ndim::Span<const uint8_t *>(v1.data(), kbig)));
}

TEST_F(UnalignedBytesCmpTests, MaskedAlignedBytesLessThan_lt) {
  MaskedAlignedBytesLessThan<LHS> Cmp{kbig};
  std::vector<uint8_t> term = {0b00110000, 0b00111111, 0b11111101, 0b00001111},
                       kmer = {0b00110001, 0b00111111, 0b00000001};
  ASSERT_TRUE(Cmp(term.data(), kmer.data()));
}

TEST_F(UnalignedBytesCmpTests, MaskedAlignedBytesLessThan_eq) {
  MaskedAlignedBytesLessThan<LHS> Cmp{kbig};
  std::vector<uint8_t> term = {0b00110001, 0b00111111, 0b11111101, 0b00001111},
                       kmer = {0b00110001, 0b00111111, 0b00000001};
  ASSERT_FALSE(Cmp(term.data(), kmer.data()));
}

TEST_F(UnalignedBytesCmpTests, MaskedAlignedBytesLessThan_gt) {
  MaskedAlignedBytesLessThan<LHS> Cmp{kbig};
  std::vector<uint8_t> term = {0b00110010, 0b01111111, 0b11111101, 0b00001111},
                       kmer = {0b00111001, 0b00111111, 0b00000001};
  ASSERT_FALSE(Cmp(term.data(), kmer.data()));
}

class TerminalDiffTests : public testing::Test {
protected:
  TerminalDiffTests() {}
};

TEST_F(TerminalDiffTests, SmallTermEqual) {
  std::vector<uint8_t> t1 = {3, 0b00101010}, t2 = {3, 0b00101010};
  auto res = TerminalDiff{1, 3}(t1.data(), t2.data());
  ASSERT_EQ(res, IS_0);
}

TEST_F(TerminalDiffTests, SmallTermSeq) {
  std::vector<uint8_t> t1 = {3, 0b00101010}, t2 = {3, 0b00100000};
  auto res = TerminalDiff{1, 3}(t1.data(), t2.data());
  ASSERT_EQ(res, BW_0_K);
}

TEST_F(TerminalDiffTests, SmallTermSize) {
  std::vector<uint8_t> t1 = {3, 0b00100000}, t2 = {1, 0b00100000};
  auto res = TerminalDiff{1, 3}(t1.data(), t2.data());
  ASSERT_EQ(res, BW_0_K);
}

TEST_F(TerminalDiffTests, SmallTermSizeClose) {
  std::vector<uint8_t> t1 = {3, 0b00101000}, t2 = {2, 0b00101000};
  auto res = TerminalDiff{1, 3}(t1.data(), t2.data());
  ASSERT_EQ(res, IS_K);
}

TEST_F(TerminalDiffTests, BigTermEqual) {
  std::vector<uint8_t> t1 = {9, 0b00001100, 0b00111100, 0b00001100},
                       t2 = {9, 0b00001100, 0b00111100, 0b00001100};
  auto res = TerminalDiff{1, 10}(t1.data(), t2.data());
  ASSERT_EQ(res, IS_0);
}

TEST_F(TerminalDiffTests, BigTermSeq) {
  std::vector<uint8_t> t1 = {10, 0b00000100, 0b00000000, 0b00001000},
                       t2 = {9, 0b00000000, 0b00000000, 0b00001000};
  auto res = TerminalDiff{1, 10}(t1.data(), t2.data());
  ASSERT_EQ(res, BW_0_K);
}

TEST_F(TerminalDiffTests, BigTermSize) {
  std::vector<uint8_t> t1 = {10, 0b00000000, 0b00000000, 0b00001100},
                       t2 = {1, 0b00000000, 0b00000000, 0b00001100};
  auto res = TerminalDiff{1, 10}(t1.data(), t2.data());
  ASSERT_EQ(res, BW_0_K);
}

TEST_F(TerminalDiffTests, BigTermSizeClose) {
  std::vector<uint8_t> t1 = {10, 0b11111100, 0b00001010, 0b00001010},
                       t2 = {9, 0b11111100, 0b00001010, 0b00001010};
  auto res = TerminalDiff{1, 10}(t1.data(), t2.data());
  ASSERT_EQ(res, IS_K);
}

TEST_F(TerminalDiffTests, BigTermSeqClose) {
  std::vector<uint8_t> t1 = {10, 0b11111100, 0b00001010, 0b00001010},
                       t2 = {10, 0b11111101, 0b00001010, 0b00001010};
  auto res = TerminalDiff{1, 10}(t1.data(), t2.data());
  ASSERT_EQ(res, IS_K);
}

// counting

/*

~seqA

>>>
ACGAGNCGGTGGCAANGAAGTTN 9
<<<
ACGAGNCGGTGGCAANGAAGTTN 9
>>>
GTGTCGGAGGNCTCCATCGAC 6
<<<
GTGTCGGAGGNCTCCATCGAC 6

~ seqB

>>>
TACACTTACTCG 3
<<<
TACACTTACTCG 3
>>>
TACTCGGACTCA 3
<<<
TACTCGGACTCA 3
>>>
GACTCAGACTCA 3
<<<
GACTCAGACTCA 3
>>>
ACGNACNACGT 6
<<<
ACGNACNACGT 6

~ seqC

>>>
AGCTTTTCATTCTGACTG 3
<<<
AGCTTTTCATTCTGACTG 3

*/

TEST(GenomeCountTests, CountGenomesGFF) {
  std::size_t k = 3;
  Colours c;
  auto genomes = parse({STRING(SEQA_GFF), STRING(SEQB_GFF), STRING(SEQC_GFF)},
                       parseGFF, c, (std::size_t)k);
  auto view = toView(genomes);
  auto blocks = planTerminalRanges(view, k);
  ASSERT_THAT(blocks, ElementsAreArray({30, 60, 66}));
}

/*

>genome1 test sequence
GTGTCGGAGGCTCCATCGACATGGAACGAGCGGTGGCAAGAAGTTACTAATGAGCTGCTG

>genome2 other test sequence
GTGTCGGAGGCTCCATCGACATGGAACGAGCGGTGGCAAGAAGTTACTAATGAGCTGCTG

>genome1 test sequence
GTGTCGGAGGCTCCATCGACATGGAACGAGCGGTGGCAAGAAGTTACTAATGAGCTGCTG

>genome2 other test sequence
GTGTCGGAGGCTCCATCGACATGGAACGAGCGGTGGCAAGAAGTTACTAATGAGCTGCTG


*/

class CountChunkedGenomeTests : public testing::Test {
protected:
  CountChunkedGenomeTests() : k(3), fna(MakeTempPath(".fna")), c() {}
  ~CountChunkedGenomeTests() { std::filesystem::remove(fna); }

  std::size_t k;
  std::string fna;
  Colours c;

  ChunkedDna4Genome readGenomeData(std::size_t granularity) {
    std::string data = ">genome1 test "
                       "sequence\nGTGTCGGAGGCTCCATCGACATGGAACGAGCGGTGGCAAGAAGTT"
                       "ACTAATGAGCTGCTG\n>genome2 other test "
                       "sequence\nGTGTCGGAGGCTCCATCGACATGGAACGAGCGGTGGCAAGAAGTT"
                       "ACTAATGAGCTGCTG\n>genome1 test "
                       "sequence\nGTGTCGGAGGCTCCATCGACATGGAACGAGCGGTGGCAAGAAGTT"
                       "ACTAATGAGCTGCTG\n>genome2 other test "
                       "sequence\nGTGTCGGAGGCTCCATCGACATGGAACGAGCGGTGGCAAGAAGTT"
                       "ACTAATGAGCTGCTG\n";
    writeToFna(data, fna);
    return parseFilterFNA(fna, c, granularity, k);
  }

  std::vector<const SequenceContainer *>
  localView(const ChunkedDna4Genome &xgenome) {
    std::vector<const SequenceContainer *> res;
    for (const auto &chunk : xgenome.chunks) {
      res.push_back(&chunk);
    }
    return res;
  }

  void checkTerminals(const TerminalBuffer &buffer) {
    auto expected =
        ndim::Matrix({{/* size */ 0, /* kmer */ 0b00000000, /* value */ 0b10},
                      {/* size */ 1, /* kmer */ 0b00100000, /* value */ 0b11},
                      {/* size */ 2, /* kmer */ 0b00111000, /* value */ 0b10},
                      {/* size */ 0, /* kmer */ 0b00000000, /* value */ 0b01},
                      {/* size */ 1, /* kmer */ 0b00010000, /* value */ 0b00},
                      {/* size */ 2, /* kmer */ 0b00000100, /* value */ 0b10},
                      {/* size */ 0, /* kmer */ 0b00000000, /* value */ 0b10},
                      {/* size */ 1, /* kmer */ 0b00100000, /* value */ 0b11},
                      {/* size */ 2, /* kmer */ 0b00111000, /* value */ 0b10},
                      {/* size */ 0, /* kmer */ 0b00000000, /* value */ 0b01},
                      {/* size */ 1, /* kmer */ 0b00010000, /* value */ 0b00},
                      {/* size */ 2, /* kmer */ 0b00000100, /* value */ 0b10},
                      {/* size */ 0, /* kmer */ 0b00000000, /* value */ 0b10},
                      {/* size */ 1, /* kmer */ 0b00100000, /* value */ 0b11},
                      {/* size */ 2, /* kmer */ 0b00111000, /* value */ 0b10},
                      {/* size */ 0, /* kmer */ 0b00000000, /* value */ 0b01},
                      {/* size */ 1, /* kmer */ 0b00010000, /* value */ 0b00},
                      {/* size */ 2, /* kmer */ 0b00000100, /* value */ 0b10},
                      {/* size */ 0, /* kmer */ 0b00000000, /* value */ 0b10},
                      {/* size */ 1, /* kmer */ 0b00100000, /* value */ 0b11},
                      {/* size */ 2, /* kmer */ 0b00111000, /* value */ 0b10},
                      {/* size */ 0, /* kmer */ 0b00000000, /* value */ 0b01},
                      {/* size */ 1, /* kmer */ 0b00010000, /* value */ 0b00},
                      {/* size */ 2, /* kmer */ 0b00000100, /* value */ 0b10}});
    ASSERT_THAT(buffer.memoryview(), MatrixEq(std::cref(expected)));
  }
};

TEST_F(CountChunkedGenomeTests, Count2) {
  auto data = readGenomeData(2);
  auto view = localView(data);
  auto blocks = planTerminalRanges(view, k);
  ASSERT_EQ(blocks.back(), 24);
}

TEST_F(CountChunkedGenomeTests, Count3) {
  auto data = readGenomeData(3);
  auto view = localView(data);
  auto blocks = planTerminalRanges(view, k);
  ASSERT_EQ(blocks.back(), 24);
}

TEST_F(CountChunkedGenomeTests, Count100) {
  auto data = readGenomeData(100);
  auto view = localView(data);
  auto blocks = planTerminalRanges(view, k);
  ASSERT_EQ(blocks.back(), 24);
}
// insert

class TerminalBufferCreationTests : public testing::Test {
public:
  TerminalBufferCreationTests()
      : k(9), c(), genomes(2), num_terminals(/* num. fragments */8 * k),
        buffer(num_terminals, k, TerminalBuffer::autofit_tag) {
    ParseFastaToGenome(genomes[0], STRING(SMALL_SEQ), c, k);
    ParseFastaToGenome(genomes[1], STRING(SMALL_SEQ), c, k);
  }

protected:
  uint8_t k;
  Colours c;
  std::vector<Dna4Genome> genomes;
  std::size_t num_terminals;
  TerminalBuffer buffer;
};

/**

Contig: GTGTCGGAGGCTCCATCGACATGGAACGAGCGGTGGCAAGAAGTTACTAATGAGCTGCTG
Encoding: (kmer)    (edge)  (index)
    $$$$ $$$$ ---$ (0b00000000, 0b00000000, 0b00000000)   G 0 (0b00000010)    0
(0b00000000)
    $$$$ $$$$ ---G (0b00000000, 0b00000000, 0b00000010)   T 0 (0b00000011)    1
(0b00000001)
    $$$$ G$$$ ---T (0b00000000, 0b10000000, 0b00000011)   G 0 (0b00000010)    2
(0b00000010)
    $$$$ TG$$ ---G (0b00000000, 0b11100000, 0b00000010)   T 0 (0b00000011)    3
(0b00000011)
    $$$$ GTG$ ---T (0b00000000, 0b10111000, 0b00000011)   C 0 (0b00000001)    4
(0b00000100)
    $$$$ TGTG ---C (0b00000000, 0b11101110, 0b00000001)   G 0 (0b00000010)    5
(0b00000101) G$$$ CTGT ---G (0b10000000, 0b01111011, 0b00000010)   G 0
(0b00000010)    6 (0b00000110) TG$$ GCTG ---G (0b11100000, 0b10011110,
0b00000010)   A 0 (0b00000000)    7 (0b00000111) GTG$ GGCT ---A (0b10111000,
0b10100111, 0b00000000)   G 0 (0b00000010)    8 (0b00001000)

 */
TEST_F(TerminalBufferCreationTests, Insert) {
  auto view = toView(genomes);

  for (auto seq : view) {
    for (auto fmt : seq->terminals()) {
      std::cout << fmt.size() << ' ';
    }
  }
  std::cout << std::endl;

  auto blocks = planTerminalRanges(view, k);
  buffer.fill(view, blocks);

  auto expected =
      ndim::Matrix({{/* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000001, /* kmer */ 0b00000000, 0b00000000,
                     0b00000010, /* value */ 0b00000011},
                    {/* size */ 0b00000010, /* kmer */ 0b00000000, 0b10000000,
                     0b00000011, /* value */ 0b00000010},
                    {/* size */ 0b00000011, /* kmer */ 0b00000000, 0b11100000,
                     0b00000010, /* value */ 0b00000011},
                    {/* size */ 0b00000100, /* kmer */ 0b00000000, 0b10111000,
                     0b00000011, /* value */ 0b00000001},
                    {/* size */ 0b00000101, /* kmer */ 0b00000000, 0b11101110,
                     0b00000001, /* value */ 0b00000010},
                    {/* size */ 0b00000110, /* kmer */ 0b10000000, 0b01111011,
                     0b00000010, /* value */ 0b00000010},
                    {/* size */ 0b00000111, /* kmer */ 0b11100000, 0b10011110,
                     0b00000010, /* value */ 0b00000000},
                    {/* size */ 0b00001000, /* kmer */ 0b10111000, 0b10100111,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000,
                     0b00000000, /* value */ 0b00000001},
                    {/* size */ 0b00000001, /* kmer */ 0b00000000, 0b00000000,
                     0b00000001, /* value */ 0b00000000},
                    {/* size */ 0b00000010, /* kmer */ 0b00000000, 0b01000000,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000011, /* kmer */ 0b00000000, 0b00010000,
                     0b00000010, /* value */ 0b00000001},
                    {/* size */ 0b00000100, /* kmer */ 0b00000000, 0b10000100,
                     0b00000001, /* value */ 0b00000000},
                    {/* size */ 0b00000101, /* kmer */ 0b00000000, 0b01100001,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000110, /* kmer */ 0b01000000, 0b00011000,
                     0b00000010, /* value */ 0b00000001},
                    {/* size */ 0b00000111, /* kmer */ 0b00010000, 0b10000110,
                     0b00000001, /* value */ 0b00000011},
                    {/* size */ 0b00001000, /* kmer */ 0b10000100, 0b01100001,
                     0b00000011, /* value */ 0b00000001},
                    {/* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000001, /* kmer */ 0b00000000, 0b00000000,
                     0b00000010, /* value */ 0b00000011},
                    {/* size */ 0b00000010, /* kmer */ 0b00000000, 0b10000000,
                     0b00000011, /* value */ 0b00000010},
                    {/* size */ 0b00000011, /* kmer */ 0b00000000, 0b11100000,
                     0b00000010, /* value */ 0b00000011},
                    {/* size */ 0b00000100, /* kmer */ 0b00000000, 0b10111000,
                     0b00000011, /* value */ 0b00000001},
                    {/* size */ 0b00000101, /* kmer */ 0b00000000, 0b11101110,
                     0b00000001, /* value */ 0b00000010},
                    {/* size */ 0b00000110, /* kmer */ 0b10000000, 0b01111011,
                     0b00000010, /* value */ 0b00000010},
                    {/* size */ 0b00000111, /* kmer */ 0b11100000, 0b10011110,
                     0b00000010, /* value */ 0b00000000},
                    {/* size */ 0b00001000, /* kmer */ 0b10111000, 0b10100111,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000,
                     0b00000000, /* value */ 0b00000001},
                    {/* size */ 0b00000001, /* kmer */ 0b00000000, 0b00000000,
                     0b00000001, /* value */ 0b00000000},
                    {/* size */ 0b00000010, /* kmer */ 0b00000000, 0b01000000,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000011, /* kmer */ 0b00000000, 0b00010000,
                     0b00000010, /* value */ 0b00000001},
                    {/* size */ 0b00000100, /* kmer */ 0b00000000, 0b10000100,
                     0b00000001, /* value */ 0b00000000},
                    {/* size */ 0b00000101, /* kmer */ 0b00000000, 0b01100001,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000110, /* kmer */ 0b01000000, 0b00011000,
                     0b00000010, /* value */ 0b00000001},
                    {/* size */ 0b00000111, /* kmer */ 0b00010000, 0b10000110,
                     0b00000001, /* value */ 0b00000011},
                    {/* size */ 0b00001000, /* kmer */ 0b10000100, 0b01100001,
                     0b00000011, /* value */ 0b00000001},
                    {/* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000001, /* kmer */ 0b00000000, 0b00000000,
                     0b00000010, /* value */ 0b00000011},
                    {/* size */ 0b00000010, /* kmer */ 0b00000000, 0b10000000,
                     0b00000011, /* value */ 0b00000010},
                    {/* size */ 0b00000011, /* kmer */ 0b00000000, 0b11100000,
                     0b00000010, /* value */ 0b00000011},
                    {/* size */ 0b00000100, /* kmer */ 0b00000000, 0b10111000,
                     0b00000011, /* value */ 0b00000001},
                    {/* size */ 0b00000101, /* kmer */ 0b00000000, 0b11101110,
                     0b00000001, /* value */ 0b00000010},
                    {/* size */ 0b00000110, /* kmer */ 0b10000000, 0b01111011,
                     0b00000010, /* value */ 0b00000010},
                    {/* size */ 0b00000111, /* kmer */ 0b11100000, 0b10011110,
                     0b00000010, /* value */ 0b00000000},
                    {/* size */ 0b00001000, /* kmer */ 0b10111000, 0b10100111,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000,
                     0b00000000, /* value */ 0b00000001},
                    {/* size */ 0b00000001, /* kmer */ 0b00000000, 0b00000000,
                     0b00000001, /* value */ 0b00000000},
                    {/* size */ 0b00000010, /* kmer */ 0b00000000, 0b01000000,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000011, /* kmer */ 0b00000000, 0b00010000,
                     0b00000010, /* value */ 0b00000001},
                    {/* size */ 0b00000100, /* kmer */ 0b00000000, 0b10000100,
                     0b00000001, /* value */ 0b00000000},
                    {/* size */ 0b00000101, /* kmer */ 0b00000000, 0b01100001,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000110, /* kmer */ 0b01000000, 0b00011000,
                     0b00000010, /* value */ 0b00000001},
                    {/* size */ 0b00000111, /* kmer */ 0b00010000, 0b10000110,
                     0b00000001, /* value */ 0b00000011},
                    {/* size */ 0b00001000, /* kmer */ 0b10000100, 0b01100001,
                     0b00000011, /* value */ 0b00000001},
                    {/* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000001, /* kmer */ 0b00000000, 0b00000000,
                     0b00000010, /* value */ 0b00000011},
                    {/* size */ 0b00000010, /* kmer */ 0b00000000, 0b10000000,
                     0b00000011, /* value */ 0b00000010},
                    {/* size */ 0b00000011, /* kmer */ 0b00000000, 0b11100000,
                     0b00000010, /* value */ 0b00000011},
                    {/* size */ 0b00000100, /* kmer */ 0b00000000, 0b10111000,
                     0b00000011, /* value */ 0b00000001},
                    {/* size */ 0b00000101, /* kmer */ 0b00000000, 0b11101110,
                     0b00000001, /* value */ 0b00000010},
                    {/* size */ 0b00000110, /* kmer */ 0b10000000, 0b01111011,
                     0b00000010, /* value */ 0b00000010},
                    {/* size */ 0b00000111, /* kmer */ 0b11100000, 0b10011110,
                     0b00000010, /* value */ 0b00000000},
                    {/* size */ 0b00001000, /* kmer */ 0b10111000, 0b10100111,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000,
                     0b00000000, /* value */ 0b00000001},
                    {/* size */ 0b00000001, /* kmer */ 0b00000000, 0b00000000,
                     0b00000001, /* value */ 0b00000000},
                    {/* size */ 0b00000010, /* kmer */ 0b00000000, 0b01000000,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000011, /* kmer */ 0b00000000, 0b00010000,
                     0b00000010, /* value */ 0b00000001},
                    {/* size */ 0b00000100, /* kmer */ 0b00000000, 0b10000100,
                     0b00000001, /* value */ 0b00000000},
                    {/* size */ 0b00000101, /* kmer */ 0b00000000, 0b01100001,
                     0b00000000, /* value */ 0b00000010},
                    {/* size */ 0b00000110, /* kmer */ 0b01000000, 0b00011000,
                     0b00000010, /* value */ 0b00000001},
                    {/* size */ 0b00000111, /* kmer */ 0b00010000, 0b10000110,
                     0b00000001, /* value */ 0b00000011},
                    {/* size */ 0b00001000, /* kmer */ 0b10000100, 0b01100001,
                     0b00000011, /* value */ 0b00000001}});

  ASSERT_THAT(buffer.memoryview(), MatrixEq(std::cref(expected)));
}

TEST_F(CountChunkedGenomeTests, ExtractTerminals2) {
  auto data = readGenomeData(2);
  auto view = localView(data);
  auto plan = planTerminalRanges(view, k);
  TerminalBuffer buffer(plan.back(), k, TerminalBuffer::autofit_tag);
  buffer.fill(view, plan);
  checkTerminals(buffer);
}

TEST_F(CountChunkedGenomeTests, ExtractTerminals3) {
  auto data = readGenomeData(3);
  auto view = localView(data);
  auto plan = planTerminalRanges(view, k);
  TerminalBuffer buffer(plan.back(), k, TerminalBuffer::autofit_tag);
  buffer.fill(view, plan);
  checkTerminals(buffer);
}

TEST_F(CountChunkedGenomeTests, ExtractTerminals100) {
  auto data = readGenomeData(100);
  auto view = localView(data);
  auto plan = planTerminalRanges(view, k);
  TerminalBuffer buffer(plan.back(), k, TerminalBuffer::autofit_tag);
  buffer.fill(view, plan);
  checkTerminals(buffer);
}

class TerminalBufferSortSubroutinesTests : public testing::Test {};

TEST_F(TerminalBufferSortSubroutinesTests, SortShort) {
  auto buffer =
      TerminalBuffer(3,
                     {
                         {/* size */ 0b00000010, /* kmer */ 0b00000000,
                          /* edge */ 0b00000010},
                         {/* size */ 0b00000010, /* kmer */ 0b00010000,
                          /* edge */ 0b00000010},
                         {/* size */ 0b00000010, /* kmer */ 0b00110000,
                          /* edge */ 0b00000010},
                         {/* size */ 0b00000010, /* kmer */ 0b00111000,
                          /* edge */ 0b00000011},
                         {/* size */ 0b00000010, /* kmer */ 0b00000000,
                          /* edge */ 0b00000000},
                         {/* size */ 0b00000001, /* kmer */ 0b00110000,
                          /* edge */ 0b00000010},
                         {/* size */ 0b00000000, /* kmer */ 0b00000000,
                          /* edge */ 0b00000000},
                         {/* size */ 0b00000010, /* kmer */ 0b00000000,
                          /* edge */ 0b00000011},
                         {/* size */ 0b00000010, /* kmer */ 0b00110100,
                          /* edge */ 0b00000010},
                         {/* size */ 0b00000010, /* kmer */ 0b00111000,
                          /* edge */ 0b00000010},
                     },
                     TerminalBuffer::autofit_tag);

  std::vector<size_t> indices(buffer.size());
  std::iota(indices.begin(), indices.end(), 0ul);
  std::sort(indices.begin(), indices.end(), TerminalsLessThan{buffer});
  ASSERT_THAT(indices, ElementsAreArray({6, 4, 0, 7, 1, 5, 2, 8, 9, 3}));
}

TEST_F(TerminalBufferSortSubroutinesTests, SortLong) {
  auto buffer =
      TerminalBuffer(9,
                     {
                         {/* size */ 0b00001000, /* kmer */ 0b00111000,
                          0b11100000, 0b00000000, /* edge */ 0b00000010},
                         {/* size */ 0b00000100, /* kmer */ 0b00000000,
                          0b11110000, 0b00000000, /* edge */ 0b00000010},
                         {/* size */ 0b00000010, /* kmer */ 0b00000000,
                          0b00000000, 0b00000010, /* edge */ 0b00000000},
                         {/* size */ 0b00000000, /* kmer */ 0b00000000,
                          0b00000000, 0b00000000, /* edge */ 0b00000010},
                         {/* size */ 0b00000001, /* kmer */ 0b00000000,
                          0b00000000, 0b00000010, /* edge */ 0b00000000},
                         {/* size */ 0b00000100, /* kmer */ 0b00000000,
                          0b11110000, 0b00000000, /* edge */ 0b00000000},
                         {/* size */ 0b00001000, /* kmer */ 0b00000000,
                          0b00000000, 0b00000010, /* edge */ 0b00000010},
                         {/* size */ 0b00000010, /* kmer */ 0b00000000,
                          0b00000000, 0b00000010, /* edge */ 0b00000011},
                         {/* size */ 0b00000011, /* kmer */ 0b00000000,
                          0b11110000, 0b00000000, /* edge */ 0b00000010},
                         {/* size */ 0b00000100, /* kmer */ 0b00000000,
                          0b11110000, 0b00000000, /* edge */ 0b00000001},
                     },
                     TerminalBuffer::autofit_tag);

  std::vector<size_t> indices(buffer.size());
  std::iota(indices.begin(), indices.end(), 0ul);

  std::sort(indices.begin(), indices.end(), TerminalsLessThan{buffer});

  ASSERT_THAT(indices, ElementsAreArray({3, 0, 8, 5, 9, 1, 4, 2, 7, 6}));
}

TEST_F(TerminalBufferSortSubroutinesTests, AdjacentDifference) {
  auto buffer =
      TerminalBuffer(9,
                     {
                         {/* size */ 0b00000000, /* kmer */ 0b00000000,
                          0b00000000, 0b00000000, /* edge */ 0b00000000},
                         {/* size */ 0b00000001, /* kmer */ 0b00000000,
                          0b00000000, 0b00000000, /* edge */ 0b00000000},
                         {/* size */ 0b00000001, /* kmer */ 0b00000000,
                          0b00000000, 0b00000000, /* edge */ 0b00000000},
                         {/* size */ 0b00000001, /* kmer */ 0b00000000,
                          0b00000000, 0b00000000, /* edge */ 0b00000000},
                         {/* size */ 0b00000010, /* kmer */ 0b00000000,
                          0b00000000, 0b00000000, /* edge */ 0b00000000},
                         {/* size */ 0b00000010, /* kmer */ 0b00000000,
                          0b00000000, 0b00000000, /* edge */ 0b00000000},
                         {/* size */ 0b00000010, /* kmer */ 0b00000000,
                          0b00000000, 0b00000100, /* edge */ 0b00000000},
                         {/* size */ 0b00000010, /* kmer */ 0b00000000,
                          0b00000000, 0b00000100, /* edge */ 0b00000010},
                         {/* size */ 0b00000010, /* kmer */ 0b00000000,
                          0b00000000, 0b00000100, /* edge */ 0b00000001},
                         {/* size */ 0b00000010, /* kmer */ 0b00000000,
                          0b00000000, 0b00000100, /* edge */ 0b00000010},
                     },
                     TerminalBuffer::autofit_tag);

  auto B = adjacentDifference(buffer);
  ASSERT_THAT(B, ElementsAreArray({BW_0_K, BW_0_K, IS_0, IS_0, BW_0_K, IS_0,
                                   BW_0_K, IS_0, IS_0, IS_0}));
}

/**
 * This test should cover cases:
 *      - duplicate terminals (including genome-edge value and size)
 *      - sorting k-mers with same encoding but less size, and then same
 * terminal but less edge
 *
 * Example 1
 * ---------
 *
 *      A$$$ TGAC ---T  ->  T
 *      AA$$ TGAC ---T  ->  A
 *      AAA$ TGAC ---T  ->  A
 *      AAA$ TGAC ---T  ->  C
 *
 *      Sequences (sequence IDs per insert)
 *      -----------------------------------
 *      TTGACATAA (0, 1)
 *      TTGACAAAA (0, 0)
 *      TTGACAAAC (0)
 *
 *      Structure
 *      ---------
 *      $$$$ $$$$ ---$ (0b00000000, 0b00000000, 0b00000000) T 0 (0b00000011) 0
 * (0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ $$$$ ---T (0b00000000, 0b00000000, 0b00000011) T 0 (0b00000011) 1
 * (0b00000001, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ T$$$ ---T (0b00000000, 0b11000000, 0b00000011) G 0 (0b00000010) 2
 * (0b00000010, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ TT$$ ---G (0b00000000, 0b11110000, 0b00000010) A 0 (0b00000000) 3
 * (0b00000011, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ GTT$ ---A (0b00000000, 0b10111100, 0b00000000) C 0 (0b00000001) 4
 * (0b00000100, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ AGTT ---C (0b00000000, 0b00101111, 0b00000001) A 0 (0b00000000) 5
 * (0b00000101, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000) T$$$ CAGT ---A (0b11000000, 0b01001011, 0b00000000) T
 * 0 (0b00000011) 6 (0b00000110, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000, 0b00000000) TT$$ ACAG ---T (0b11110000, 0b00010010,
 * 0b00000011) A 0 (0b00000000) 7 (0b00000111, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000) GTT$ TACA ---A
 * (0b10111100, 0b11000100, 0b00000000) A 0 (0b00000000) 8 (0b00001000,
 * 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000)
 *
 *      $$$$ $$$$ ---$ (0b00000000, 0b00000000, 0b00000000) T 1 (0b00001011) 0
 * (0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ $$$$ ---T (0b00000000, 0b00000000, 0b00000011) T 1 (0b00001011) 1
 * (0b00000001, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ T$$$ ---T (0b00000000, 0b11000000, 0b00000011) G 1 (0b00001010) 2
 * (0b00000010, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ TT$$ ---G (0b00000000, 0b11110000, 0b00000010) A 1 (0b00001000) 3
 * (0b00000011, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ GTT$ ---A (0b00000000, 0b10111100, 0b00000000) C 1 (0b00001001) 4
 * (0b00000100, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ AGTT ---C (0b00000000, 0b00101111, 0b00000001) A 1 (0b00001000) 5
 * (0b00000101, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000) T$$$ CAGT ---A (0b11000000, 0b01001011, 0b00000000) T
 * 1 (0b00001011) 6 (0b00000110, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000, 0b00000000) TT$$ ACAG ---T (0b11110000, 0b00010010,
 * 0b00000011) A 1 (0b00001000) 7 (0b00000111, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000) GTT$ TACA ---A
 * (0b10111100, 0b11000100, 0b00000000) A 1 (0b00001000) 8 (0b00001000,
 * 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000)
 *
 *      $$$$ $$$$ ---$ (0b00000000, 0b00000000, 0b00000000) T 0 (0b00000011) 0
 * (0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ $$$$ ---T (0b00000000, 0b00000000, 0b00000011) T 0 (0b00000011) 1
 * (0b00000001, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ T$$$ ---T (0b00000000, 0b11000000, 0b00000011) G 0 (0b00000010) 2
 * (0b00000010, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ TT$$ ---G (0b00000000, 0b11110000, 0b00000010) A 0 (0b00000000) 3
 * (0b00000011, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ GTT$ ---A (0b00000000, 0b10111100, 0b00000000) C 0 (0b00000001) 4
 * (0b00000100, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ AGTT ---C (0b00000000, 0b00101111, 0b00000001) A 0 (0b00000000) 5
 * (0b00000101, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000) T$$$ CAGT ---A (0b11000000, 0b01001011, 0b00000000) A
 * 0 (0b00000000) 6 (0b00000110, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000, 0b00000000) TT$$ ACAG ---A (0b11110000, 0b00010010,
 * 0b00000000) A 0 (0b00000000) 7 (0b00000111, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000) GTT$ AACA ---A
 * (0b10111100, 0b00000100, 0b00000000) A 0 (0b00000000) 8 (0b00001000,
 * 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000)
 *
 *      $$$$ $$$$ ---$ (0b00000000, 0b00000000, 0b00000000) T 0 (0b00000011) 0
 * (0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ $$$$ ---T (0b00000000, 0b00000000, 0b00000011) T 0 (0b00000011) 1
 * (0b00000001, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ T$$$ ---T (0b00000000, 0b11000000, 0b00000011) G 0 (0b00000010) 2
 * (0b00000010, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ TT$$ ---G (0b00000000, 0b11110000, 0b00000010) A 0 (0b00000000) 3
 * (0b00000011, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ GTT$ ---A (0b00000000, 0b10111100, 0b00000000) C 0 (0b00000001) 4
 * (0b00000100, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000)
 *      $$$$ AGTT ---C (0b00000000, 0b00101111, 0b00000001) A 0 (0b00000000) 5
 * (0b00000101, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000) T$$$ CAGT ---A (0b11000000, 0b01001011, 0b00000000) A
 * 0 (0b00000000) 6 (0b00000110, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000, 0b00000000) TT$$ ACAG ---A (0b11110000, 0b00010010,
 * 0b00000000) A 0 (0b00000000) 7 (0b00000111, 0b00000000, 0b00000000,
 * 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000) GTT$ AACA ---A
 * (0b10111100, 0b00000100, 0b00000000) C 0 (0b00000001) 8 (0b00001000,
 * 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000, 0b00000000,
 * 0b00000000)
 *
 */
class TerminalBufferSortTests : public testing::Test {
protected:
  TerminalBufferSortTests()
      : k(9), seq1("TTGACATAA"_dna4), seq2("TTGACAAAA"_dna4),
        seq3("TTGACAAAC"_dna4), num_terminals(45),
        buffer(num_terminals, k, TerminalBuffer::autofit_tag),
        expected({{/* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000,
                   0b00000000, /* edge */ 0b00000011},
                  {/* size */ 0b00001000, /* kmer */ 0b10111100, 0b00000100,
                   0b00000000, /* edge */ 0b00000000},
                  {/* size */ 0b00001000, /* kmer */ 0b10111100, 0b00000100,
                   0b00000000, /* edge */ 0b00000001},
                  {/* size */ 0b00000111, /* kmer */ 0b11110000, 0b00010010,
                   0b00000000, /* edge */ 0b00000000},
                  {/* size */ 0b00000110, /* kmer */ 0b11000000, 0b01001011,
                   0b00000000, /* edge */ 0b00000000},
                  {/* size */ 0b00000110, /* kmer */ 0b11000000, 0b01001011,
                   0b00000000, /* edge */ 0b00000011},
                  {/* size */ 0b00000100, /* kmer */ 0b00000000, 0b10111100,
                   0b00000000, /* edge */ 0b00000001},
                  {/* size */ 0b00001000, /* kmer */ 0b10111100, 0b11000100,
                   0b00000000, /* edge */ 0b00000000},
                  {/* size */ 0b00000101, /* kmer */ 0b00000000, 0b00101111,
                   0b00000001, /* edge */ 0b00000000},
                  {/* size */ 0b00000011, /* kmer */ 0b00000000, 0b11110000,
                   0b00000010, /* edge */ 0b00000000},
                  {/* size */ 0b00000001, /* kmer */ 0b00000000, 0b00000000,
                   0b00000011, /* edge */ 0b00000011},
                  {/* size */ 0b00000111, /* kmer */ 0b11110000, 0b00010010,
                   0b00000011, /* edge */ 0b00000000},
                  {/* size */ 0b00000010, /* kmer */ 0b00000000, 0b11000000,
                   0b00000011, /* edge */ 0b00000010}}) {}

  uint8_t k;
  Dna4Sequence seq1, seq2, seq3;
  std::size_t num_terminals;
  TerminalBuffer buffer;
  ndim::Matrix<> expected;
};

TEST_F(TerminalBufferSortTests, Sort) {
  auto p = buffer.begin();
  p = buffer.insert(p, seq1.cbegin(), seq1.cbegin() + k, false);
  p = buffer.insert(p, seq1.cbegin(), seq1.cbegin() + k, false);
  p = buffer.insert(p, seq2.cbegin(), seq2.cbegin() + k, false);
  p = buffer.insert(p, seq2.cbegin(), seq2.cbegin() + k, false);
  p = buffer.insert(p, seq3.cbegin(), seq3.cbegin() + k, false);
  ASSERT_EQ(p, buffer.end());

  auto result = buffer.OOPsort();
  CheckMatrixEq(result, expected);
}

class TerminalBufferSearchTests : public testing::Test {
protected:
  TerminalBufferSearchTests()
      : k(9), num_genomes(2),
        buffer(k,
               {{/* size */ 0b00000001, /* kmer */ 0b00000000, 0b00000000,
                 0b00000000, /* edge */ 0b00000000},
                {/* size */ 0b00000110, /* kmer */ 0b00000000, 0b10101000,
                 0b00000010, /* edge */ 0b00000000},
                {/* size */ 0b00000110, /* kmer */ 0b10000000, 0b10101010,
                 0b00000010, /* edge */ 0b00000000},
                {/* size */ 0b00000110, /* kmer */ 0b10000000, 0b10101010,
                 0b00000010, /* edge */ 0b00000001},
                {/* size */ 0b00000110, /* kmer */ 0b10000000, 0b10101010,
                 0b00000010, /* edge */ 0b00000010},
                {/* size */ 0b00000110, /* kmer */ 0b10000000, 0b10101010,
                 0b00000010, /* edge */ 0b00000011},
                {/* size */ 0b00000110, /* kmer */ 0b10000000, 0b10101010,
                 0b00000010, /* edge */ 0b00000100},
                {/* size */ 0b00000110, /* kmer */ 0b10000000, 0b10101010,
                 0b00000010, /* edge */ 0b00000101},
                {/* size */ 0b00001000, /* kmer */ 0b10001000, 0b10101010,
                 0b00000010, /* edge */ 0b00000110},
                {/* size */ 0b00001000, /* kmer */ 0b10001000, 0b10101010,
                 0b00000010, /* edge */ 0b00000111},
                {/* size */ 0b00001000, /* kmer */ 0b11111000, 0b11111111,
                 0b00000011, /* edge */ 0b00000111}},
               TerminalBuffer::autofit_tag),
        range(buffer.asRange()),
        subrange(buffer.getEffK(), buffer.lengthBytes(), buffer.constAt(2),
                 buffer.constAt(10)) {}

  uint8_t k;
  std::size_t num_genomes;
  TerminalBuffer buffer;
  TerminalRange range, subrange;
};

TEST_F(TerminalBufferSearchTests, LowerBound_CommonTerminal_Seen) {
  std::vector<uint8_t> qry = {0b10000000, 0b10101010, 0b00000010};
  auto it = range.lowerBound({qry.data(), 3}, k);
  ASSERT_EQ(it, buffer.constAt(2));
}

TEST_F(TerminalBufferSearchTests, LowerBound_CommonTerminal_NotSeen) {
  std::vector<uint8_t> qry = {0b00000000, 0b10101010, 0b00000010};
  auto it = range.lowerBound({qry.data(), 3}, k);
  ASSERT_EQ(it, buffer.constAt(2));
}

TEST_F(TerminalBufferSearchTests, LowerBound_Begin) {
  std::vector<uint8_t> qry = {0b00000000, 0b00000000, 0b00000000};
  auto it = range.lowerBound({qry.data(), 3}, k);
  ASSERT_EQ(it, buffer.constBegin());
}

TEST_F(TerminalBufferSearchTests, LowerBound_End) {
  std::vector<uint8_t> qry = {0b11111100, 0b11111111, 0b00000011};
  auto it = range.lowerBound({qry.data(), 3}, k);
  ASSERT_EQ(it, buffer.constEnd());
}

TEST_F(TerminalBufferSearchTests, LowerBound_Mask_Short) {
  std::vector<uint8_t> qry = {0b10000000, 0b00001010};
  auto it = subrange.lowerBound({qry.data(), 2}, 6);
  ASSERT_EQ(it, buffer.constAt(2));
}

TEST_F(TerminalBufferSearchTests, LowerBound_Mask_VeryShort) {
  std::vector<uint8_t> qry = {0b10000001};
  auto it = subrange.lowerBound({qry.data(), 1}, 4);
  ASSERT_EQ(it, buffer.constAt(8));
}

TEST_F(TerminalBufferSearchTests, LowerBound_Mask_VeryShort_Truncated) {
  std::vector<uint8_t> qry = {0b00000001};
  auto it = subrange.lowerBound({qry.data(), 1}, 3);
  ASSERT_EQ(it, buffer.constAt(8));
}

TEST_F(TerminalBufferSearchTests, LowerBound_Mask_VeryShortEnd) {
  std::vector<uint8_t> qry = {0b11111100};
  auto it = subrange.lowerBound({qry.data(), 1}, 4);
  ASSERT_EQ(it, buffer.constAt(10));
}

TEST_F(TerminalBufferSearchTests, UpperBound_CommonTerminal_Seen) {
  std::vector<uint8_t> qry = {0b10000000, 0b10101010, 0b00000010};
  auto it = range.upperBound({qry.data(), 3}, k);
  ASSERT_EQ(it, buffer.constAt(8));
}

TEST_F(TerminalBufferSearchTests, UpperBound_CommonTerminal_NotSeen) {
  std::vector<uint8_t> qry = {0b10000011, 0b10101010, 0b00000010};
  auto it = range.upperBound({qry.data(), 3}, k);
  ASSERT_EQ(it, buffer.constAt(8));
}

TEST_F(TerminalBufferSearchTests, UpperBound_Begin) {
  std::vector<uint8_t> qry = {0b00000000, 0b00000000, 0b00000000};
  auto it = range.upperBound({qry.data(), 3}, k);
  ASSERT_EQ(it, buffer.constAt(1));
}

TEST_F(TerminalBufferSearchTests, UpperBound_End) {
  std::vector<uint8_t> qry = {0b11111001, 0b11111111, 0b00000011};
  auto it = range.upperBound({qry.data(), 3}, k);
  ASSERT_EQ(it, buffer.constEnd());
}

TEST_F(TerminalBufferSearchTests, UpperBound_Mask_VeryShort) {
  std::vector<uint8_t> qry = {0b10000000};
  auto it = subrange.upperBound({qry.data(), 1}, 4);
  ASSERT_EQ(it, buffer.constAt(8));
}

TEST_F(TerminalBufferSearchTests, UpperBound_Mask_VeryShort_Truncated) {
  std::vector<uint8_t> qry = {0b00000001};
  auto it = subrange.upperBound({qry.data(), 1}, 3);
  ASSERT_EQ(it, buffer.constAt(8));
}

TEST_F(TerminalBufferSearchTests, UpperBound_Mask_VeryShort_TruncatedSeen) {
  std::vector<uint8_t> qry = {0b00001000};
  auto it = subrange.upperBound({qry.data(), 1}, 3);
  ASSERT_EQ(it, buffer.constAt(10));
}

TEST_F(TerminalBufferSearchTests, EndsWith_Exists_Partial) {
  ShortSuffix sfx(5, 0b0000001010101010u);
  auto res = range.endsWith(sfx);
  ASSERT_EQ(res.constBegin(), buffer.constAt(2));
  ASSERT_EQ(res.constEnd(), buffer.constAt(10));
  ASSERT_EQ(res.size(), 8);
}

TEST_F(TerminalBufferSearchTests, EndsWith_Exists_Full) {
  ShortSuffix sfx(8, 0b0000001010101010100010u);
  auto res = range.endsWith(sfx);
  ASSERT_EQ(res.constBegin(), buffer.constAt(8));
  ASSERT_EQ(res.constEnd(), buffer.constAt(10));
  ASSERT_EQ(res.size(), 2);
}

TEST_F(TerminalBufferSearchTests, EndsWith_DoesntExists) {
  ShortSuffix sfx(2, 0b1011);
  auto res = range.endsWith(sfx);
  ASSERT_TRUE(res.empty());
}

TEST_F(TerminalBufferSearchTests, Retrieve_Exists) {
  ShortSuffix sfx(6, 0b101010100000u);
  auto res = range.retrieve(sfx);
  ASSERT_EQ(res.constBegin(), buffer.constAt(1));
  ASSERT_EQ(res.constEnd(), buffer.constAt(2));
  ASSERT_EQ(res.size(), 1);
}

TEST_F(TerminalBufferSearchTests, Retrieve_DoesntExists_Key) {
  ShortSuffix sfx(6, 0b101010100011u);
  auto res = range.retrieve(sfx);
  ASSERT_TRUE(res.empty());
}

TEST_F(TerminalBufferSearchTests, Retrieve_DoesntExists_Size) {
  ShortSuffix sfx(5, 0b1010101000u);
  auto res = range.retrieve(sfx);
  ASSERT_TRUE(res.empty());
}

TEST_F(TerminalBufferSearchTests, Retrieve_Null) {
  ShortSuffix sfx(0);
  auto res = range.retrieve(sfx);
  ASSERT_TRUE(res.empty());
}

class NullTerminalSearchTests : public testing::Test {
protected:
  NullTerminalSearchTests()
      : k(9), num_genomes(2),
        buffer(k,
               {{/* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000,
                 0b00000000, /* edge */ 0b00000000},
                {/* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000,
                 0b00000000, /* edge */ 0b00000001},
                {/* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000,
                 0b00000000, /* edge */ 0b00000001},
                {/* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000,
                 0b00000000, /* edge */ 0b00000001},
                {/* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000,
                 0b00000000, /* edge */ 0b00000011},
                {/* size */ 0b00000001, /* kmer */ 0b00000000, 0b00000000,
                 0b00000000, /* edge */ 0b00000000},
                {/* size */ 0b00000001, /* kmer */ 0b00000000, 0b00000000,
                 0b00000000, /* edge */ 0b00000001},
                {/* size */ 0b00000001, /* kmer */ 0b00000000, 0b00000000,
                 0b00000000, /* edge */ 0b00000011},
                {/* size */ 0b00000110, /* kmer */ 0b00000000, 0b10101000,
                 0b00000000, /* edge */ 0b00000000},
                {/* size */ 0b00000110, /* kmer */ 0b10000000, 0b10101010,
                 0b00000000, /* edge */ 0b00000000},
                {/* size */ 0b00000110, /* kmer */ 0b10000000, 0b10101010,
                 0b00000000, /* edge */ 0b00000001},
                {/* size */ 0b00000110, /* kmer */ 0b10000000, 0b10101010,
                 0b00000000, /* edge */ 0b00000010},
                {/* size */ 0b00001000, /* kmer */ 0b11111000, 0b11111111,
                 0b00000011, /* edge */ 0b00000111}},
               TerminalBuffer::autofit_tag),
        range(buffer.asRange()),
        subrange(buffer.getEffK(), buffer.lengthBytes(), buffer.constAt(2),
                 buffer.constAt(10)) {}

  uint8_t k;
  std::size_t num_genomes;
  TerminalBuffer buffer;
  TerminalRange range, subrange;
};

TEST_F(NullTerminalSearchTests, RetrieveNull) {
  auto res = range.retrieve(ShortSuffix(0));
  ASSERT_EQ(res.constBegin(), buffer.constAt(0));
  ASSERT_EQ(res.constEnd(), buffer.constAt(5));
}

TEST_F(NullTerminalSearchTests, RetrieveA) {
  auto res = range.retrieve(ShortSuffix(1u, 0b00u));
  ASSERT_EQ(res.constBegin(), buffer.constAt(5));
  ASSERT_EQ(res.constEnd(), buffer.constAt(8));
}

TEST_F(NullTerminalSearchTests, EndsWith_A) {
  ShortSuffix sfx(1, 0b00);
  auto res = range.endsWith(sfx);
  ASSERT_EQ(res.constBegin(), buffer.constAt(5));
  ASSERT_EQ(res.constEnd(), buffer.constAt(12));
}
