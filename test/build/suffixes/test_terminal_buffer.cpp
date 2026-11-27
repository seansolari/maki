
#include "maki/build/io/fasta.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/build/kmers/construct_terminals.hpp"
#include "maki/build/io/fastq.hpp"
#include "test_common.hpp"
#include <cstdint>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>

using ::testing::ElementsAreArray;

using namespace seqan3::literals;

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
                       t2 = {01, 0b00000000, 0b00000000, 0b00001100};
  auto res = TerminalDiff{1, 10}(t1.data(), t2.data());
  ASSERT_EQ(res, BW_0_K);
}

TEST_F(TerminalDiffTests, BigTermSeeminglyKm1) {
  std::vector<uint8_t> t1 = {3, 0b01100000, 0b00000011},
                       t2 = {5, 0b01100001, 0b00000011};
  auto res = TerminalDiff{1, 5}(t1.data(), t2.data());
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
