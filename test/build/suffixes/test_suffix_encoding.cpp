
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/core/seq/seq_io.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using ::testing::ElementsAreArray;
using namespace seqan3::literals;

// encoding ------------------------------------------------------------------------------------

class WriteKmers : public testing::Test {
protected:
    WriteKmers() :
        _seq("ACGTACGTACGTACGT"_dna4)
    {}

    Dna4Sequence _seq;
};

TEST_F(WriteKmers, WriteSeven) {
    uint8_t arr[2] = { 0, 0 };
    writeMerTo(_seq.cbegin(), arr, 7);
    ASSERT_THAT(arr, ElementsAreArray({ 0b11100100, 0b00100100 }));
}

TEST_F(WriteKmers, WriteSixteen) {
    uint8_t arr[4] = { 0, 0, 0, 0 };
    writeMerTo(_seq.cbegin(), arr, 16);
    ASSERT_THAT(arr, ElementsAreArray({ 0b11100100, 0b11100100, 0b11100100, 0b11100100 }));
}

class SerDeserTests : public testing::Test {
protected:
    SerDeserTests() {}
};

TEST_F(SerDeserTests, FullSer) {
    size_t in = (size_t)-1;
    uint8_t arr[sizeof(size_t)];
    ser(in, arr, sizeof(size_t));
    size_t out = deser(arr, sizeof(size_t));
    ASSERT_EQ(in, out);
}

TEST_F(SerDeserTests, PartialSer) {
    size_t in = 0xFF10FF;
    uint8_t arr[3];
    ser(in, arr, 3);
    size_t out = deser(arr, 3);
    ASSERT_EQ(in, out);
}

// rolling sequences ---------------------------------------------------------------------------

class RollingSequences : public testing::Test {
protected:
    RollingSequences() :
        _seq("ACGTACGTACGTACGT"_dna4)
    {}
    
    Dna4Sequence _seq;
};

// ShortSuffix

TEST_F(RollingSequences, Small_SizeItCons) {
    ShortSuffix x(4, _seq.cbegin());
    EXPECT_EQ(x.size(), 4);
    EXPECT_EQ(x.data(), 0b11100100);
}

TEST_F(RollingSequences, Small_RangeCons) {
    ShortSuffix x(_seq.cbegin(), _seq.cend());
    EXPECT_EQ(x.size(), 16);
    EXPECT_EQ(x.data(), 0b11100100111001001110010011100100);
}

TEST_F(RollingSequences, Small_RollBack) {
    ShortSuffix x(4, _seq.cbegin());
    x.roll(0b00);
    x.roll(0b11);
    EXPECT_EQ(x.size(), 4);
    EXPECT_EQ(x.data(), 0b11001110);
}

TEST_F(RollingSequences, Small_OperatorPlus) {
    ShortSuffix x(4, _seq.cbegin());
    x = x + 0b00;
    x = x + 0b11;
    EXPECT_EQ(x.size(), 6);
    EXPECT_EQ(x.data(), 0b111001000011);
}

TEST_F(RollingSequences, Small_OperatorLeftBitShift) {
    ShortSuffix x(4, _seq.cbegin());
    auto y = x << 2;
    EXPECT_EQ(y.size(), 6);
    EXPECT_EQ(y.data(), 0b111001000000);
}

TEST_F(RollingSequences, Small_ToString) {
    ShortSuffix x(4, _seq.cbegin());
    EXPECT_EQ(x.toString(), std::string { "ACGT" });
}

TEST_F(RollingSequences, Small_MSB) {
    ShortSuffix x(4, _seq.cbegin());
    EXPECT_EQ(x.msb(), 0b00000011u);
}

// long rolling sequence

TEST_F(RollingSequences, Long_SizeItCons) {
    Kmer x(15, _seq.cbegin());
    EXPECT_EQ(x.size(), 15);
    ASSERT_THAT(x.view(), ElementsAreArray({ 0b11100100, 0b11100100, 0b11100100, 0b00100100 }));
}

TEST_F(RollingSequences, Long_RangeCons) {
    Kmer x(_seq.cbegin(), _seq.cbegin() + 15);
    EXPECT_EQ(x.size(), 15);
    ASSERT_THAT(x.view(), ElementsAreArray({ 0b11100100, 0b11100100, 0b11100100, 0b00100100 }));
}

TEST_F(RollingSequences, Long_RollBack) {
    Kmer x(15, _seq.cbegin());
    x.roll(0b11);
    EXPECT_EQ(x.size(), 15);
    ASSERT_THAT(x.view(), ElementsAreArray({ 0b00111001, 0b00111001, 0b00111001, 0b00111001 }));
}

TEST_F(RollingSequences, Long_ToString) {
    Kmer x(5, _seq.cbegin() + 5);
    std::string v = x.toString();
    EXPECT_EQ(v, std::string{ "CGTAC" });
}

// rolling suffix ------------------------------------------------------------------

class RollingTerminalTests : public testing::Test {
protected:
    RollingTerminalTests() : _sfx(7) {}

    LongSuffix _sfx;
};

TEST_F(RollingTerminalTests, Initial) {
    EXPECT_EQ(_sfx.terminalLength(), 0);
    ASSERT_THAT(_sfx.view(), ElementsAreArray({ 0, 0 }));
}

TEST_F(RollingTerminalTests, RollBackNotFull) {
    _sfx.push(0b11);
    _sfx.push(0b10);
    _sfx.push(0b01);
    _sfx.push(0b00);
    _sfx.push(0b11);
    _sfx.push(0b10);
    EXPECT_EQ(_sfx.terminalLength(), 6);
    ASSERT_THAT(_sfx.view(), ElementsAreArray({ 0b01101100, 0b00101100 }));
}

TEST_F(RollingTerminalTests, RollBackFull) {
    for (size_t i = 0; i < 7; ++i)
        _sfx.push(0b10);
    EXPECT_EQ(_sfx.terminalLength(), 7);
    ASSERT_THAT(_sfx.view(), ElementsAreArray({ 0b10101010, 0b00101010 }));
}

TEST_F(RollingTerminalTests, RollBackExceed) {
    for (size_t i = 0; i < 14; ++i)
        _sfx.push(0b10);
    EXPECT_EQ(_sfx.terminalLength(), 7);
    ASSERT_THAT(_sfx.view(), ElementsAreArray({ 0b10101010, 0b00101010 }));
}
