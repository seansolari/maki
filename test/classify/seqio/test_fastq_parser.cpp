#include <algorithm>
#include <filesystem>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/utility/views/zip.hpp>

#include "maki/classify/io/fastq.hpp"
#include "maki/maki.h"

namespace fs = std::filesystem;
using namespace seqan3::literals;
using ::testing::ContainerEq;

class FilePairingTests : public testing::Test {
protected:
  FilePairingTests()
      : paths({STRING(GCF_000005845_FWD), STRING(GCF_000005845_REV),
               STRING(GCF_000006765_FWD), STRING(GCF_000006765_REV),
               STRING(GCF_000006925_FWD), STRING(GCF_000006925_REV)}) {}

  std::vector<fs::path> paths;
};

TEST_F(FilePairingTests, PairFastQ) {
  auto pairs = pairInputFiles(paths);
  EXPECT_EQ(pairs.size(), 3);
  EXPECT_EQ(pairs[0],
            DataFilePair(STRING(GCF_000005845_FWD), STRING(GCF_000005845_REV)));
  EXPECT_EQ(pairs[1],
            DataFilePair(STRING(GCF_000006765_FWD), STRING(GCF_000006765_REV)));
  EXPECT_EQ(pairs[2],
            DataFilePair(STRING(GCF_000006925_FWD), STRING(GCF_000006925_REV)));
}

class CharBufferTests : public testing::Test {
protected:
  CharBufferTests() : bufferSize(30), buffer(bufferSize) {}

  const char data[13] = "Hello World!";
  size_t bufferSize;
  detail::CharBuffer buffer;
};

TEST_F(CharBufferTests, EmptyBuffer) {
  EXPECT_EQ(buffer.size(), 0);
  EXPECT_EQ(buffer.capacity(), bufferSize);
  EXPECT_EQ(buffer.bytesRemaining(), bufferSize);
  ASSERT_TRUE(buffer.empty());
}

TEST_F(CharBufferTests, Write) {
  buffer.write(data, 13);
  ASSERT_TRUE(std::equal(buffer.begin(), buffer.back(), data));

  EXPECT_EQ(buffer.size(), 13);
  EXPECT_EQ(buffer.capacity(), bufferSize);
  EXPECT_EQ(buffer.bytesRemaining(), 17);
  ASSERT_FALSE(buffer.empty());
}

TEST_F(CharBufferTests, Flush) {
  buffer.write(data, 13);

  detail::CharBuffer other(bufferSize);
  buffer.flushTo(other);

  EXPECT_EQ(buffer.size(), 0);
  EXPECT_EQ(buffer.capacity(), bufferSize);
  EXPECT_EQ(buffer.bytesRemaining(), bufferSize);
  ASSERT_TRUE(buffer.empty());

  ASSERT_TRUE(std::equal(other.begin(), other.back(), data));
  EXPECT_EQ(other.size(), 13);
  EXPECT_EQ(other.capacity(), bufferSize);
  EXPECT_EQ(other.bytesRemaining(), 17);
  ASSERT_FALSE(other.empty());
}

TEST_F(CharBufferTests, FlushOverflow) {
  buffer.write(data, 13);

  detail::CharBuffer other(bufferSize);
  buffer.flushOverflow(buffer.begin() + 5, other);

  ASSERT_TRUE(std::equal(buffer.begin(), buffer.back(), "Hello"));
  EXPECT_EQ(buffer.size(), 5);
  EXPECT_EQ(buffer.capacity(), bufferSize);
  EXPECT_EQ(buffer.bytesRemaining(), 25);
  ASSERT_FALSE(buffer.empty());

  ASSERT_TRUE(std::equal(other.begin(), other.back(), " World!"));
  EXPECT_EQ(other.size(), 8);
  EXPECT_EQ(other.capacity(), bufferSize);
  EXPECT_EQ(other.bytesRemaining(), 22);
  ASSERT_FALSE(other.empty());
}

TEST_F(CharBufferTests, Swap) {
  buffer.write(data, 13);

  detail::CharBuffer other(bufferSize);
  buffer.flushTo(other);
  other.swap(buffer);

  EXPECT_EQ(other.size(), 0);
  EXPECT_EQ(other.capacity(), bufferSize);
  EXPECT_EQ(other.bytesRemaining(), bufferSize);
  ASSERT_TRUE(other.empty());

  ASSERT_TRUE(std::equal(buffer.begin(), buffer.back(), data));
  EXPECT_EQ(buffer.size(), 13);
  EXPECT_EQ(buffer.capacity(), bufferSize);
  EXPECT_EQ(buffer.bytesRemaining(), 17);
  ASSERT_FALSE(buffer.empty());
}

TEST_F(CharBufferTests, Count) {
  buffer.write(data, 13);
  EXPECT_EQ(buffer.count('o'), 2);
  EXPECT_EQ(buffer.count('\0'), 1);
}

TEST_F(CharBufferTests, CountTo) {
  buffer.write(data, 13);
  EXPECT_EQ(buffer.countTo(buffer.begin() + 7, 'o'), 1);
}

TEST_F(CharBufferTests, Prev) {
  buffer.write(data, 13);

  auto *p = buffer.prev('W');
  EXPECT_EQ(buffer[6], 'W');
  EXPECT_EQ(p, buffer.begin() + 6);
}

TEST_F(CharBufferTests, PrevStart) {
  buffer.write(data, 13);

  auto *p = buffer.prev('H');
  EXPECT_EQ(p, buffer.begin());
}

TEST_F(CharBufferTests, PrevEnd) {
  buffer.write(data, 13);

  auto *p = buffer.prev('\0');
  EXPECT_EQ(p, buffer.back() - 1);
}

TEST_F(CharBufferTests, PrevDNE) {
  buffer.write(data, 13);

  auto *p = buffer.prev('w');
  EXPECT_EQ(p, buffer.back());
}

TEST_F(CharBufferTests, NPrev) {
  buffer.write(data, 13);

  auto *p = buffer.nPrevFrom(buffer.begin() + 7, 1, 'o');
  EXPECT_EQ(p, buffer.begin() + 4);
}

TEST_F(CharBufferTests, NPrevDNE) {
  buffer.write(data, 13);

  auto *p = buffer.nPrevFrom(buffer.begin() + 7, 2, 'o');
  EXPECT_EQ(p, buffer.back());
}

TEST_F(CharBufferTests, NPrevNotPossible) {
  buffer.write(data, 13);

  auto *p = buffer.nPrevFrom(buffer.begin(), 100, 'o');
  EXPECT_EQ(p, buffer.back());
}

class ReadTests : public testing::Test {
protected:
  ReadTests() : read("test_read") {}
  Read read;
};

TEST_F(ReadTests, ReadId) { EXPECT_EQ(read.id, std::string("test_read")); }

TEST_F(ReadTests, Push) {
  seqan3::dna5_vector data = {
      'A'_dna5, 'C'_dna5, 'C'_dna5, 'A'_dna5, 'G'_dna5, 'T'_dna5, 'A'_dna5,
      'G'_dna5, 'T'_dna5, 'N'_dna5, 'N'_dna5, 'G'_dna5, 'C'_dna5, 'N'_dna5,
      'C'_dna5, 'A'_dna5, 'C'_dna5, 'T'_dna5, 'T'_dna5, 'N'_dna5, 'G'_dna5,
      'A'_dna5, 'C'_dna5, 'A'_dna5, 'G'_dna5, 'N'_dna5, 'T'_dna5, 'A'_dna5,
      'C'_dna5, 'G'_dna5, 'T'_dna5, 'T'_dna5, 'G'_dna5, 'T'_dna5, 'G'_dna5,
      'T'_dna5, 'G'_dna5, 'N'_dna5, 'G'_dna5, 'N'_dna5};

  read.push(std::move(data), 2);

  EXPECT_EQ(read.sequences[0].size(), 5);
  EXPECT_EQ(read.sequences[1].size(), 5);

  std::vector expectedRaw = {Dna4Sequence("ACCAGTAGT"_dna4),
                             Dna4Sequence("GC"_dna4),
                             Dna4Sequence("CACTT"_dna4),
                             Dna4Sequence("GACAG"_dna4),
                             Dna4Sequence("TACGTTGTGTG"_dna4)},
              expectedRcomp = {
                  Dna4Sequence("ACTACTGGT"_dna4), Dna4Sequence("GC"_dna4),
                  Dna4Sequence("AAGTG"_dna4), Dna4Sequence("CTGTC"_dna4),
                  Dna4Sequence("CACACAACGTA"_dna4)};

  ASSERT_THAT(read.sequences[0], ContainerEq(expectedRaw));
  ASSERT_THAT(read.sequences[1], ContainerEq(expectedRcomp));

  ASSERT_EQ(read.numFragments(), 10);
  ASSERT_EQ(detail::countBps(ReadVector{read}), 64);
  ASSERT_EQ(detail::numFragments(ReadVector{read}), 10);
}

class ReadFakeChunkTests : public testing::Test {
protected:
  ReadFakeChunkTests()
      : fwdData("@fwd1 test read\nACGACNTAGCTGACGTNACGTAGC\n@fwd2 test "
                "read\nAGTCGACTTCAGTNCATCGGCTA\n"),
        revData(
            "@rev1 test read\nAGTCAGTCGACGTACTNNAGCTNAGCTAGTCAGTCGTAC\n@rev2 "
            "test read\nATATANTATTAGCTANATG\n"),
        ifwd(fwdData), irev(revData), buffersize(64), buffer(buffersize),
        overflow(buffersize) {}

  std::string fwdData, revData;
  std::stringstream ifwd, irev;
  size_t buffersize;
  detail::SharedBufferPair buffer;
  detail::BufferPair overflow;
};

TEST_F(ReadFakeChunkTests, Read) {
  detail::Operator_ReadChunk Op(&ifwd, &irev, '@', &buffer, 1, &overflow);
  detail::SharedBufferPair *p = Op.readChunk();

  EXPECT_EQ(p, &buffer);

  // parsed data

  ASSERT_TRUE(std::equal(buffer.forward().begin(), buffer.forward().back(),
                         "@fwd1 test read\nACGACNTAGCTGACGTNACGTAGC\n"));
  EXPECT_EQ(buffer.forward().size(), 41);
  ASSERT_TRUE(
      std::equal(buffer.reverse().begin(), buffer.reverse().back(),
                 "@rev1 test read\nAGTCAGTCGACGTACTNNAGCTNAGCTAGTCAGTCGTAC\n"));
  EXPECT_EQ(buffer.reverse().size(), 56);

  // overflow

  ASSERT_TRUE(std::equal(overflow.forward().begin(), overflow.forward().back(),
                         "@fwd2 test read\nAGTCGAC"));
  EXPECT_EQ(overflow.forward().size(), 23);
  ASSERT_TRUE(std::equal(overflow.reverse().begin(), overflow.reverse().back(),
                         "@rev2 te"));
  EXPECT_EQ(overflow.reverse().size(), 8);

  // parse remaining

  buffer.clear();
  buffer.setAvail();

  Op.readChunk();

  // stream state

  ASSERT_TRUE(ifwd.eof() && irev.eof());
  ASSERT_TRUE(Op.completedRec(&ifwd));
  ASSERT_TRUE(Op.completedRec(&irev));
  ASSERT_TRUE(overflow.empty());

  // parsed data

  ASSERT_TRUE(std::equal(buffer.forward().begin(), buffer.forward().back(),
                         "@fwd2 test read\nAGTCGACTTCAGTNCATCGGCTA\n"));
  EXPECT_EQ(buffer.forward().size(), 40);
  ASSERT_TRUE(std::equal(buffer.reverse().begin(), buffer.reverse().back(),
                         "@rev2 test read\nATATANTATTAGCTANATG\n"));
  EXPECT_EQ(buffer.reverse().size(), 36);
}

class ReadParsingTests : public testing::Test {
protected:
  ReadParsingTests()
      : rp(STRING(GCF_000005845_FWD), STRING(GCF_000005845_REV)),
        minReadLength(1), minFragLength(1) {}

  ReadVector parseSerial() const {
    ReadVector out;

    seqan3::sequence_file_input ifwd{rp.forwardFile}, irev{rp.reverseFile};

    for (auto &&[rec1, rec2] : seqan3::views::zip(ifwd, irev)) {
      if ((rec1.sequence().size() >= minReadLength) &&
          (rec2.sequence().size() >= minReadLength)) {
        Read &v = out.emplace_back(std::move(rec1.id()));
        v.push(std::move(rec1).sequence(), minFragLength);
        v.push(std::move(rec2).sequence(), minFragLength);
      }
    }

    return out;
  }

  void sort(ReadVector &v) const {
    std::sort(v.begin(), v.end(), [](Read const &lhs, Read const &rhs) -> bool {
      return lhs.id < rhs.id;
    });
  }

  DataFilePair rp;
  size_t minReadLength, minFragLength;
};

TEST_F(ReadParsingTests, Parse) {
  ReadVector raw = flattenReads(
                 detail::parsePairedFastq(rp, minReadLength, minFragLength, 2)),
             actual = parseSerial();

  sort(raw);
  sort(actual);

  ASSERT_THAT(raw, ContainerEq(actual));
}

class ReadChunkingTests : public testing::Test {
protected:
  ReadChunkingTests()
      : rp(STRING(GCF_000005845_SMFWD), STRING(GCF_000005845_SMREV)), k(3),
        threads(3) {}

  DataFilePair rp;
  uint8_t k;
  size_t threads;
};

TEST_F(ReadChunkingTests, Chunks) {
  auto chunks = detail::parsePairedFastq(rp, k, k, threads, 400);
  ASSERT_GE(chunks.data.size(), threads);
  ASSERT_EQ(chunks.numFragments(), 20);
  ASSERT_EQ(chunks.numTerminals(k), k * 20);
}

TEST_F(ReadChunkingTests, Counting) {
  auto chunks = detail::parsePairedFastq(rp, k, k, threads, 400);
  ASSERT_EQ(chunks.numFragments(), 5 * 2 * 2);
  ASSERT_EQ(chunks.countBps(), 5 * 2 * 2 * 150);
}