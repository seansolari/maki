
#include "maki/build/graph/build_colours.hpp"
#include "maki/build/io/filter.hpp"
#include "maki/core/utils/tempfile.hpp"
#include "maki/maki.h"
#include "test_common.hpp"
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using ::testing::ContainerEq;

const size_t minContigSize = 3;

TEST(ParseFasta, ParseBitpackedFasta) {
  Colours c;
  auto genome = parseGenome(STRING(SEQA_FNA), c, minContigSize);

  EXPECT_EQ(genome.length(), 80);
  EXPECT_EQ(genome.contigs.size(), 2);
  EXPECT_EQ(genome.numFragments(), 10);

  const auto &contig1 = genome.contigs[0];
  EXPECT_EQ(contig1.accn, "genome1");
  EXPECT_EQ(contig1.minFragmentSize, minContigSize);
  EXPECT_EQ(contig1.totalLength, 40);
  EXPECT_EQ(contig1.numAnnotations, 0);
  EXPECT_EQ(contig1.numFragments(), 6);

  EXPECT_EQ(contig1.sequences[0].size(), 3);
  EXPECT_EQ(contig1.sequences[0][0].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences[0][0].sequence, Dna4Sequence("ACGAG"_dna4));
  EXPECT_EQ(contig1.sequences[0][0].annotations.size(), 0);
  EXPECT_EQ(contig1.sequences[0][1].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences[0][1].sequence,
            Dna4Sequence("CGGTGGCAA"_dna4));
  EXPECT_EQ(contig1.sequences[0][1].annotations.size(), 0);
  EXPECT_EQ(contig1.sequences[0][2].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences[0][2].sequence, Dna4Sequence("GAAGTT"_dna4));
  EXPECT_EQ(contig1.sequences[0][2].annotations.size(), 0);

  EXPECT_EQ(contig1.sequences[1].size(), 3);
  EXPECT_EQ(contig1.sequences[1][0].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences[1][0].sequence, Dna4Sequence("CTCGT"_dna4));
  EXPECT_EQ(contig1.sequences[1][0].annotations.size(), 0);
  EXPECT_EQ(contig1.sequences[1][1].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences[1][1].sequence,
            Dna4Sequence("TTGCCACCG"_dna4));
  EXPECT_EQ(contig1.sequences[1][1].annotations.size(), 0);
  EXPECT_EQ(contig1.sequences[1][2].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences[1][2].sequence, Dna4Sequence("AACTTC"_dna4));
  EXPECT_EQ(contig1.sequences[1][2].annotations.size(), 0);

  const auto &contig2 = genome.contigs[1];
  EXPECT_EQ(contig2.accn, "genome2");
  EXPECT_EQ(contig2.minFragmentSize, minContigSize);
  EXPECT_EQ(contig2.totalLength, 40);
  EXPECT_EQ(contig2.numAnnotations, 0);
  EXPECT_EQ(contig2.numFragments(), 4);

  EXPECT_EQ(contig2.sequences[0].size(), 2);
  EXPECT_EQ(contig2.sequences[0][0].nullFeatureId, 2);
  EXPECT_EQ(contig2.sequences[0][0].sequence,
            Dna4Sequence("GTGTCGGAGG"_dna4));
  EXPECT_EQ(contig2.sequences[0][0].annotations.size(), 0);
  EXPECT_EQ(contig2.sequences[0][1].nullFeatureId, 2);
  EXPECT_EQ(contig2.sequences[0][1].sequence,
            Dna4Sequence("CTCCATCGAC"_dna4));
  EXPECT_EQ(contig2.sequences[0][1].annotations.size(), 0);

  EXPECT_EQ(contig2.sequences[1].size(), 2);
  EXPECT_EQ(contig2.sequences[1][0].nullFeatureId, 2);
  EXPECT_EQ(contig2.sequences[1][0].sequence,
            Dna4Sequence("CCTCCGACAC"_dna4));
  EXPECT_EQ(contig2.sequences[1][0].annotations.size(), 0);
  EXPECT_EQ(contig2.sequences[1][1].nullFeatureId, 2);
  EXPECT_EQ(contig2.sequences[1][1].sequence,
            Dna4Sequence("GTCGATGGAG"_dna4));
  EXPECT_EQ(contig2.sequences[1][1].annotations.size(), 0);
}

inline std::string MakeTempPath(const char *ext) {
  return tempio::create_temporary_file(std::filesystem::temp_directory_path(),
                                       ext);
}

class ParseChunks : public testing::Test {
protected:
  ParseChunks() : fna(MakeTempPath(".fna")), c() {
    std::string data =
        ">seq1\nACGTACGTACGTACGTACGTACGTACGATCAGTCAGTCAGTCGTA\n>"
        "seq2\nAGTACGTCGTACGATCGTC\n>seq3\nATCGTCGATGCTAGCTAGCTAGCTAGCTACGT\n>"
        "seq4\nGTACGTGCTAGCTAGCTGACTCGATGCATTAA\n>seq5\nTAGTATATATAGTAGTAGT\n";
    writeToFna(data, fna);
  }

  ChunkedDna4Genome parse() {
    auto genome = parseGenome(fna, c, kmer_size);
    return chunkFNA(std::move(genome), 16, kmer_size);
  }

  ~ParseChunks() { std::filesystem::remove(fna); }

  std::string fna;
  uint8_t kmer_size = 3;
  size_t threads = 2;
  Colours c;
};

TEST_F(ParseChunks, CheckStructure) {
  auto xgenome = parse();

  // check IDs
  ASSERT_EQ(c.getOrAssign("seq1"), 1);
  ASSERT_EQ(c.getOrAssign("seq2"), 2);
  ASSERT_EQ(c.getOrAssign("seq3"), 3);
  ASSERT_EQ(c.getOrAssign("seq4"), 4);
  ASSERT_EQ(c.getOrAssign("seq5"), 5);

  // check sequence metadata
  ASSERT_EQ(xgenome.size(), 18);

  ASSERT_FALSE(xgenome[0].endIsTerminal());
  ASSERT_EQ(xgenome[0].numTerminals(kmer_size), kmer_size);
  ASSERT_EQ(xgenome[0].id(), 1);
  ASSERT_FALSE(xgenome[1].endIsTerminal());
  ASSERT_EQ(xgenome[1].numTerminals(kmer_size), 0);
  ASSERT_EQ(xgenome[1].id(), 1);
  ASSERT_TRUE(xgenome[2].endIsTerminal());
  ASSERT_EQ(xgenome[2].numTerminals(kmer_size), 0);
  ASSERT_EQ(xgenome[2].id(), 1);
  ASSERT_FALSE(xgenome[3].endIsTerminal());
  ASSERT_EQ(xgenome[3].numTerminals(kmer_size), kmer_size);
  ASSERT_EQ(xgenome[3].id(), 1);
  ASSERT_FALSE(xgenome[4].endIsTerminal());
  ASSERT_EQ(xgenome[4].numTerminals(kmer_size), 0);
  ASSERT_EQ(xgenome[4].id(), 1);
  ASSERT_TRUE(xgenome[5].endIsTerminal());
  ASSERT_EQ(xgenome[5].numTerminals(kmer_size), 0);
  ASSERT_EQ(xgenome[5].id(), 1);
  ASSERT_TRUE(xgenome[6].endIsTerminal());
  ASSERT_EQ(xgenome[6].numTerminals(kmer_size), kmer_size);
  ASSERT_EQ(xgenome[6].id(), 2);
  ASSERT_TRUE(xgenome[7].endIsTerminal());
  ASSERT_EQ(xgenome[7].numTerminals(kmer_size), kmer_size);
  ASSERT_EQ(xgenome[7].id(), 2);
  ASSERT_FALSE(xgenome[8].endIsTerminal());
  ASSERT_EQ(xgenome[8].numTerminals(kmer_size), kmer_size);
  ASSERT_EQ(xgenome[8].id(), 3);
  ASSERT_TRUE(xgenome[9].endIsTerminal());
  ASSERT_EQ(xgenome[9].numTerminals(kmer_size), 0);
  ASSERT_EQ(xgenome[9].id(), 3);
  ASSERT_FALSE(xgenome[10].endIsTerminal());
  ASSERT_EQ(xgenome[10].numTerminals(kmer_size), kmer_size);
  ASSERT_EQ(xgenome[10].id(), 3);
  ASSERT_TRUE(xgenome[11].endIsTerminal());
  ASSERT_EQ(xgenome[11].numTerminals(kmer_size), 0);
  ASSERT_EQ(xgenome[11].id(), 3);
  ASSERT_FALSE(xgenome[12].endIsTerminal());
  ASSERT_EQ(xgenome[12].numTerminals(kmer_size), kmer_size);
  ASSERT_EQ(xgenome[12].id(), 4);
  ASSERT_TRUE(xgenome[13].endIsTerminal());
  ASSERT_EQ(xgenome[13].numTerminals(kmer_size), 0);
  ASSERT_EQ(xgenome[13].id(), 4);
  ASSERT_FALSE(xgenome[14].endIsTerminal());
  ASSERT_EQ(xgenome[14].numTerminals(kmer_size), kmer_size);
  ASSERT_EQ(xgenome[14].id(), 4);
  ASSERT_TRUE(xgenome[15].endIsTerminal());
  ASSERT_EQ(xgenome[15].numTerminals(kmer_size), 0);
  ASSERT_EQ(xgenome[15].id(), 4);
  ASSERT_TRUE(xgenome[16].endIsTerminal());
  ASSERT_EQ(xgenome[16].numTerminals(kmer_size), kmer_size);
  ASSERT_EQ(xgenome[16].id(), 5);
  ASSERT_TRUE(xgenome[17].endIsTerminal());
  ASSERT_EQ(xgenome[17].numTerminals(kmer_size), kmer_size);
  ASSERT_EQ(xgenome[17].id(), 5);
}

TEST_F(ParseChunks, NumKmers) {
  auto xgenome = parse();

  // predicted sequence lengths
  std::vector<std::size_t> actual;
  for (std::size_t i = 0; i < xgenome.size(); ++i) {
    actual.push_back(xgenome[i].numKmers(kmer_size));
  }

  // check
  std::vector<std::size_t> expected = {17, 17, 11, 17, 17, 11, 17, 17, 17,
                                       14, 17, 14, 17, 14, 17, 14, 17, 17};
  ASSERT_THAT(actual, ContainerEq(expected));
}

TEST_F(ParseChunks, NumTerminals) {
  auto xgenome = parse();

  // predicted sequence lengths
  std::vector<std::size_t> actual;
  for (std::size_t i = 0; i < xgenome.size(); ++i) {
    actual.push_back(xgenome[i].numTerminals(kmer_size));
  }

  // check
  std::vector<std::size_t> expected = {
      kmer_size, 0,         0,         kmer_size, 0,         0,
      kmer_size, kmer_size, kmer_size, 0,         kmer_size, 0,
      kmer_size, 0,         kmer_size, 0,         kmer_size, kmer_size};
  ASSERT_THAT(actual, ContainerEq(expected));
}

TEST_F(ParseChunks, FragmentView) {
  auto xgenome = parse();

  // output sequences
  std::vector<Dna4Sequence> sequenceChunks;
  for (std::size_t i = 0; i < xgenome.size(); ++i) {
    sequenceChunks.emplace_back(xgenome[i].fragments(kmer_size));
  }

  // check

  std::vector<Dna4Sequence> expectedChunks = {
      Dna4Sequence("ACGTACGTACGTACGTACG"_dna4),
      Dna4Sequence("ACGTACGTACGATCAGTCA"_dna4),
      Dna4Sequence("TCAGTCAGTCGTA"_dna4),
      Dna4Sequence("TACGACTGACTGACTGATC"_dna4),
      Dna4Sequence("ATCGTACGTACGTACGTAC"_dna4),
      Dna4Sequence("TACGTACGTACGT"_dna4),
      Dna4Sequence("AGTACGTCGTACGATCGTC"_dna4),
      Dna4Sequence("GACGATCGTACGACGTACT"_dna4),
      Dna4Sequence("ATCGTCGATGCTAGCTAGC"_dna4),
      Dna4Sequence("AGCTAGCTAGCTACGT"_dna4),
      Dna4Sequence("ACGTAGCTAGCTAGCTAGC"_dna4),
      Dna4Sequence("AGCTAGCATCGACGAT"_dna4),
      Dna4Sequence("GTACGTGCTAGCTAGCTGA"_dna4),
      Dna4Sequence("TGACTCGATGCATTAA"_dna4),
      Dna4Sequence("TTAATGCATCGAGTCAGCT"_dna4),
      Dna4Sequence("GCTAGCTAGCACGTAC"_dna4),
      Dna4Sequence("TAGTATATATAGTAGTAGT"_dna4),
      Dna4Sequence("ACTACTACTATATATACTA"_dna4)};

  ASSERT_THAT(sequenceChunks, ContainerEq(expectedChunks))
      << "actual: " << toString(sequenceChunks, ", ")
      << "\nexpected: " << toString(expectedChunks, ", ");
}

TEST_F(ParseChunks, TerminalView) {
  auto xgenome = parse();

  // output sequences
  std::vector<Dna4Sequence> sequenceChunks;
  for (std::size_t i = 0; i < xgenome.size(); ++i) {
    if (xgenome[i].numTerminals(kmer_size) > 0)
      sequenceChunks.emplace_back(xgenome[i].terminals());
  }

  // check
  std::vector<Dna4Sequence> expectedChunks = {
      Dna4Sequence("ACGTACGTACGTACGTACG"_dna4),
      Dna4Sequence("TACGACTGACTGACTGATC"_dna4),
      Dna4Sequence("AGTACGTCGTACGATCGTC"_dna4),
      Dna4Sequence("GACGATCGTACGACGTACT"_dna4),
      Dna4Sequence("ATCGTCGATGCTAGCTAGC"_dna4),
      Dna4Sequence("ACGTAGCTAGCTAGCTAGC"_dna4),
      Dna4Sequence("GTACGTGCTAGCTAGCTGA"_dna4),
      Dna4Sequence("TTAATGCATCGAGTCAGCT"_dna4),
      Dna4Sequence("TAGTATATATAGTAGTAGT"_dna4),
      Dna4Sequence("ACTACTACTATATATACTA"_dna4)};

  ASSERT_THAT(sequenceChunks, ContainerEq(expectedChunks))
      << "actual: " << toString(sequenceChunks, ", ")
      << "\nexpected: " << toString(expectedChunks, ", ");
}
