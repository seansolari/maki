
#include "maki/build/graph/build_colours.hpp"
#include "maki/build/io/fasta.hpp"
#include "maki/core/seq/concepts.hpp"
#include "maki/core/utils/tempfile.hpp"
#include "maki/maki.h"
#include "test_common.hpp"
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using ::testing::ContainerEq;

Dna4Genome parseGenome(const std::string &fastaFile, Colours &colours,
                       std::size_t k) {
  // base input file stream that reads bytes
  zstr::ifstream zis(fastaFile);
  std::string fastaData(std::istreambuf_iterator<char>(zis), {});

  Dna4Genome obj;
  std::istringstream fs(fastaData);
  parseFastaStream(obj, fs, colours, k);

  return obj;
}

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
  ParseChunks() : fna(MakeTempPath(".fna")), c(), xgenome() {
    std::string data =
        ">seq1\nACGTACGTACGTACGTACGTACGTACGATCAGTCAGTCAGTCGTA\n>"
        "seq2\nAGTACGTCGTACGATCGTC\n>seq3\nATCGTCGATGCTAGCTAGCTAGCTAGCTACGT\n>"
        "seq4\nGTACGTGCTAGCTAGCTGACTCGATGCATTAA\n>seq5\nTAGTATATATAGTAGTAGT\n";
    writeToFna(data, fna);
    xgenome = parseFilterFNA(fna, c, threads, kmer_size);
  }
  ~ParseChunks() { std::filesystem::remove(fna); }

  inline std::vector<const SequenceContainer *> getView() const {
    std::vector<const SequenceContainer *> ptrs;
    for (const auto &chunk : xgenome.chunks) {
      ptrs.push_back(&chunk);
    }
    return ptrs;
  }

  std::string fna;
  uint8_t kmer_size = 3;
  size_t threads = 2;
  Colours c;
  ChunkedDna4Genome xgenome;
};

TEST_F(ParseChunks, CheckStructure) {
  // check IDs
  ASSERT_EQ(c.getOrAssign("seq1"), 1);
  ASSERT_EQ(c.getOrAssign("seq2"), 2);
  ASSERT_EQ(c.getOrAssign("seq3"), 3);
  ASSERT_EQ(c.getOrAssign("seq4"), 4);
  ASSERT_EQ(c.getOrAssign("seq5"), 5);

  ASSERT_EQ(xgenome.genome.contigs[0].sequences[0][0].nullFeatureId, 1);
  ASSERT_EQ(xgenome.genome.contigs[0].sequences[1][0].nullFeatureId, 1);
  ASSERT_EQ(xgenome.genome.contigs[1].sequences[0][0].nullFeatureId, 2);
  ASSERT_EQ(xgenome.genome.contigs[1].sequences[1][0].nullFeatureId, 2);
  ASSERT_EQ(xgenome.genome.contigs[2].sequences[0][0].nullFeatureId, 3);
  ASSERT_EQ(xgenome.genome.contigs[2].sequences[1][0].nullFeatureId, 3);
  ASSERT_EQ(xgenome.genome.contigs[3].sequences[0][0].nullFeatureId, 4);
  ASSERT_EQ(xgenome.genome.contigs[3].sequences[1][0].nullFeatureId, 4);
  ASSERT_EQ(xgenome.genome.contigs[4].sequences[0][0].nullFeatureId, 5);
  ASSERT_EQ(xgenome.genome.contigs[4].sequences[1][0].nullFeatureId, 5);

  // check sequence metadata
  ASSERT_EQ(xgenome.chunks.size(), 18);

  ASSERT_FALSE(xgenome.chunks[0].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[0].numTerminals(kmer_size), kmer_size);
  ASSERT_FALSE(xgenome.chunks[1].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[1].numTerminals(kmer_size), 0);
  ASSERT_TRUE(xgenome.chunks[2].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[2].numTerminals(kmer_size), 0);
  ASSERT_FALSE(xgenome.chunks[3].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[3].numTerminals(kmer_size), kmer_size);
  ASSERT_FALSE(xgenome.chunks[4].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[4].numTerminals(kmer_size), 0);
  ASSERT_TRUE(xgenome.chunks[5].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[5].numTerminals(kmer_size), 0);
  ASSERT_TRUE(xgenome.chunks[6].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[6].numTerminals(kmer_size), kmer_size);
  ASSERT_TRUE(xgenome.chunks[7].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[7].numTerminals(kmer_size), kmer_size);
  ASSERT_FALSE(xgenome.chunks[8].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[8].numTerminals(kmer_size), kmer_size);
  ASSERT_TRUE(xgenome.chunks[9].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[9].numTerminals(kmer_size), 0);
  ASSERT_FALSE(xgenome.chunks[10].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[10].numTerminals(kmer_size), kmer_size);
  ASSERT_TRUE(xgenome.chunks[11].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[11].numTerminals(kmer_size), 0);
  ASSERT_FALSE(xgenome.chunks[12].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[12].numTerminals(kmer_size), kmer_size);
  ASSERT_TRUE(xgenome.chunks[13].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[13].numTerminals(kmer_size), 0);
  ASSERT_FALSE(xgenome.chunks[14].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[14].numTerminals(kmer_size), kmer_size);
  ASSERT_TRUE(xgenome.chunks[15].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[15].numTerminals(kmer_size), 0);
  ASSERT_TRUE(xgenome.chunks[16].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[16].numTerminals(kmer_size), kmer_size);
  ASSERT_TRUE(xgenome.chunks[17].endIsTerminal());
  ASSERT_EQ(xgenome.chunks[17].numTerminals(kmer_size), kmer_size);
}

TEST_F(ParseChunks, NumKmers) {
  auto view = getView();

  // predicted sequence lengths
  std::vector<std::size_t> actual;
  for (auto rng : view) {
    actual.push_back(rng->numKmers(kmer_size));
  }

  // check
  std::vector<std::size_t> expected = {17, 17, 11, 17, 17, 11, 17, 17, 17,
                                       14, 17, 14, 17, 14, 17, 14, 17, 17};
  ASSERT_THAT(actual, ContainerEq(expected));
}

TEST_F(ParseChunks, NumTerminals) {
  auto view = getView();

  // predicted sequence lengths
  std::vector<std::size_t> actual;
  for (auto rng : view) {
    actual.push_back(rng->numTerminals(kmer_size));
  }

  // check
  std::vector<std::size_t> expected = {
      kmer_size, 0,         0,         kmer_size, 0,         0,
      kmer_size, kmer_size, kmer_size, 0,         kmer_size, 0,
      kmer_size, 0,         kmer_size, 0,         kmer_size, kmer_size};
  ASSERT_THAT(actual, ContainerEq(expected));
}

TEST_F(ParseChunks, FragmentView) {
  auto view = getView();

  // output sequences
  std::vector<Dna4Sequence> sequenceChunks;
  for (auto rng : view) {
    for (auto seq : rng->fragments(kmer_size)) {
      sequenceChunks.emplace_back(seq.begin(), seq.end());
    }
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
  auto view = getView();

  // output sequences
  std::vector<Dna4Sequence> sequenceChunks;
  for (auto rng : view) {
    for (auto seq : rng->terminals()) {
      sequenceChunks.emplace_back(seq.begin(), seq.end());
    }
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
