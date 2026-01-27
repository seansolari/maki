
#include "maki/build/graph/build_colours.hpp"
#include "maki/build/io/fasta.hpp"
#include "maki/core/utils/tempfile.hpp"
#include "maki/maki.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

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

  EXPECT_EQ(contig1.sequences.first.size(), 3);
  EXPECT_EQ(contig1.sequences.first[0].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences.first[0].sequence, Dna4Sequence("ACGAG"_dna4));
  EXPECT_EQ(contig1.sequences.first[0].annotations.size(), 0);
  EXPECT_EQ(contig1.sequences.first[1].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences.first[1].sequence,
            Dna4Sequence("CGGTGGCAA"_dna4));
  EXPECT_EQ(contig1.sequences.first[1].annotations.size(), 0);
  EXPECT_EQ(contig1.sequences.first[2].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences.first[2].sequence, Dna4Sequence("GAAGTT"_dna4));
  EXPECT_EQ(contig1.sequences.first[2].annotations.size(), 0);

  EXPECT_EQ(contig1.sequences.second.size(), 3);
  EXPECT_EQ(contig1.sequences.second[0].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences.second[0].sequence, Dna4Sequence("CTCGT"_dna4));
  EXPECT_EQ(contig1.sequences.second[0].annotations.size(), 0);
  EXPECT_EQ(contig1.sequences.second[1].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences.second[1].sequence,
            Dna4Sequence("TTGCCACCG"_dna4));
  EXPECT_EQ(contig1.sequences.second[1].annotations.size(), 0);
  EXPECT_EQ(contig1.sequences.second[2].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences.second[2].sequence, Dna4Sequence("AACTTC"_dna4));
  EXPECT_EQ(contig1.sequences.second[2].annotations.size(), 0);

  const auto &contig2 = genome.contigs[1];
  EXPECT_EQ(contig2.accn, "genome2");
  EXPECT_EQ(contig2.minFragmentSize, minContigSize);
  EXPECT_EQ(contig2.totalLength, 40);
  EXPECT_EQ(contig2.numAnnotations, 0);
  EXPECT_EQ(contig2.numFragments(), 4);

  EXPECT_EQ(contig2.sequences.first.size(), 2);
  EXPECT_EQ(contig2.sequences.first[0].nullFeatureId, 2);
  EXPECT_EQ(contig2.sequences.first[0].sequence,
            Dna4Sequence("GTGTCGGAGG"_dna4));
  EXPECT_EQ(contig2.sequences.first[0].annotations.size(), 0);
  EXPECT_EQ(contig2.sequences.first[1].nullFeatureId, 2);
  EXPECT_EQ(contig2.sequences.first[1].sequence,
            Dna4Sequence("CTCCATCGAC"_dna4));
  EXPECT_EQ(contig2.sequences.first[1].annotations.size(), 0);

  EXPECT_EQ(contig2.sequences.second.size(), 2);
  EXPECT_EQ(contig2.sequences.second[0].nullFeatureId, 2);
  EXPECT_EQ(contig2.sequences.second[0].sequence,
            Dna4Sequence("CCTCCGACAC"_dna4));
  EXPECT_EQ(contig2.sequences.second[0].annotations.size(), 0);
  EXPECT_EQ(contig2.sequences.second[1].nullFeatureId, 2);
  EXPECT_EQ(contig2.sequences.second[1].sequence,
            Dna4Sequence("GTCGATGGAG"_dna4));
  EXPECT_EQ(contig2.sequences.second[1].annotations.size(), 0);
}

inline std::string MakeTempPath(const char *ext) {
  return tempio::create_temporary_file(std::filesystem::temp_directory_path(),
                                       ext);
}

class ParseChunks : public testing::Test {
protected:
  ParseChunks() : fna(MakeTempPath(".fna")) {}
  ~ParseChunks() { std::filesystem::remove(fna); }
  std::string fna;

  int writeToFna(const std::string &data) {
    std::ofstream output_file(fna);

    // Check if the file was successfully opened
    if (output_file.is_open()) {
      // Write the string data to the file using the insertion operator (<<)
      output_file << data;

      // Close the file
      output_file.close();
    } else {
      std::cerr << "Error: Unable to open the file." << std::endl;
      return 1; // Return an error code
    }

    return 0;
  }
};

/*

ACGTACGTACGTACGT ACGTACGTACGATCAG TCAGTCAGTCGTA

TACGACTGACTGACTG ATCGTACGTACGTACG TACGTACGTACGT

AGTACGTCGTACGATC GTC

GACGATCGTACGACGT ACT

ATCGTCGATGCTAGCT AGCTAGCTAGCTACGT

ACGTAGCTAGCTAGCT AGCTAGCATCGACGAT

GTACGTGCTAGCTAGC TGACTCGATGCATTAA

TTAATGCATCGAGTCA GCTAGCTAGCACGTAC

TAGTATATATAGTAGT AGT

ACTACTACTATATATA CTA

*/
TEST_F(ParseChunks, ParseChunkedGenome) {
  // dummy params
  uint8_t kmer_size = 3;
  size_t threads = 2;

  // create data
  std::string data =
      ">seq1\nACGTACGTACGTACGTACGTACGTACGATCAGTCAGTCAGTCGTA\n>"
      "seq2\nAGTACGTCGTACGATCGTC\n>seq3\nATCGTCGATGCTAGCTAGCTAGCTAGCTACGT\n>"
      "seq4\nGTACGTGCTAGCTAGCTGACTCGATGCATTAA\n>seq5\nTAGTATATATAGTAGTAGT\n";
  writeToFna(data);

  // parse data
  Colours c;
  auto xgenome = parseFilterFNA(fna, c, threads, kmer_size);

  // check IDs
  ASSERT_EQ(c.getOrAssign("seq1"), 1);
  ASSERT_EQ(c.getOrAssign("seq2"), 2);
  ASSERT_EQ(c.getOrAssign("seq3"), 3);
  ASSERT_EQ(c.getOrAssign("seq4"), 4);
  ASSERT_EQ(c.getOrAssign("seq5"), 5);

  ASSERT_EQ(xgenome.genome.contigs[0].sequences.first[0].nullFeatureId, 1);
  ASSERT_EQ(xgenome.genome.contigs[0].sequences.second[0].nullFeatureId, 1);
  ASSERT_EQ(xgenome.genome.contigs[1].sequences.first[0].nullFeatureId, 2);
  ASSERT_EQ(xgenome.genome.contigs[1].sequences.second[0].nullFeatureId, 2);
  ASSERT_EQ(xgenome.genome.contigs[2].sequences.first[0].nullFeatureId, 3);
  ASSERT_EQ(xgenome.genome.contigs[2].sequences.second[0].nullFeatureId, 3);
  ASSERT_EQ(xgenome.genome.contigs[3].sequences.first[0].nullFeatureId, 4);
  ASSERT_EQ(xgenome.genome.contigs[3].sequences.second[0].nullFeatureId, 4);
  ASSERT_EQ(xgenome.genome.contigs[4].sequences.first[0].nullFeatureId, 5);
  ASSERT_EQ(xgenome.genome.contigs[4].sequences.second[0].nullFeatureId, 5);

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

  // output sequences

  std::vector<Dna4Sequence> sequenceChunks;

  for (auto rng : xgenome.chunks) {
    for (auto seq : rng.fragments(kmer_size)) {
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
