#include <cstring>
#include <fstream>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <maki/fasta.hpp>
#include <maki/maki.h>
#include <maki/utils.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <sstream>
#include <string>
#include <utility>
#include <zstr.hpp>

using namespace seqan3::literals;

using ::testing::ContainerEq;
using ::testing::ElementsAreArray;

using GenomeFeatureIdGenerator =
    colour::encoding::FeatureIdGenerator<GenomeStatsSummary, Gff3Record>;
using GenomeFeatureIdAssigner =
    colour::encoding::FeatureIdAssigner<GenomeStatsSummary, Gff3Record>;

#define STR_PAIR(a, b) std::make_pair<std::string, std::string>(a, b)

// File parsing tests

class FileTypeTests : public testing::Test {
protected:
  FileTypeTests() {}

  std::string tarFile = "/path/to/tar.gz.file.tar.gz",
              gzFile = "/path/to/gz.file.gz",
              noCompFile = "/test/path/file.txt",
              fnaCompFile = "/test/path/file.fna.gz",
              fnaNoCompFile = "/test/path/otherFile.fna",
              fastaCompFile = "/new/path/seqFile.fasta.gz",
              fastaNoCompFile = "seqFile.fasta",
              fastaFalseTestFile = "seqFile.fasta.txt",
              gffCompFile = "/home/new/genome.gff3.gz",
              gffNoCompFile = "/home/new/genome.gff3",
              fastqCompFile = "/other/test.fq.gz";
};

TEST_F(FileTypeTests, DetectTar) {
  EXPECT_EQ(fileutils::removeCompressedExtensions(tarFile),
            "/path/to/tar.gz.file");
}

TEST_F(FileTypeTests, DetectGz) {
  EXPECT_EQ(fileutils::removeCompressedExtensions(gzFile), "/path/to/gz.file");
}

TEST_F(FileTypeTests, NoDetect) {
  EXPECT_EQ(fileutils::removeCompressedExtensions(noCompFile),
            "/test/path/file.txt");
}

TEST_F(FileTypeTests, DetectFnaGz) {
  EXPECT_EQ(fileutils::detectFileType(fnaCompFile), fileutils::FastaFileType);
}

TEST_F(FileTypeTests, DetectFna) {
  EXPECT_EQ(fileutils::detectFileType(fnaNoCompFile), fileutils::FastaFileType);
}

TEST_F(FileTypeTests, DetectFastaGz) {
  EXPECT_EQ(fileutils::detectFileType(fastaCompFile), fileutils::FastaFileType);
}

TEST_F(FileTypeTests, DetectFasta) {
  EXPECT_EQ(fileutils::detectFileType(fastaNoCompFile),
            fileutils::FastaFileType);
}

TEST_F(FileTypeTests, DetectUnknown) {
  EXPECT_EQ(fileutils::detectFileType(fastaFalseTestFile),
            fileutils::UnknownFileType);
}

TEST_F(FileTypeTests, DetectGffGz) {
  EXPECT_EQ(fileutils::detectFileType(gffCompFile), fileutils::Gff3FileType);
}

TEST_F(FileTypeTests, DetectGff) {
  EXPECT_EQ(fileutils::detectFileType(gffNoCompFile), fileutils::Gff3FileType);
}

TEST_F(FileTypeTests, DetectFastqGz) {
  EXPECT_EQ(fileutils::detectFileType(fastqCompFile), fileutils::FastQFileType);
}

TEST_F(FileTypeTests, ExtractGzName) {
  EXPECT_EQ(fileutils::extractSequenceName(fnaCompFile), "file");
}

TEST_F(FileTypeTests, ExtractName) {
  EXPECT_EQ(fileutils::extractSequenceName(fnaNoCompFile), "otherFile");
}

TEST_F(FileTypeTests, ExtractNameNoDir) {
  EXPECT_EQ(fileutils::extractSequenceName(fastaNoCompFile), "seqFile");
}

// Parse fasta

constexpr size_t minContigSize = 3;

class ParseFasta : public testing::Test {
protected:
  ParseFasta() {}
};

TEST_F(ParseFasta, ParseBitpackedFasta) {
  GenomeFeatureIdGenerator idAssigner(32); // log2(32) = 5
  auto genome =
      parseGenome(STRING(SEQA_FNA), idAssigner.specify(1), minContigSize);

  EXPECT_EQ(genome.stats.genomeName, "seqA");
  EXPECT_EQ(genome.stats.sourceFile, STRING(SEQA_FNA));
  EXPECT_EQ(genome.stats.minObsContigLength, 21);
  EXPECT_EQ(genome.stats.maxObsContigLength, 23);
  EXPECT_EQ(genome.stats.numAnnotations, 0);
  EXPECT_EQ(genome.length(), 80);
  EXPECT_EQ(genome.numContigs(), 2);
  EXPECT_EQ(genome.numFragments(), 10);

  const auto &contig1 = genome.contigs[0];
  EXPECT_EQ(contig1.contigName, "genome1");
  EXPECT_EQ(contig1.minFragmentSize, minContigSize);
  EXPECT_EQ(contig1.totalLength, 40);
  EXPECT_EQ(contig1.numAnnotations, 0);
  EXPECT_EQ(contig1.gcCount, 22);
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
  EXPECT_EQ(contig2.contigName, "genome2");
  EXPECT_EQ(contig2.minFragmentSize, minContigSize);
  EXPECT_EQ(contig2.totalLength, 40);
  EXPECT_EQ(contig2.numAnnotations, 0);
  EXPECT_EQ(contig2.gcCount, 26);
  EXPECT_EQ(contig2.numFragments(), 4);

  EXPECT_EQ(contig2.sequences.first.size(), 2);
  EXPECT_EQ(contig2.sequences.first[0].nullFeatureId, 1);
  EXPECT_EQ(contig2.sequences.first[0].sequence,
            Dna4Sequence("GTGTCGGAGG"_dna4));
  EXPECT_EQ(contig2.sequences.first[0].annotations.size(), 0);
  EXPECT_EQ(contig2.sequences.first[1].nullFeatureId, 1);
  EXPECT_EQ(contig2.sequences.first[1].sequence,
            Dna4Sequence("CTCCATCGAC"_dna4));
  EXPECT_EQ(contig2.sequences.first[1].annotations.size(), 0);

  EXPECT_EQ(contig2.sequences.second.size(), 2);
  EXPECT_EQ(contig2.sequences.second[0].nullFeatureId, 1);
  EXPECT_EQ(contig2.sequences.second[0].sequence,
            Dna4Sequence("CCTCCGACAC"_dna4));
  EXPECT_EQ(contig2.sequences.second[0].annotations.size(), 0);
  EXPECT_EQ(contig2.sequences.second[1].nullFeatureId, 1);
  EXPECT_EQ(contig2.sequences.second[1].sequence,
            Dna4Sequence("GTCGATGGAG"_dna4));
  EXPECT_EQ(contig2.sequences.second[1].annotations.size(), 0);
}

TEST_F(ParseFasta, ParseFastaParallel) {
  auto genomes = loadGenomes(
      {STRING(SEQA_FNA), STRING(SEQB_FNA), STRING(SEQC_FNA)}, minContigSize);

  EXPECT_EQ(genomes.numGenomes(), 3);
  EXPECT_EQ(genomes.getFeatureWidth(), 2);
  EXPECT_EQ(genomes.totalLength(), 202);
}

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
TEST_F(ParseFasta, ParseChunkedGenome) {
  // dummy params

  uint8_t kmer_size = 3;
  size_t threads = 2;

  // parse data

  std::string data =
      ">seq1\nACGTACGTACGTACGTACGTACGTACGATCAGTCAGTCAGTCGTA\n>"
      "seq2\nAGTACGTCGTACGATCGTC\n>seq3\nATCGTCGATGCTAGCTAGCTAGCTAGCTACGT\n>"
      "seq4\nGTACGTGCTAGCTAGCTGACTCGATGCATTAA\n>seq5\nTAGTATATATAGTAGTAGT\n";
  std::istringstream datastream(data);

  Dna4Genome genome;
  parseFastaStream(genome, datastream, 0, kmer_size);
  ASSERT_EQ(genome.medianContigSize(), 32);

  ChunkedDna4Genome xgenome(std::move(genome), threads, kmer_size);

  // check sequence metadata

  ASSERT_EQ(xgenome.chunks(), 18);

  ASSERT_FALSE(xgenome[0].endIsTerminal());
  ASSERT_EQ(xgenome[0].numTerminals(kmer_size), kmer_size);
  ASSERT_FALSE(xgenome[1].endIsTerminal());
  ASSERT_EQ(xgenome[1].numTerminals(kmer_size), 0);
  ASSERT_TRUE(xgenome[2].endIsTerminal());
  ASSERT_EQ(xgenome[2].numTerminals(kmer_size), 0);
  ASSERT_FALSE(xgenome[3].endIsTerminal());
  ASSERT_EQ(xgenome[3].numTerminals(kmer_size), kmer_size);
  ASSERT_FALSE(xgenome[4].endIsTerminal());
  ASSERT_EQ(xgenome[4].numTerminals(kmer_size), 0);
  ASSERT_TRUE(xgenome[5].endIsTerminal());
  ASSERT_EQ(xgenome[5].numTerminals(kmer_size), 0);
  ASSERT_TRUE(xgenome[6].endIsTerminal());
  ASSERT_EQ(xgenome[6].numTerminals(kmer_size), kmer_size);
  ASSERT_TRUE(xgenome[7].endIsTerminal());
  ASSERT_EQ(xgenome[7].numTerminals(kmer_size), kmer_size);
  ASSERT_FALSE(xgenome[8].endIsTerminal());
  ASSERT_EQ(xgenome[8].numTerminals(kmer_size), kmer_size);
  ASSERT_TRUE(xgenome[9].endIsTerminal());
  ASSERT_EQ(xgenome[9].numTerminals(kmer_size), 0);
  ASSERT_FALSE(xgenome[10].endIsTerminal());
  ASSERT_EQ(xgenome[10].numTerminals(kmer_size), kmer_size);
  ASSERT_TRUE(xgenome[11].endIsTerminal());
  ASSERT_EQ(xgenome[11].numTerminals(kmer_size), 0);
  ASSERT_FALSE(xgenome[12].endIsTerminal());
  ASSERT_EQ(xgenome[12].numTerminals(kmer_size), kmer_size);
  ASSERT_TRUE(xgenome[13].endIsTerminal());
  ASSERT_EQ(xgenome[13].numTerminals(kmer_size), 0);
  ASSERT_FALSE(xgenome[14].endIsTerminal());
  ASSERT_EQ(xgenome[14].numTerminals(kmer_size), kmer_size);
  ASSERT_TRUE(xgenome[15].endIsTerminal());
  ASSERT_EQ(xgenome[15].numTerminals(kmer_size), 0);
  ASSERT_TRUE(xgenome[16].endIsTerminal());
  ASSERT_EQ(xgenome[16].numTerminals(kmer_size), kmer_size);
  ASSERT_TRUE(xgenome[17].endIsTerminal());
  ASSERT_EQ(xgenome[17].numTerminals(kmer_size), kmer_size);

  // output sequences

  std::vector<Dna4Sequence> sequenceChunks;

  for (SequenceRange const &rng : xgenome) {
    sequenceChunks.emplace_back(rng.begin(), rng.end());
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

// Parse gff3

class ParseGff3Tests : public testing::Test {
protected:
  ParseGff3Tests() {}
};

TEST_F(ParseGff3Tests, ParseFromString) {
  std::string gffData =
      "genome1\tPyrodigal\tCDS\t1\t5\t.\t+\t0\tID=OOLAGM_00001;Name=TetR "
      "family transcriptional regulator;locus_tag=OOLAGM_00001;product=TetR "
      "family transcriptional "
      "regulator;Dbxref=SO:0001217,UniRef:UniRef50_A0A1Z4C077;gene=tetR";
  Gff3Record record(gffData);

  EXPECT_EQ(record.ref, "genome1");
  EXPECT_EQ(record.beginPos, 0);
  EXPECT_EQ(record.endPos, 5);
  EXPECT_EQ(record.strand, '+');

  std::vector expectedValues = {
      STR_PAIR("ID", "OOLAGM_00001"),
      STR_PAIR("Name", "TetR family transcriptional regulator"),
      STR_PAIR("locus_tag", "OOLAGM_00001"),
      STR_PAIR("product", "TetR family transcriptional regulator"),
      STR_PAIR("Dbxref", "SO:0001217,UniRef:UniRef50_A0A1Z4C077"),
      STR_PAIR("gene", "tetR")};
}

TEST_F(ParseGff3Tests, ParseFromStringNoTags) {
  std::string gffData = "genome1\tPyrodigal\tCDS\t1\t5\t1e-6\t+\t0\t";
  Gff3Record record(gffData);

  EXPECT_EQ(record.ref, "genome1");
  EXPECT_EQ(record.beginPos, 0);
  EXPECT_EQ(record.endPos, 5);
  EXPECT_EQ(record.strand, '+');
}

class TestAnnotsFileTests : public testing::Test {
protected:
  TestAnnotsFileTests() {
    zstr::ifstream gffStream(STRING(TEST_ANNOTS));
    parseAnnotationStream(annots, gffStream);
    annots.sort();
  }

  GenomeAnnotationList annots;
};

TEST_F(TestAnnotsFileTests, ParseAndSortFromFile) {
  std::vector expectedAnnotations = {
      Gff3Record(std::string{"genome1"}, 0, 5, '+'),
      Gff3Record(std::string{"genome1"}, 0, 999, '+'),
      Gff3Record(std::string{"genome1"}, 1, 49, '+'),
      Gff3Record(std::string{"genome1"}, 99, 2199, '+'),
      Gff3Record(std::string{"genome1"}, 2199, 2229, '+'),
      Gff3Record(std::string{"genome1"}, 2214, 16999, '+'),
      Gff3Record(std::string{"genome1"}, 2214, 16999, '-'),
      Gff3Record(std::string{"genome1"}, 2264, 17049, '-'),
      Gff3Record(std::string{"genome2"}, 0, 21, '+'),
      Gff3Record(std::string{"genome2"}, 0, 21, '-'),
      Gff3Record(std::string{"genome2"}, 199, 221, '-'),
      Gff3Record(std::string{"genome2"}, 209, 49999, '-')};
  ASSERT_THAT(annots.records, ElementsAreArray(expectedAnnotations));
}

TEST_F(TestAnnotsFileTests, SearchForContig) {
  auto range = annots.findContig("genome2");
  EXPECT_EQ(range.begin, annots.records.cbegin() + 8);
  EXPECT_EQ(range.end, annots.records.cend());
  EXPECT_EQ(range.size(), 4);
}

TEST_F(TestAnnotsFileTests, SearchForContigStrand) {
  auto range = annots.findContig("genome1").getStrand('-');
  EXPECT_EQ(range.begin, annots.records.cbegin() + 6);
  EXPECT_EQ(range.end, annots.records.cbegin() + 8);
  EXPECT_EQ(range.size(), 2);
}

TEST_F(TestAnnotsFileTests, SearchForContigStrandNotFound) {
  auto range = annots.findContig("genome1").getStrand('x');
  EXPECT_EQ(range.begin, annots.records.cbegin() + 8);
  EXPECT_EQ(range.end, annots.records.cbegin() + 8);
  EXPECT_EQ(range.size(), 0);
}

class ParseGffFastaTests : public testing::Test {
protected:
  ParseGffFastaTests() {}
};

TEST_F(ParseGffFastaTests, ParseGff) {
  GenomeFeatureIdGenerator idAssigner(32); // 5 < log2(32 + 1) < 6
  EXPECT_EQ(idAssigner.getBaseFeatureWidth(), 6);
  auto genome = parseAnnotatedGenome(STRING(SEQA_GFF), idAssigner.specify(1),
                                     minContigSize);

  EXPECT_EQ(genome.stats.genomeName, "seqA");
  EXPECT_EQ(genome.stats.sourceFile, STRING(SEQA_GFF));
  EXPECT_EQ(genome.stats.minObsContigLength, 21);
  EXPECT_EQ(genome.stats.maxObsContigLength, 23);
  EXPECT_EQ(genome.stats.numAnnotations, 6);
  EXPECT_EQ(genome.length(), 80);
  EXPECT_EQ(genome.numContigs(), 2);
  EXPECT_EQ(genome.numFragments(), 10);

  const auto &contig1 = genome.contigs[0];
  EXPECT_EQ(contig1.contigName, "genome1");
  EXPECT_EQ(contig1.minFragmentSize, minContigSize);
  EXPECT_EQ(contig1.totalLength, 40);
  EXPECT_EQ(contig1.numAnnotations, 3);
  EXPECT_EQ(contig1.gcCount, 22);
  EXPECT_EQ(contig1.numFragments(), 6);

  // Contig 1, Forward
  EXPECT_EQ(contig1.sequences.first.size(), 3);
  // Contig 1, Forward, Fragment 1
  EXPECT_EQ(contig1.sequences.first[0].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences.first[0].sequence, Dna4Sequence("ACGAG"_dna4));
  std::vector contig1Frag1Annots = {Gff3Token{0, 5, 0b001000001}};
  ASSERT_THAT(contig1.sequences.first[0].annotations,
              ElementsAreArray(contig1Frag1Annots));
  // Contig 1, Forward, Fragment 2
  EXPECT_EQ(contig1.sequences.first[1].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences.first[1].sequence,
            Dna4Sequence("CGGTGGCAA"_dna4));
  EXPECT_EQ(contig1.sequences.first[1].annotations.size(), 0);
  // Contig 1, Forward, Fragment 3
  EXPECT_EQ(contig1.sequences.first[2].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences.first[2].sequence, Dna4Sequence("GAAGTT"_dna4));
  std::vector contig1Frag3Annots = {Gff3Token{0, 6, 0b010000001}};
  ASSERT_THAT(contig1.sequences.first[2].annotations,
              ElementsAreArray(contig1Frag3Annots));

  // Contig 1, Reverse
  EXPECT_EQ(contig1.sequences.second.size(), 3);
  // Contig 1, Reverse, Fragment 1
  EXPECT_EQ(contig1.sequences.second[0].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences.second[0].sequence, Dna4Sequence("AACTTC"_dna4));
  EXPECT_EQ(contig1.sequences.second[0].annotations.size(), 0);
  // Contig 1, Reverse, Fragment 2
  EXPECT_EQ(contig1.sequences.second[1].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences.second[1].sequence,
            Dna4Sequence("TTGCCACCG"_dna4));
  std::vector contig1Frag2Annots = {Gff3Token{0, 5, 0b011000001}};
  ASSERT_THAT(contig1.sequences.second[1].annotations,
              ElementsAreArray(contig1Frag2Annots));

  // Contig 1, Reverse, Fragment 3
  EXPECT_EQ(contig1.sequences.second[2].nullFeatureId, 1);
  EXPECT_EQ(contig1.sequences.second[2].sequence, Dna4Sequence("CTCGT"_dna4));
  EXPECT_EQ(contig1.sequences.second[2].annotations.size(), 0);

  const auto &contig2 = genome.contigs[1];
  EXPECT_EQ(contig2.contigName, "genome2");
  EXPECT_EQ(contig2.minFragmentSize, minContigSize);
  EXPECT_EQ(contig2.totalLength, 40);
  EXPECT_EQ(contig2.numAnnotations, 3);
  EXPECT_EQ(contig2.gcCount, 26);
  EXPECT_EQ(contig2.numFragments(), 4);

  // Contig 2, Forward
  EXPECT_EQ(contig2.sequences.first.size(), 2);
  // Contig 2, Forward, Fragment 1
  EXPECT_EQ(contig2.sequences.first[0].nullFeatureId, 1);
  EXPECT_EQ(contig2.sequences.first[0].sequence,
            Dna4Sequence("GTGTCGGAGG"_dna4));
  std::vector contig2Frag1Annots = {Gff3Token{0, 6, 0b100000001},
                                    Gff3Token{1, 10, 0b101000001}};
  ASSERT_THAT(contig2.sequences.first[0].annotations,
              ElementsAreArray(contig2Frag1Annots));
  // Contig 2, Forward, Fragment 2
  EXPECT_EQ(contig2.sequences.first[1].nullFeatureId, 1);
  EXPECT_EQ(contig2.sequences.first[1].sequence,
            Dna4Sequence("CTCCATCGAC"_dna4));
  EXPECT_EQ(contig2.sequences.first[1].annotations.size(), 0);

  // Contig 2, Reverse
  EXPECT_EQ(contig2.sequences.second.size(), 2);
  // Contig 2, Reverse, Fragment 1
  EXPECT_EQ(contig2.sequences.second[0].nullFeatureId, 1);
  EXPECT_EQ(contig2.sequences.second[0].sequence,
            Dna4Sequence("GTCGATGGAG"_dna4));
  EXPECT_EQ(contig2.sequences.second[0].annotations.size(), 0);
  // Contig 2, Reverse, Fragment 2
  EXPECT_EQ(contig2.sequences.second[1].nullFeatureId, 1);
  EXPECT_EQ(contig2.sequences.second[1].sequence,
            Dna4Sequence("CCTCCGACAC"_dna4));
  std::vector contig2Frag2RevAnnots = {Gff3Token{0, 10, 0b110000001}};
  ASSERT_THAT(contig2.sequences.second[1].annotations,
              ElementsAreArray(contig2Frag2RevAnnots));
}

TEST_F(ParseGffFastaTests, ParseGffParallel) {
  auto genomes = loadGenomes(
      {STRING(SEQA_GFF), STRING(SEQB_GFF), STRING(SEQC_GFF)}, minContigSize);

  EXPECT_EQ(genomes.numGenomes(), 3);
  EXPECT_EQ(genomes.getFeatureWidth(), 6);
  EXPECT_EQ(genomes.totalLength(), 202);

  EXPECT_EQ(genomes[0].stats.genomeName, "seqA");
  EXPECT_EQ(genomes[0].stats.sourceFile, STRING(SEQA_GFF));
  EXPECT_EQ(genomes[0].stats.minObsContigLength, 21);
  EXPECT_EQ(genomes[0].stats.maxObsContigLength, 23);
  EXPECT_EQ(genomes[0].stats.numAnnotations, 6);
  EXPECT_EQ(genomes[0].length(), 80);
  EXPECT_EQ(genomes[0].numContigs(), 2);
  EXPECT_EQ(genomes[0].numFragments(), 10);

  EXPECT_EQ(genomes[1].stats.genomeName, "seqB");
  EXPECT_EQ(genomes[1].stats.sourceFile, STRING(SEQB_GFF));
  EXPECT_EQ(genomes[1].stats.minObsContigLength, 11);
  EXPECT_EQ(genomes[1].stats.maxObsContigLength, 12);
  EXPECT_EQ(genomes[1].stats.numAnnotations, 8);
  EXPECT_EQ(genomes[1].length(), 86);
  EXPECT_EQ(genomes[1].numContigs(), 4);
  EXPECT_EQ(genomes[1].numFragments(), 10);

  EXPECT_EQ(genomes[2].stats.genomeName, "seqC");
  EXPECT_EQ(genomes[2].stats.sourceFile, STRING(SEQC_GFF));
  EXPECT_EQ(genomes[2].stats.minObsContigLength, 18);
  EXPECT_EQ(genomes[2].stats.maxObsContigLength, 18);
  EXPECT_EQ(genomes[2].stats.numAnnotations, 2);
  EXPECT_EQ(genomes[2].length(), 36);
  EXPECT_EQ(genomes[2].numContigs(), 1);
  EXPECT_EQ(genomes[2].numFragments(), 2);
}
