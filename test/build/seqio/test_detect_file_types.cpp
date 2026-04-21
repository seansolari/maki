
#include "maki/core/seq/io.hpp"
#include <gtest/gtest.h>

TEST(FileTypeTests, DetectTar) {
  EXPECT_EQ(removeCompressedExtensions("/path/to/tar.gz.file.tar.gz"),
            "/path/to/tar.gz.file");
}

TEST(FileTypeTests, DetectGz) {
  EXPECT_EQ(removeCompressedExtensions("/path/to/gz.file.gz"),
            "/path/to/gz.file");
}

TEST(FileTypeTests, NoDetect) {
  EXPECT_EQ(removeCompressedExtensions("/test/path/file.txt"),
            "/test/path/file.txt");
}

TEST(FileTypeTests, DetectFnaGz) {
  EXPECT_EQ(detectFileType("/test/path/file.fna.gz"), FastaFileType);
}

TEST(FileTypeTests, DetectFna) {
  EXPECT_EQ(detectFileType("/test/path/otherFile.fna"), FastaFileType);
}

TEST(FileTypeTests, DetectFastaGz) {
  EXPECT_EQ(detectFileType("/new/path/seqFile.fasta.gz"), FastaFileType);
}

TEST(FileTypeTests, DetectFasta) {
  EXPECT_EQ(detectFileType("seqFile.fasta"), FastaFileType);
}

TEST(FileTypeTests, DetectUnknown) {
  EXPECT_EQ(detectFileType("seqFile.fasta.txt"), UnknownFileType);
}

TEST(FileTypeTests, DetectGffGz) {
  EXPECT_EQ(detectFileType("/home/new/genome.gff3.gz"), Gff3FileType);
}

TEST(FileTypeTests, DetectGff) {
  EXPECT_EQ(detectFileType("/home/new/genome.gff3"), Gff3FileType);
}

TEST(FileTypeTests, DetectFastqGz) {
  EXPECT_EQ(detectFileType("/other/test.fq.gz"), FastQFileType);
}

TEST(FileTypeTests, ExtractGzName) {
  EXPECT_EQ(extractSequenceName("/test/path/file.fna.gz"), "file");
}

TEST(FileTypeTests, ExtractName) {
  EXPECT_EQ(extractSequenceName("/test/path/otherFile.fna"), "otherFile");
}

TEST(FileTypeTests, ExtractNameNoDir) {
  EXPECT_EQ(extractSequenceName("/new/path/seqFile.fasta.gz"), "seqFile");
}
