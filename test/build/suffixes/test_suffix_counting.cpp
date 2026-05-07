#include "maki/build/io/fasta.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/core/utils/tempfile.hpp"
#include <cstddef>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <vector>

using ::testing::ElementsAre;
using namespace seqan3::literals;

constexpr size_t k = 3;

std::string MakeTempPath(const char *ext) {
  return tempio::create_temporary_file(std::filesystem::temp_directory_path(),
                                       ext);
}

int writeToFna(const std::string &data, const std::string &file) {
  std::ofstream output_file(file);

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

class CountGffSuffixTests : public testing::Test {
protected:
  CountGffSuffixTests() : c(), genome(parseGFF(STRING(SEQA_GFF), c, k)) {}

  Colours c;
  Dna4Genome genome;
};

TEST_F(CountGffSuffixTests, Contig1F1S1Forward) {
  uint8_t s = 1;
  SuffixTable table;
  table.resize(s);

  auto &fragment = genome.contigs[0].sequences[0][0];
  table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

  ASSERT_THAT(table.cdata(), ElementsAre(1, 0, 2, 0));
}

TEST_F(CountGffSuffixTests, Contig1F2S1Forward) {
  uint8_t s = 1;
  SuffixTable table;
  table.resize(s);

  auto &fragment = genome.contigs[0].sequences[0][1];
  table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

  ASSERT_THAT(table.cdata(), ElementsAre(2, 1, 3, 1));
}

TEST_F(CountGffSuffixTests, Contig1F3S1Forward) {
  uint8_t s = 1;
  SuffixTable table;
  table.resize(s);

  auto &fragment = genome.contigs[0].sequences[0][2];
  table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

  ASSERT_THAT(table.cdata(), ElementsAre(1, 0, 1, 2));
}

TEST_F(CountGffSuffixTests, Contig1F1S1Reverse) {
  uint8_t s = 1;
  SuffixTable table;
  table.resize(s);

  auto &fragment = genome.contigs[0].sequences[1][0];
  table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

  ASSERT_THAT(table.cdata(), ElementsAre(0, 2, 0, 2));
}

TEST_F(CountGffSuffixTests, Contig1F2S1Reverse) {
  uint8_t s = 1;
  SuffixTable table;
  table.resize(s);

  auto &fragment = genome.contigs[0].sequences[1][1];
  table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

  ASSERT_THAT(table.cdata(), ElementsAre(1, 4, 2, 0));
}

TEST_F(CountGffSuffixTests, Contig1F3S1Reverse) {
  uint8_t s = 1;
  SuffixTable table;
  table.resize(s);

  auto &fragment = genome.contigs[0].sequences[1][2];
  table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

  ASSERT_THAT(table.cdata(), ElementsAre(0, 1, 1, 1));
}

TEST_F(CountGffSuffixTests, Contig2F1S1Forward) {
  uint8_t s = 1;
  SuffixTable table;
  table.resize(s);

  auto &fragment = genome.contigs[1].sequences[0][0];
  table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

  ASSERT_THAT(table.cdata(), ElementsAre(1, 1, 5, 1));
}

TEST_F(CountGffSuffixTests, Contig2F1S2Forward) {
  uint8_t s = 2;
  SuffixTable table;
  table.resize(s);

  auto &fragment = genome.contigs[1].sequences[0][0];
  table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

  ASSERT_THAT(table.cdata(),
              ElementsAre(0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 2, 1, 0, 0, 1, 0));
}

TEST_F(CountGffSuffixTests, Contig2F2S1Forward) {
  uint8_t s = 1;
  SuffixTable table;
  table.resize(s);

  auto &fragment = genome.contigs[1].sequences[0][1];
  table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

  ASSERT_THAT(table.cdata(), ElementsAre(2, 4, 1, 1));
}

TEST_F(CountGffSuffixTests, Contig2F2S2Forward) {
  uint8_t s = 2;
  SuffixTable table;
  table.resize(s);

  auto &fragment = genome.contigs[1].sequences[0][1];
  table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

  ASSERT_THAT(table.cdata(),
              ElementsAre(0, 1, 1, 0, 1, 1, 0, 2, 0, 1, 0, 0, 1, 0, 0, 0));
}

TEST_F(CountGffSuffixTests, Contig2F1S1Reverse) {
  uint8_t s = 1;
  SuffixTable table;
  table.resize(s);

  auto &fragment = genome.contigs[1].sequences[1][0];
  table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

  ASSERT_THAT(table.cdata(), ElementsAre(2, 1, 4, 1));
}

TEST_F(CountGffSuffixTests, Contig2F1S2Reverse) {
  uint8_t s = 2;
  SuffixTable table;
  table.resize(s);

  auto &fragment = genome.contigs[1].sequences[1][0];
  table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

  ASSERT_THAT(table.cdata(),
              ElementsAre(0, 0, 2, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0));
}

TEST_F(CountGffSuffixTests, Contig2F2S1Reverse) {
  uint8_t s = 1;
  SuffixTable table;
  table.resize(s);

  auto &fragment = genome.contigs[1].sequences[1][1];
  table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

  ASSERT_THAT(table.cdata(), ElementsAre(2, 4, 1, 1));
}

TEST_F(CountGffSuffixTests, Contig2F2S2Reverse) {
  uint8_t s = 2;
  SuffixTable table;
  table.resize(s);

  auto &fragment = genome.contigs[1].sequences[1][1];
  table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

  ASSERT_THAT(table.cdata(),
              ElementsAre(0, 1, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0));
}

class ChunkedGenomeSuffixTests : public testing::Test {
protected:
  ChunkedGenomeSuffixTests()
      : c(), xgenome(parseFilterFNA(STRING(SEQA_FNA), c, 2, k)) {}

  Colours c;
  ChunkedDna4Genome xgenome;
};

TEST_F(ChunkedGenomeSuffixTests, CountChunkedGenome) {
  constexpr uint8_t s = 1;

  std::vector<const SequenceContainer *> view;
  view.reserve(xgenome.chunks.size());
  for (const auto &chunk : xgenome.chunks) {
    view.push_back(&chunk);
  }

  auto counts = createSuffixPlan(view, k, s, k - s);
  ASSERT_THAT(counts.back().cdata(), ElementsAre(12, 18, 20, 10));
}

class ChunkedGenomeChunkTests : public testing::Test {
protected:
  ChunkedGenomeChunkTests()
      : fna(MakeTempPath(".fna")), c(), xgenome(makeGenome()) {}
  ~ChunkedGenomeChunkTests() { std::filesystem::remove(fna); }

  ChunkedDna4Genome makeGenome() {
    // create data
    std::string data =
        ">seq1\nACGTACGTACGTACGTACGTACGTACGATCAGTCAGTCAGTCGTA\n>"
        "seq2\nAGTACGTCGTACGATCGTC\n>seq3\nATCGTCGATGCTAGCTAGCTAGCTAGCTACGT\n>"
        "seq4\nGTACGTGCTAGCTAGCTGACTCGATGCATTAA\n>seq5\nTAGTATATATAGTAGTAGT\n";
    writeToFna(data, fna);

    // parse data
    return parseFilterFNA(fna, c, 2, k);
  }

  std::vector<const SequenceContainer *> getView() const {
    std::vector<const SequenceContainer *> view;
    view.reserve(xgenome.chunks.size());
    for (const auto &chunk : xgenome.chunks) {
      view.push_back(&chunk);
    }
    return view;
  }

  std::string fna;
  Colours c;
  ChunkedDna4Genome xgenome;
};

TEST_F(ChunkedGenomeChunkTests, CountChunks) {
  constexpr uint8_t s = 2;
  auto view = getView();
  auto counts = createSuffixPlan(view, k, s, k - s, false);
  ASSERT_EQ(counts.size(), 18);

  ASSERT_THAT(counts[0].cdata(),
              ElementsAre(0, 0, 0, 4, 4, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0));
  ASSERT_THAT(counts[1].cdata(),
              ElementsAre(0, 1, 1, 2, 2, 0, 0, 2, 1, 3, 0, 0, 1, 0, 3, 0));
  ASSERT_THAT(counts[2].cdata(),
              ElementsAre(0, 2, 0, 1, 0, 0, 0, 2, 2, 1, 0, 0, 0, 0, 3, 0));
  ASSERT_THAT(counts[3].cdata(),
              ElementsAre(0, 0, 4, 0, 4, 0, 0, 0, 0, 1, 0, 3, 1, 3, 0, 0));
  ASSERT_THAT(counts[4].cdata(),
              ElementsAre(0, 0, 0, 4, 3, 0, 0, 1, 0, 4, 0, 0, 0, 0, 4, 0));
  ASSERT_THAT(counts[5].cdata(),
              ElementsAre(0, 0, 0, 2, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 0));
  ASSERT_THAT(counts[6].cdata(),
              ElementsAre(0, 0, 1, 2, 2, 0, 0, 3, 0, 4, 0, 0, 1, 0, 4, 0));
  ASSERT_THAT(counts[7].cdata(),
              ElementsAre(0, 0, 2, 2, 4, 0, 0, 1, 0, 4, 0, 0, 1, 1, 2, 0));
  ASSERT_THAT(counts[8].cdata(),
              ElementsAre(0, 0, 1, 2, 0, 0, 2, 2, 2, 2, 0, 1, 1, 2, 1, 0));
  ASSERT_THAT(counts[9].cdata(),
              ElementsAre(0, 0, 0, 3, 1, 0, 3, 0, 2, 1, 0, 0, 0, 3, 1, 0));
  ASSERT_THAT(counts[10].cdata(),
              ElementsAre(0, 0, 0, 4, 0, 0, 3, 0, 4, 1, 0, 0, 0, 3, 1, 0));
  ASSERT_THAT(counts[11].cdata(),
              ElementsAre(0, 1, 2, 1, 1, 0, 2, 1, 1, 2, 0, 0, 2, 1, 0, 0));
  ASSERT_THAT(counts[12].cdata(),
              ElementsAre(0, 0, 0, 3, 1, 0, 3, 0, 2, 1, 0, 2, 0, 3, 1, 0));
  ASSERT_THAT(counts[13].cdata(),
              ElementsAre(1, 1, 2, 1, 1, 0, 1, 1, 0, 1, 0, 1, 2, 1, 0, 1));
  ASSERT_THAT(counts[14].cdata(),
              ElementsAre(1, 2, 1, 1, 0, 0, 2, 2, 2, 1, 0, 1, 2, 0, 1, 0));
  ASSERT_THAT(counts[15].cdata(),
              ElementsAre(0, 1, 0, 3, 2, 0, 2, 0, 2, 1, 0, 0, 0, 2, 1, 0));
  ASSERT_THAT(counts[16].cdata(),
              ElementsAre(0, 0, 0, 6, 0, 0, 0, 0, 4, 0, 0, 0, 3, 0, 4, 0));
  ASSERT_THAT(counts[17].cdata(),
              ElementsAre(0, 0, 0, 7, 3, 0, 0, 0, 0, 0, 0, 0, 3, 4, 0, 0));
}

TEST_F(ChunkedGenomeChunkTests, CountTotal) {
  constexpr uint8_t s = 2;
  auto view = getView();
  auto counts = createSuffixPlan(view, k, s, k - s);
  ASSERT_THAT(counts.back().cdata(), ElementsAre(2, 8, 14, 48, 31, 0, 18, 15,
                                                 22, 34, 0, 8, 17, 23, 33, 1));
}
