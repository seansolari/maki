
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <string_view>

#include "maki/build/io/fasta.hpp"
#include "maki/maki.h"

using ::testing::ElementsAreArray;

// Helper to produce a valid 9-field GFF line.
// GFF fields: 0:seqid 1:source 2:type 3:start 4:end 5:score 6:strand 7:phase
// 8:attributes
static std::string MakeGffLine(std::string_view seqid, std::string_view start,
                               std::string_view end, std::string_view strand,
                               std::string_view attrs,
                               std::string_view source = "src",
                               std::string_view type = "feature",
                               std::string_view score = ".",
                               std::string_view phase = ".") {
  std::string line;
  line.reserve(seqid.size() + source.size() + type.size() + start.size() +
               end.size() + score.size() + strand.size() + phase.size() +
               attrs.size() + 16);
  line.append(seqid).push_back('\t');
  line.append(source).push_back('\t');
  line.append(type).push_back('\t');
  line.append(start).push_back('\t');
  line.append(end).push_back('\t');
  line.append(score).push_back('\t');
  line.append(strand).push_back('\t');
  line.append(phase).push_back('\t');
  line.append(attrs);
  return line;
}

// --- Convenience assertions ---
static void
ExpectSuccessWithId(const std::pair<gffRecord, std::optional<std::string>> &res,
                    std::string_view exp_seqid, int64_t exp_begin,
                    int64_t exp_end, char exp_strand, std::string_view exp_id) {
  ASSERT_TRUE(res.second.has_value())
      << "Expected optional to be engaged (ID present).";
  const auto &[rec, id] = res;
  EXPECT_EQ(rec.accn, exp_seqid);
  EXPECT_EQ(rec.begin + 1, exp_begin);
  EXPECT_EQ(rec.end, exp_end);
  EXPECT_EQ(rec.strand, exp_strand);
  EXPECT_EQ(*id, exp_id);
}

static void
ExpectNoValue(const std::pair<gffRecord, std::optional<std::string>> &res) {
  EXPECT_FALSE(res.second.has_value())
      << "Expected optional to be disengaged (no usable ID / ignorable).";
}

// ======================= Success Cases =======================

TEST(ParseGff_Success, BasicValidLineWithID) {
  const std::string line = MakeGffLine("chr1", "11869", "14412", "+",
                                       "ID=gene:ENSG00000223972;Name=DDX11L1");
  auto res = parseGffLine(line);
  ExpectSuccessWithId(res, "chr1", 11869, 14412, '+',
                      "chr1-gene:ENSG00000223972");
}

TEST(ParseGff_Success, TrimmingWhitespaceAndCRLF) {
  std::string line =
      std::string(" \t") +
      MakeGffLine("chr2", "1", "2", "-", "  ID = G2  ; Other=foo  ") + "\r\n";
  auto res = parseGffLine(line);
  ExpectSuccessWithId(res, "chr2", 1, 2, '-', "chr2-G2");
}

TEST(ParseGff_Success, StrandVariantsAllowed) {
  {
    auto res = parseGffLine(MakeGffLine("chr3", "10", "20", ".", "ID=X"));
    ExpectSuccessWithId(res, "chr3", 10, 20, '.', "chr3-X");
  }
  {
    auto res = parseGffLine(MakeGffLine("chr3", "10", "20", "?", "ID=Y"));
    ExpectSuccessWithId(res, "chr3", 10, 20, '?', "chr3-Y");
  }
  {
    auto res = parseGffLine(MakeGffLine("chr3", "10", "20", "+", "ID=Z"));
    ExpectSuccessWithId(res, "chr3", 10, 20, '+', "chr3-Z");
  }
  {
    auto res = parseGffLine(MakeGffLine("chr3", "10", "20", "-", "ID=W"));
    ExpectSuccessWithId(res, "chr3", 10, 20, '-', "chr3-W");
  }
}

TEST(ParseGff_Success, StartEqualsEndIsAllowed) {
  auto res = parseGffLine(MakeGffLine("chrX", "42", "42", "+", "ID=SINGLE"));
  ExpectSuccessWithId(res, "chrX", 42, 42, '+', "chrX-SINGLE");
}

TEST(ParseGff_Success, IDAtEndOfAttributes) {
  auto res = parseGffLine(
      MakeGffLine("chr7", "100", "200", "+",
                  "Name=foo;Note=something;Dbxref=db:123;ID=ABC123"));
  ExpectSuccessWithId(res, "chr7", 100, 200, '+', "chr7-ABC123");
}

TEST(ParseGff_Success, MultipleAttributesOrderIrrelevant) {
  auto res = parseGffLine(MakeGffLine("scaffold_1", "5", "9", "-",
                                      "Note=foo;ID=X-1;Parent=P;Name=name"));
  ExpectSuccessWithId(res, "scaffold_1", 5, 9, '-', "scaffold_1-X-1");
}

TEST(ParseGff_Success, MultipleIDFirstWins) {
  auto res = parseGffLine(
      MakeGffLine("chr9", "11", "12", "+", "ID=FIRST;ID=SECOND;Name=bar"));
  ExpectSuccessWithId(res, "chr9", 11, 12, '+', "chr9-FIRST");
}

TEST(ParseGff_Success, PreservesPercentEncoding) {
  auto res = parseGffLine(
      MakeGffLine("chrM", "900", "950", "+", "ID=gene%3AABC%2F123;Name=raw"));
  // Expect NO decoding; ID string should contain percent escapes as-is.
  ExpectSuccessWithId(res, "chrM", 900, 950, '+', "chrM-gene%3AABC%2F123");
}

TEST(ParseGff_Success, QuotedIDIsPreservedAsIs) {
  // GFF3 doesn't require quotes, but if present in your data, we keep them.
  auto res = parseGffLine(MakeGffLine("tig0001", "1000", "2000", "-",
                                      "ID=\"weird quoted id\";Name=x"));
  ExpectSuccessWithId(res, "tig0001", 1000, 2000, '-',
                      "tig0001-\"weird quoted id\"");
}

TEST(ParseGff_Success, LeadingZerosInCoordinatesAreAccepted) {
  auto res = parseGffLine(MakeGffLine("chr10", "0005", "0010", "+", "ID=Z05"));
  ExpectSuccessWithId(res, "chr10", 5, 10, '+', "chr10-Z05");
}

// ======================= Optional (No ID / Ignorable) =======================

TEST(ParseGff_Optional, MissingID_YieldsNullopt) {
  auto res =
      parseGffLine(MakeGffLine("chr1", "1", "10", "+", "Name=noid;Parent=p"));
  ExpectNoValue(res);
}

TEST(ParseGff_Optional, ParentOnly_YieldsNullopt) {
  auto res = parseGffLine(MakeGffLine("chr1", "1", "10", "+", "Parent=p1,p2"));
  ExpectNoValue(res);
}

TEST(ParseGff_Optional, GtfLikeGeneIdButNoID_YieldsNullopt) {
  auto res = parseGffLine(
      MakeGffLine("chr1", "1", "10", "+", "gene_id=G1;transcript_id=T1"));
  ExpectNoValue(res);
}

TEST(ParseGff_Optional, EmptyIDValue_TreatedAsMissing) {
  // If your implementation treats empty value as present, change to
  // EXPECT_TRUE(res.has_value())
  auto res = parseGffLine(MakeGffLine("chr1", "1", "10", "+", "ID=;Name=foo"));
  ExpectNoValue(res);
}

TEST(ParseGff_Optional, CommentLine_YieldsNullopt) {
  auto res = parseGffLine("# A comment line that should be ignored");
  ExpectNoValue(res);
}

TEST(ParseGff_Optional, BlankLine_YieldsNullopt) {
  auto res = parseGffLine(" \t\r\n");
  ExpectNoValue(res);
}

// ======================= Structural Errors (Throw) =======================

TEST(ParseGff_Errors, NotNineFields_Throws) {
  std::string line = "chr1\tsrc\ttype\t1\t100\t.\t+\t.\n"; // only 8 fields
  EXPECT_THROW({ auto _ = parseGffLine(line); }, std::invalid_argument);
}

TEST(ParseGff_Errors, NonIntegerStart_Throws) {
  auto line = MakeGffLine("chr1", "abc", "100", "+", "ID=X");
  EXPECT_THROW({ auto _ = parseGffLine(line); }, std::invalid_argument);
}

TEST(ParseGff_Errors, NonIntegerEnd_Throws) {
  auto line = MakeGffLine("chr1", "1", "100x", "+", "ID=X");
  EXPECT_THROW({ auto _ = parseGffLine(line); }, std::invalid_argument);
}

TEST(ParseGff_Errors, NegativeStart_Throws) {
  auto line = MakeGffLine("chr1", "-5", "10", "+", "ID=X");
  EXPECT_THROW({ auto _ = parseGffLine(line); }, std::invalid_argument);
}

TEST(ParseGff_Errors, ZeroStart_Throws) {
  auto line = MakeGffLine("chr1", "0", "10", "+", "ID=X");
  EXPECT_THROW({ auto _ = parseGffLine(line); }, std::invalid_argument);
}

TEST(ParseGff_Errors, StartGreaterThanEnd_Throws) {
  auto line = MakeGffLine("chr1", "200", "100", "+", "ID=X");
  EXPECT_THROW({ auto _ = parseGffLine(line); }, std::invalid_argument);
  // If you swap instead of throwing, replace with:
  // auto res = parseGffLine(line);
  // ExpectSuccessWithId(res, "chr1", 100, 200, '+', "chr1-X");
}

TEST(ParseGff_Errors, InvalidStrandCharacter_Throws) {
  auto line = MakeGffLine("chr1", "1", "10", "*", "ID=X");
  EXPECT_THROW({ auto _ = parseGffLine(line); }, std::invalid_argument);
  // If you map invalid to '.', assert success with '.' instead.
}

TEST(ParseGff_Errors, EmptySeqId_Throws) {
  auto line = MakeGffLine("", "1", "10", "+", "ID=X");
  EXPECT_THROW({ auto _ = parseGffLine(line); }, std::invalid_argument);
}

TEST(ParseGff_Errors, Int64Overflow_Throws) {
  // 19+ digits to exceed int64_t max
  auto line = MakeGffLine("chr1", "9223372036854775808", "9223372036854775809",
                          "+", "ID=BIG");
  EXPECT_THROW({ auto _ = parseGffLine(line); }, std::invalid_argument);
}

// ======================= Robustness / Corner Behaviors =======================

TEST(ParseGff_Robustness, AttributesWithExtraSemicolonsAndSpaces) {
  auto res = parseGffLine(
      MakeGffLine("chr2", "5", "6", "+", ";;  ID = A1  ; ; Name = n ;;"));
  ExpectSuccessWithId(res, "chr2", 5, 6, '+', "chr2-A1");
}

TEST(ParseGff_Robustness, IDWithHyphensAndColons) {
  auto res = parseGffLine(
      MakeGffLine("scaf-01", "10", "20", "-", "ID=gene:ABC-123;foo=bar"));
  ExpectSuccessWithId(res, "scaf-01", 10, 20, '-', "scaf-01-gene:ABC-123");
}

TEST(ParseGff_Robustness, LargeAccessionsAndAttributes) {
  std::string acc(1000, 'A');
  std::string idv(1000, 'Z');
  auto res = parseGffLine(MakeGffLine(acc, "1", "2", "+", "ID=" + idv));
  // ID is "<accession>-<ID>"
  std::string expected = acc + "-" + idv;
  ASSERT_TRUE(res.second.has_value());
  EXPECT_EQ(res.first.accn, acc);
  EXPECT_EQ(*res.second, expected);
}

TEST(ParseGff_Robustness, TrailingWhitespaceInAttributes) {
  auto line = MakeGffLine("chr3", "7", "9", "+", "ID=JJJ   ") + "   ";
  auto res = parseGffLine(line);
  ExpectSuccessWithId(res, "chr3", 7, 9, '+', "chr3-JJJ");
}

TEST(ParseGff_Robustness, CROnlyOrLFOnlyAreTrimmed) {
  auto line = MakeGffLine("chr4", "8", "12", "+", "ID=A") + "\r";
  auto res = parseGffLine(line);
  ExpectSuccessWithId(res, "chr4", 8, 12, '+', "chr4-A");
}

// ======================= Regression-like checks =======================

TEST(ParseGff_Regression, IDEarlyWithOtherAttrsContainingEquals) {
  // Ensure '=' in other attrs doesn't break ID extraction
  auto res = parseGffLine(
      MakeGffLine("chr5", "50", "60", "+", "ID=K;Note=a=b=c;Other=x=y"));
  ExpectSuccessWithId(res, "chr5", 50, 60, '+', "chr5-K");
}

TEST(ParseGff_Regression, MultipleSeparatorsAndEmptyAttributes) {
  auto res = parseGffLine(
      MakeGffLine("chr6", "3", "4", "+", ";; ;ID=I1; ;; Name=n ;; ;"));
  ExpectSuccessWithId(res, "chr6", 3, 4, '+', "chr6-I1");
}

class TestAnnotsFileTests : public testing::Test {
protected:
  TestAnnotsFileTests() {
    zstr::ifstream gffStream(STRING(TEST_ANNOTS));
    Colours c;
    parseAnnotationStream(annots, c, gffStream);
    annots.sort();
  }

  GenomeAnnotationList annots;
};

TEST_F(TestAnnotsFileTests, ParseAndSortFromFile) {
  std::vector expectedAnnotations = {
      gffRecord(std::string{"genome1"}, 0, 5, '+', 1),
      gffRecord(std::string{"genome1"}, 0, 999, '+', 1),
      gffRecord(std::string{"genome1"}, 1, 49, '+', 1),
      gffRecord(std::string{"genome1"}, 99, 2199, '+', 1),
      gffRecord(std::string{"genome1"}, 2199, 2229, '+', 1),
      gffRecord(std::string{"genome1"}, 2214, 16999, '+', 1),
      gffRecord(std::string{"genome1"}, 2214, 16999, '-', 1),
      gffRecord(std::string{"genome1"}, 2264, 17049, '-', 1),
      gffRecord(std::string{"genome2"}, 0, 21, '+', 2),
      gffRecord(std::string{"genome2"}, 0, 21, '-', 2),
      gffRecord(std::string{"genome2"}, 199, 221, '-', 2),
      gffRecord(std::string{"genome2"}, 209, 49999, '-', 2)};
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

  const size_t minContigSize = 3;
};

TEST_F(ParseGffFastaTests, ParseGff) {
  Colours c;
  auto genome = parseGFF(STRING(SEQA_GFF), c, minContigSize);

  EXPECT_EQ(genome.length(), 80);
  EXPECT_EQ(genome.numFragments(), 10);

  const auto &contig1 = genome.contigs[0];
  EXPECT_EQ(contig1.accn, "genome1");
  EXPECT_EQ(contig1.minFragmentSize, minContigSize);
  EXPECT_EQ(contig1.totalLength, 40);
  EXPECT_EQ(contig1.numAnnotations, 3);
  EXPECT_EQ(contig1.numFragments(), 6);

  // Contig 1, Forward
  EXPECT_EQ(contig1.sequences.first.size(), 3);
  // Contig 1, Forward, Fragment 1
  EXPECT_EQ(contig1.sequences.first[0].nullFeatureId, 0);
  EXPECT_EQ(contig1.sequences.first[0].sequence, Dna4Sequence("ACGAG"_dna4));
  std::vector contig1Frag1Annots = {gffToken{1, 0, 5}};
  ASSERT_THAT(contig1.sequences.first[0].annotations,
              ElementsAreArray(contig1Frag1Annots));
  // Contig 1, Forward, Fragment 2
  EXPECT_EQ(contig1.sequences.first[1].nullFeatureId, 0);
  EXPECT_EQ(contig1.sequences.first[1].sequence,
            Dna4Sequence("CGGTGGCAA"_dna4));
  EXPECT_EQ(contig1.sequences.first[1].annotations.size(), 0);
  // Contig 1, Forward, Fragment 3
  EXPECT_EQ(contig1.sequences.first[2].nullFeatureId, 0);
  EXPECT_EQ(contig1.sequences.first[2].sequence, Dna4Sequence("GAAGTT"_dna4));
  std::vector contig1Frag3Annots = {gffToken{2, 0, 6}};
  ASSERT_THAT(contig1.sequences.first[2].annotations,
              ElementsAreArray(contig1Frag3Annots));

  // Contig 1, Reverse
  EXPECT_EQ(contig1.sequences.second.size(), 3);
  // Contig 1, Reverse, Fragment 1
  EXPECT_EQ(contig1.sequences.second[0].nullFeatureId, 0);
  EXPECT_EQ(contig1.sequences.second[0].sequence, Dna4Sequence("AACTTC"_dna4));
  EXPECT_EQ(contig1.sequences.second[0].annotations.size(), 0);
  // Contig 1, Reverse, Fragment 2
  EXPECT_EQ(contig1.sequences.second[1].nullFeatureId, 0);
  EXPECT_EQ(contig1.sequences.second[1].sequence,
            Dna4Sequence("TTGCCACCG"_dna4));
  std::vector contig1Frag2Annots = {gffToken{3, 0, 5}};
  ASSERT_THAT(contig1.sequences.second[1].annotations,
              ElementsAreArray(contig1Frag2Annots));

  // Contig 1, Reverse, Fragment 3
  EXPECT_EQ(contig1.sequences.second[2].nullFeatureId, 0);
  EXPECT_EQ(contig1.sequences.second[2].sequence, Dna4Sequence("CTCGT"_dna4));
  EXPECT_EQ(contig1.sequences.second[2].annotations.size(), 0);

  const auto &contig2 = genome.contigs[1];
  EXPECT_EQ(contig2.accn, "genome2");
  EXPECT_EQ(contig2.minFragmentSize, minContigSize);
  EXPECT_EQ(contig2.totalLength, 40);
  EXPECT_EQ(contig2.numAnnotations, 3);
  EXPECT_EQ(contig2.numFragments(), 4);

  // Contig 2, Forward
  EXPECT_EQ(contig2.sequences.first.size(), 2);
  // Contig 2, Forward, Fragment 1
  EXPECT_EQ(contig2.sequences.first[0].nullFeatureId, 0);
  EXPECT_EQ(contig2.sequences.first[0].sequence,
            Dna4Sequence("GTGTCGGAGG"_dna4));
  std::vector contig2Frag1Annots = {gffToken{4, 0, 6}, gffToken{5, 1, 10}};
  ASSERT_THAT(contig2.sequences.first[0].annotations,
              ElementsAreArray(contig2Frag1Annots));
  // Contig 2, Forward, Fragment 2
  EXPECT_EQ(contig2.sequences.first[1].nullFeatureId, 0);
  EXPECT_EQ(contig2.sequences.first[1].sequence,
            Dna4Sequence("CTCCATCGAC"_dna4));
  EXPECT_EQ(contig2.sequences.first[1].annotations.size(), 0);

  // Contig 2, Reverse
  EXPECT_EQ(contig2.sequences.second.size(), 2);
  // Contig 2, Reverse, Fragment 1
  EXPECT_EQ(contig2.sequences.second[0].nullFeatureId, 0);
  EXPECT_EQ(contig2.sequences.second[0].sequence,
            Dna4Sequence("GTCGATGGAG"_dna4));
  EXPECT_EQ(contig2.sequences.second[0].annotations.size(), 0);
  // Contig 2, Reverse, Fragment 2
  EXPECT_EQ(contig2.sequences.second[1].nullFeatureId, 0);
  EXPECT_EQ(contig2.sequences.second[1].sequence,
            Dna4Sequence("CCTCCGACAC"_dna4));
  std::vector contig2Frag2RevAnnots = {gffToken{6, 0, 10}};
  ASSERT_THAT(contig2.sequences.second[1].annotations,
              ElementsAreArray(contig2Frag2RevAnnots));
}
