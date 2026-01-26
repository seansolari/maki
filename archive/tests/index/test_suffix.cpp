#include <cstddef>
#include <vector>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <maki/fasta.hpp>
#include <maki/mers.hpp>
#include <maki/reads.hpp>

using ::testing::ElementsAre;
using ::testing::ElementsAreArray;

using namespace seqan3::literals;

constexpr size_t k = 3;

class CountFnaSuffixTests : public testing::Test {
protected:
    CountFnaSuffixTests() : genomes(loadGenomes({ STRING(SEQB_FNA) }, k)) {}

    Dna4GenomeVector genomes;
};

TEST_F(CountFnaSuffixTests, Contig1S1Forward) {
    uint8_t s = 1;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[0].sequences.first[0];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 2, 4, 1, 3 ));
}

TEST_F(CountFnaSuffixTests, Contig1S1Reverse) {
    uint8_t s = 1;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[0].sequences.second[0];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 4, 0, 3, 3 ));
}

TEST_F(CountFnaSuffixTests, Contig1S2Forward) {
    uint8_t s = 2;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[0].sequences.first[0];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 0, 1, 0, 1, 3, 0, 0, 1, 0, 1, 0, 0, 0, 2, 0, 1 ));
}

TEST_F(CountFnaSuffixTests, Contig1S2Reverse) {
    uint8_t s = 2;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[0].sequences.second[0];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 1, 0, 1, 2, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 3, 0 ));
}

class CountGffSuffixTests : public testing::Test {
protected:
    CountGffSuffixTests() : genomes(loadGenomes({ STRING(SEQA_GFF) }, k)) {}

    Dna4GenomeVector genomes;
};

TEST_F(CountGffSuffixTests, Contig1F1S1Forward) {
    uint8_t s = 1;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[0].sequences.first[0];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 1, 0, 2, 0 ));
}

TEST_F(CountGffSuffixTests, Contig1F2S1Forward) {
    uint8_t s = 1;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[0].sequences.first[1];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 2, 1, 3, 1 ));
}

TEST_F(CountGffSuffixTests, Contig1F3S1Forward) {
    uint8_t s = 1;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[0].sequences.first[2];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 1, 0, 1, 2 ));
}

TEST_F(CountGffSuffixTests, Contig1F1S1Reverse) {
    uint8_t s = 1;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[0].sequences.second[0];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 0, 2, 0, 2 ));
}

TEST_F(CountGffSuffixTests, Contig1F2S1Reverse) {
    uint8_t s = 1;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[0].sequences.second[1];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 1, 4, 2, 0 ));
}

TEST_F(CountGffSuffixTests, Contig1F3S1Reverse) {
    uint8_t s = 1;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[0].sequences.second[2];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 0, 1, 1, 1 ));
}

TEST_F(CountGffSuffixTests, Contig2F1S1Forward) {
    uint8_t s = 1;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[1].sequences.first[0];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 1, 1, 5, 1 ));
}

TEST_F(CountGffSuffixTests, Contig2F1S2Forward) {
    uint8_t s = 2;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[1].sequences.first[0];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 2, 1, 0, 0, 1, 0 ));
}

TEST_F(CountGffSuffixTests, Contig2F2S1Forward) {
    uint8_t s = 1;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[1].sequences.first[1];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 2, 4, 1, 1 ));
}

TEST_F(CountGffSuffixTests, Contig2F2S2Forward) {
    uint8_t s = 2;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[1].sequences.first[1];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 0, 1, 1, 0, 1, 1, 0, 2, 0, 1, 0, 0, 1, 0, 0, 0 ));
}

TEST_F(CountGffSuffixTests, Contig2F1S1Reverse) {
    uint8_t s = 1;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[1].sequences.second[0];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 2, 1, 4, 1 ));
}

TEST_F(CountGffSuffixTests, Contig2F1S2Reverse) {
    uint8_t s = 2;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[1].sequences.second[0];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 0, 0, 2, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0 ));
}

TEST_F(CountGffSuffixTests, Contig2F2S1Reverse) {
    uint8_t s = 1;
    suffix::SuffixTable table(s);
    
    auto &fragment = genomes[0].contigs[1].sequences.second[1];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 2, 4, 1, 1 ));
}

TEST_F(CountGffSuffixTests, Contig2F2S2Reverse) {
    uint8_t s = 2;
    suffix::SuffixTable table(s);

    auto &fragment = genomes[0].contigs[1].sequences.second[1];
    table.count(fragment.sequence.cbegin() + k - s, fragment.sequence.cend());

    ASSERT_THAT(table.cdata(), ElementsAre( 0, 1, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0 ));
}

TEST_F(CountGffSuffixTests, CombinedNoAnnot) {
    constexpr uint8_t s = 1;
    auto counts = suffix::countSuffixes(genomes, UnwindGenome{}, s, k - s);
    ASSERT_EQ(counts.size(), 1);
    ASSERT_THAT(counts[0].cdata(), ElementsAre( 12, 18, 20, 10 ));
}

class ChunkedGenomeSuffixTests : public testing::Test {
protected:
    ChunkedGenomeSuffixTests()
        : xgenome(std::move(loadGenomes({ STRING(SEQA_FNA) }, k)[0]), 2, k)
    {
    }
    
    ChunkedDna4Genome xgenome;
};

TEST_F(ChunkedGenomeSuffixTests, CountChunkedGenome) {
    constexpr uint8_t s = 1;
    auto counts = suffix::accumulate(suffix::countSuffixes(xgenome, UnwindChunkedGenome{}, s, k - s));
    ASSERT_THAT(counts.cdata(), ElementsAre( 12, 18, 20, 10 ));
}

class ChunkedGenomeChunkTests : public testing::Test {
protected:
    ChunkedGenomeChunkTests()
        : xgenome(makeGenome())
    {
    }
    
    static ChunkedDna4Genome makeGenome() {
        std::string data = ">seq1\nACGTACGTACGTACGTACGTACGTACGATCAGTCAGTCAGTCGTA\n>seq2\nAGTACGTCGTACGATCGTC\n>seq3\nATCGTCGATGCTAGCTAGCTAGCTAGCTACGT\n>seq4\nGTACGTGCTAGCTAGCTGACTCGATGCATTAA\n>seq5\nTAGTATATATAGTAGTAGT\n";
        std::istringstream datastream(data);

        Dna4Genome genome;
        parseFastaStream(genome, datastream, 0, k);

        return ChunkedDna4Genome(std::move(genome), 2, k);
    }

    ChunkedDna4Genome xgenome;
};

TEST_F(ChunkedGenomeChunkTests, CountChunks) {
    constexpr uint8_t s = 2;
    auto counts = suffix::countSuffixes(xgenome, UnwindChunkedGenome{}, s, k - s);
    ASSERT_EQ(counts.size(), 18);

    ASSERT_THAT(counts[ 0].cdata(), ElementsAre( 0, 0, 0, 4, 4, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0 ));
    ASSERT_THAT(counts[ 1].cdata(), ElementsAre( 0, 1, 1, 2, 2, 0, 0, 2, 1, 3, 0, 0, 1, 0, 3, 0 ));
    ASSERT_THAT(counts[ 2].cdata(), ElementsAre( 0, 2, 0, 1, 0, 0, 0, 2, 2, 1, 0, 0, 0, 0, 3, 0 ));
    ASSERT_THAT(counts[ 3].cdata(), ElementsAre( 0, 0, 4, 0, 4, 0, 0, 0, 0, 1, 0, 3, 1, 3, 0, 0 ));
    ASSERT_THAT(counts[ 4].cdata(), ElementsAre( 0, 0, 0, 4, 3, 0, 0, 1, 0, 4, 0, 0, 0, 0, 4, 0 ));
    ASSERT_THAT(counts[ 5].cdata(), ElementsAre( 0, 0, 0, 2, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 0 ));
    ASSERT_THAT(counts[ 6].cdata(), ElementsAre( 0, 0, 1, 2, 2, 0, 0, 3, 0, 4, 0, 0, 1, 0, 4, 0 ));
    ASSERT_THAT(counts[ 7].cdata(), ElementsAre( 0, 0, 2, 2, 4, 0, 0, 1, 0, 4, 0, 0, 1, 1, 2, 0 ));
    ASSERT_THAT(counts[ 8].cdata(), ElementsAre( 0, 0, 1, 2, 0, 0, 2, 2, 2, 2, 0, 1, 1, 2, 1, 0 ));
    ASSERT_THAT(counts[ 9].cdata(), ElementsAre( 0, 0, 0, 3, 1, 0, 3, 0, 2, 1, 0, 0, 0, 3, 1, 0 ));
    ASSERT_THAT(counts[10].cdata(), ElementsAre( 0, 0, 0, 4, 0, 0, 3, 0, 4, 1, 0, 0, 0, 3, 1, 0 ));
    ASSERT_THAT(counts[11].cdata(), ElementsAre( 0, 1, 2, 1, 1, 0, 2, 1, 1, 2, 0, 0, 2, 1, 0, 0 ));
    ASSERT_THAT(counts[12].cdata(), ElementsAre( 0, 0, 0, 3, 1, 0, 3, 0, 2, 1, 0, 2, 0, 3, 1, 0 ));
    ASSERT_THAT(counts[13].cdata(), ElementsAre( 1, 1, 2, 1, 1, 0, 1, 1, 0, 1, 0, 1, 2, 1, 0, 1 ));
    ASSERT_THAT(counts[14].cdata(), ElementsAre( 1, 2, 1, 1, 0, 0, 2, 2, 2, 1, 0, 1, 2, 0, 1, 0 ));
    ASSERT_THAT(counts[15].cdata(), ElementsAre( 0, 1, 0, 3, 2, 0, 2, 0, 2, 1, 0, 0, 0, 2, 1, 0 ));
    ASSERT_THAT(counts[16].cdata(), ElementsAre( 0, 0, 0, 6, 0, 0, 0, 0, 4, 0, 0, 0, 3, 0, 4, 0 ));
    ASSERT_THAT(counts[17].cdata(), ElementsAre( 0, 0, 0, 7, 3, 0, 0, 0, 0, 0, 0, 0, 3, 4, 0, 0 ));
}

TEST_F(ChunkedGenomeChunkTests, CountTotal) {
    constexpr uint8_t s = 2;
    auto counts = suffix::accumulate(suffix::countSuffixes(xgenome, UnwindChunkedGenome{}, s, k - s));
    ASSERT_THAT(counts.cdata(), ElementsAre( 2, 8, 14, 48, 31, 0, 18, 15, 22, 34, 0, 8, 17, 23, 33, 1 ));
}

class ChunkedReadsSuffixTests : public testing::Test {
protected:
    ChunkedReadsSuffixTests()
        : _data(5)
    {
        _data[0]
            .emplace_back()
            .emplace_back("NC_000913.3-61888-1")
            .push("GCAACGTGTTCCTAAAGACGTATTTATGGGCGTTGATGAACTGCAGGTAGGTATGCGTTTCCTGGCTGAAACCGACCAGGGTCCGGTACCGGTTGAAATCACTGCGGTTGAAGACGATCACGTCGTGGTTGATGGTAACCACATGCTGGC"_dna5, k);
        _data[0][0]
            .emplace_back("NC_000913.3-61888-2")
            .push("GTCGCTTCGCGAATCGCCACAACTTCAACGTTGAATTTCAGGTTCTGACCGGCCAGCATGTGGTTACCATCAACCACGACGTGATCGTCTTCAACCGCAGTGATTTCAACCGGTACCGGACCCTGGTCGGTTTCAGCCAGGAAACGCATA"_dna5, k);
        _data[1]
            .emplace_back()
            .emplace_back("NC_000913.3-61886")
            .push("GGAACTGGATGTGGTAATCCACATAGCTGACGCCTGCTTTCATGGTGGCGATCAGCGCCAGTGGTTCATCATTAACGTCCTTCACCAGCTGCGCGTAGTCGTTGTCACTTTTTGCCGACCAGGTACGGGTCAGGTCAGCCGCATAGCCGT"_dna5, k);
        _data[2]
            .emplace_back()
            .emplace_back("NC_000913.3-61884")
            .push("GCTTCACTGAAACGGGCTGGCGGCTTGGTAAAGTGCTGGGCTGGTGTAAGTTCAACGAGCGTCAGAGCATCGCCTTTATTAACTGCTGGTAAGATGCGATCTTCATCGCCTTTACGCAACGCAGGCATCACTTTTGTCCAGCCATCAAAA"_dna5, k);
        _data[3]
            .emplace_back()
            .emplace_back("NC_000913.3-61882")
            .push("TGTCCTGAAAACGGATGATTATCGCATAACCATCGATTATCTCCGAACCCGTCCTGTTGATTTAATCATTATGGATATAGACTTGCCCGGAACAGACGGTTTTACCTTCCTGAAAAGGATCAAACAAATCCAGAGCACAGTGAAAGTGTT"_dna5, k);
        _data[4]
            .emplace_back()
            .emplace_back("NC_000913.3-61880")
            .push("TAAAAGACAAACGCGAGGCTAAGACCTCGCGTTTTGCTTTAATCAACCAGATGATATTTTTCTGAAAGCACATGGGCCAGGTGTTTGAACATATTAAACACCGCGGTGCTTTTGGCTGTTGGCAATCCTTGTTCATCTAAAAAGTAGTCG"_dna5, k);
    }

    std::vector<reads::ReadChunks> _data;
};

TEST_F(ChunkedReadsSuffixTests, CountChunkedReads) {
constexpr uint8_t s = 2;
    auto counts = suffix::countSuffixes(_data, reads::UnwindReads{}, s, k - s);
    ASSERT_EQ(counts.size(), _data.size());

    ASSERT_THAT(counts[ 0].cdata(), ElementsAre( 40, 45, 37, 19, 50, 39, 30, 37, 21, 50, 39, 45, 30, 21, 49, 40 ));
    ASSERT_THAT(counts[ 1].cdata(), ElementsAre( 13, 25, 16, 12, 21, 19, 26, 16, 16, 22, 18, 25, 16, 16, 22, 13 ));
    ASSERT_THAT(counts[ 2].cdata(), ElementsAre( 24, 23, 15, 12, 15, 14, 31, 15, 20, 18, 14, 23, 14, 20, 15, 23 ));
    ASSERT_THAT(counts[ 3].cdata(), ElementsAre( 27, 19, 24, 16, 16, 16,  6, 24, 13, 14, 16, 18, 30, 13, 16, 28 ));
    ASSERT_THAT(counts[ 4].cdata(), ElementsAre( 36, 21, 14, 15, 16, 12, 20, 14, 16, 13, 12, 21, 18, 16, 16, 36 ));
}

TEST_F(ChunkedReadsSuffixTests, CountTotal) {
    constexpr uint8_t s = 2;
    auto counts = suffix::accumulate(suffix::countSuffixes(_data, reads::UnwindReads{}, s, k - s));
    ASSERT_THAT(counts.cdata(), ElementsAre( 140, 133, 106, 74, 118, 100, 113, 106, 86, 117, 99, 132, 108, 86, 118, 140 ));
}
