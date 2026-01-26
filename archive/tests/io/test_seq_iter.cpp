#include <ostream>
#include <vector>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <maki/maki.h>
#include <maki/fasta.hpp>

using ::testing::ElementsAreArray;

#define SV(...) std::vector<size_t>({ __VA_ARGS__ })

constexpr size_t minContigSize = 3;

std::vector<std::vector<uint64_t>> unwindColours(const AnnotatedSequence &frag, uint8_t windowSize) {
    std::vector<std::vector<uint64_t>> annotTracker;
    for (auto it = frag.colouredSegmentsEdgeBegin(windowSize);
         it != frag.colouredSegmentsEdgeEnd();
         ++it)
    {
        auto &currentColours = annotTracker.emplace_back();

        if (!it.positionIsAnnotated())
        {
            currentColours.push_back(it.getNullFeatureId());
        }
        else
        {
            for (const auto& annot : it.viewCurrentAnnots())
            {
                currentColours.push_back(annot.featureId);
            }
        }
    }
    return annotTracker;
}

struct _SegmentToken {
    size_t beginPos, endPos;
    uint64_t featureId;
    bool endIsTerminal;

    inline bool operator==(const _SegmentToken&) const =default;

    inline bool operator!=(const _SegmentToken&) const =default;

    friend std::ostream& operator<<(std::ostream &os, const _SegmentToken &obj) {
        return os
            << "SegmentToken("
            << obj.beginPos << ", "
            << obj.endPos << ", "
            << obj.featureId << ", "
            << (obj.endIsTerminal ? "yes" : "no") << ")";
    }

    friend void PrintTo(const _SegmentToken& obj, std::ostream *os) { (*os) << obj; return; }
};

_SegmentToken operator-(const AnnotatedSequenceSegment &lhs, Dna4SequenceConstIter rhs) {
    return { static_cast<size_t>(lhs.seqBegin - rhs), static_cast<size_t>(lhs.seqEnd - rhs), lhs.featureId, lhs.endIsTerminal };
}

std::vector<_SegmentToken> unwindSegments(const AnnotatedSequence &frag, uint8_t windowSize) {
    std::vector<_SegmentToken> observedTokens;
    for (auto it = frag.colouredSegmentsBegin(windowSize);
         it != frag.colouredSegmentsEnd();
         ++it)
    {
        observedTokens.emplace_back(*it - frag.sequence.cbegin());
    }
    return observedTokens;
}

class SegmentIteratorTests : public testing::Test {
protected:
    SegmentIteratorTests() : genomes(loadGenomes({ STRING(SEQA_GFF), STRING(SEQB_GFF) }, minContigSize)) {}

    Dna4GenomeVector genomes;
};

// remember that annotation IDs are established based on alphabetical order of annotations in GFF file

// seqA tests

/*

Contig 1
========

Annotations
-----------
0b0101 -> genome1	1	6	+
0b1001 -> genome1	17	24	+
0b1101 -> genome1	8	14	-

Forward Sequence
----------------

1 1 1 1 1
                                2 2 2 2 2 2 2

A C G A G N C G G T G G C A A N G A A G T T N

0                   1                   2
0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2

Reverse Sequence
----------------

              3 3 3 3 3 3

N A A C T T C N T T G C C A C C G N C T C G T

0                   1                   2
0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2

*/

TEST_F(SegmentIteratorTests, SeqAContig1SegmentEdges) {
    const Dna4Contig &contig = genomes[0].contigs[0];

    // forward sequences

    std::vector<std::vector<size_t>>
        observedColours = unwindColours(contig.sequences.first[0], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b0101 ), SV( 0b0101 ) }));

    observedColours = unwindColours(contig.sequences.first[1], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ) }));

    observedColours = unwindColours(contig.sequences.first[2], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b1001 ), SV( 0b1001 ), SV( 0b1001 ) }));

    // reverse sequences

    observedColours = unwindColours(contig.sequences.second[0], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ) }));

    observedColours = unwindColours(contig.sequences.second[1], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b1101 ), SV( 0b1101 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ) }));

    observedColours = unwindColours(contig.sequences.second[2], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ), SV( 0b0001 ) }));
}

TEST_F(SegmentIteratorTests, SeqAContig1Segments) {
    const Dna4Contig &contig = genomes[0].contigs[0];

    // forward sequences

    std::vector<_SegmentToken>
        observedTokens = unwindSegments(contig.sequences.first[0], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 5, 0b0101, false }, _SegmentToken{ 2, 5, 0b0001, true } }));

    observedTokens = unwindSegments(contig.sequences.first[1], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 9, 0b0001, true } }));

    observedTokens = unwindSegments(contig.sequences.first[2], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 6, 0b1001, false }, _SegmentToken{ 3, 6, 0b0001, true } }));

    // reverse sequences

    observedTokens = unwindSegments(contig.sequences.second[0], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 6, 0b0001, true } }));

    observedTokens = unwindSegments(contig.sequences.second[1], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 5, 0b1101, false }, _SegmentToken{ 2, 9, 0b0001, true } }));

    observedTokens = unwindSegments(contig.sequences.second[2], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 5, 0b0001, true } }));
}

/*

Contig 2
========

Annotations
-----------
0b00001
0b10001 -> genome2	1	7	+
0b10101 -> genome2	2	11	+
0b11001 -> genome2	8	22	-

Forward Sequence
----------------

4 4 4 4 4 4
  5 5 5 5 5 5 5 5 5

G T G T C G G A G G N C T C C A T C G A C

0                   1                   2
0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0

Reverse Sequence
----------------

              6 6 6 6 6 6 6 6 6 6 6 6 6 6

G T C G A T G G A G N C C T C C G A C A C

0                   1                   2
0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0

*/

TEST_F(SegmentIteratorTests, SeqAContig2SegmentEdges) {
    const Dna4Contig &contig = genomes[0].contigs[1];

    // forward sequences

    std::vector<std::vector<size_t>>
        observedColours = unwindColours(contig.sequences.first[0], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b10001 ), SV( 0b10001, 0b10101 ), SV( 0b10001, 0b10101 ), SV( 0b10101 ), SV( 0b10101 ), SV( 0b10101 ), SV( 0b10101 ) }));

    observedColours = unwindColours(contig.sequences.first[1], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ) }));

    // reverse sequences

    observedColours = unwindColours(contig.sequences.second[0], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ) }));

    observedColours = unwindColours(contig.sequences.second[1], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b00001 ), SV( 0b00001 ), SV( 0b00001 ), SV( 0b11001 ), SV( 0b11001 ), SV( 0b11001 ), SV( 0b11001 ), SV( 0b11001 ), SV( 0b11001 ), SV( 0b11001 ) }));
}

TEST_F(SegmentIteratorTests, SeqAContig2Segments) {
    const Dna4Contig &contig = genomes[0].contigs[1];

    // forward sequences

    std::vector<_SegmentToken>
        observedTokens = unwindSegments(contig.sequences.first[0], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 6, 0b10001, false }, _SegmentToken{ 1, 10, 0b10101, false }, _SegmentToken{ 7, 10, 0b00001, true } }));

    observedTokens = unwindSegments(contig.sequences.first[1], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 10, 0b00001, true } }));

    // reverse sequences

    observedTokens = unwindSegments(contig.sequences.second[0], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 10, 0b00001, true } }));

    observedTokens = unwindSegments(contig.sequences.second[1], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 10, 0b11001, false }, _SegmentToken{ 7, 10, 0b00001, true } }));
}

// seqB tests

/*

Contig 1
========

Annotations
-----------

0b00010
0b01110 -> contig_one	0	12	+
0b10010 -> contig_one	0	12	-

Forward Sequence
----------------

T A C A C T T A C T C G

0                   1
0 1 2 3 4 5 6 7 8 9 0 1

Reverse Sequence
----------------

C G A G T A A G T G T A

0                   1
0 1 2 3 4 5 6 7 8 9 0 1

*/

TEST_F(SegmentIteratorTests, SeqBContig1SegmentEdges) {
    const Dna4Contig &contig = genomes[1].contigs[0];

    // forward sequences

    std::vector<std::vector<size_t>>
        observedColours = unwindColours(contig.sequences.first[0], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b0010 ), SV( 0b0010 ), SV( 0b0010 ), SV( 0b01110 ), SV( 0b01110 ), SV( 0b01110 ), SV( 0b01110 ), SV( 0b01110 ), SV( 0b01110 ), SV( 0b01110 ), SV( 0b01110 ), SV( 0b01110 ) }));

    // reverse sequences

    observedColours = unwindColours(contig.sequences.second[0], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b0010 ), SV( 0b0010 ), SV( 0b0010 ), SV( 0b10010 ), SV( 0b10010 ), SV( 0b10010 ), SV( 0b10010 ), SV( 0b10010 ), SV( 0b10010 ), SV( 0b10010 ), SV( 0b10010 ), SV( 0b10010 ) }));
}

TEST_F(SegmentIteratorTests, SeqBContig1Segments) {
    const Dna4Contig &contig = genomes[1].contigs[0];

    // forward sequences

    std::vector<_SegmentToken>
        observedTokens = unwindSegments(contig.sequences.first[0], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 12, 0b01110, false }, _SegmentToken{ 9, 12, 0b0010, true } }));

    // reverse sequences

    observedTokens = unwindSegments(contig.sequences.second[0], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 12, 0b10010, false }, _SegmentToken{ 9, 12, 0b0010, true } }));
}

/*

Contig 2
========

Annotations
-----------

0b000010
0b011110 -> contig_two	0	12	+
0b100010 -> contig_two	0	12	-

Forward Sequence
----------------

T A C T C G G A C T C A

0                   1
0 1 2 3 4 5 6 7 8 9 0 1

Reverse Sequence
----------------

T G A G T C C G A G T A

0                   1
0 1 2 3 4 5 6 7 8 9 0 1

*/

TEST_F(SegmentIteratorTests, SeqBContig2SegmentEdges) {
    const Dna4Contig &contig = genomes[1].contigs[1];

    // forward sequences

    std::vector<std::vector<size_t>>
        observedColours = unwindColours(contig.sequences.first[0], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b000010 ), SV( 0b000010 ), SV( 0b000010 ), SV( 0b011110 ), SV( 0b011110 ), SV( 0b011110 ), SV( 0b011110 ), SV( 0b011110 ), SV( 0b011110 ), SV( 0b011110 ), SV( 0b011110 ), SV( 0b011110 ) }));

    // reverse sequences

    observedColours = unwindColours(contig.sequences.second[0], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b000010 ), SV( 0b000010 ), SV( 0b000010 ), SV( 0b100010 ), SV( 0b100010 ), SV( 0b100010 ), SV( 0b100010 ), SV( 0b100010 ), SV( 0b100010 ), SV( 0b100010 ), SV( 0b100010 ), SV( 0b100010 ) }));
}

TEST_F(SegmentIteratorTests, SeqBContig2Segments) {
    const Dna4Contig &contig = genomes[1].contigs[1];

    // forward sequences

    std::vector<_SegmentToken>
        observedTokens = unwindSegments(contig.sequences.first[0], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 12, 0b011110, false }, _SegmentToken{ 9, 12, 0b000010, true } }));

    // reverse sequences

    observedTokens = unwindSegments(contig.sequences.second[0], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 12, 0b100010, false }, _SegmentToken{ 9, 12, 0b000010, true } }));
}

/*

Contig 3
========

Annotations
-----------

0b00010
0b10110 -> contig_three	0	12	+
0b11010 -> contig_three	0	12	-

Forward Sequence
----------------

G A C T C A G A C T C A

0                   1
0 1 2 3 4 5 6 7 8 9 0 1

Reverse Sequence
----------------

T G A G T C T G A G T C

0                   1
0 1 2 3 4 5 6 7 8 9 0 1

*/

TEST_F(SegmentIteratorTests, SeqBContig3SegmentEdges) {
    const Dna4Contig &contig = genomes[1].contigs[2];

    // forward sequences

    std::vector<std::vector<size_t>>
        observedColours = unwindColours(contig.sequences.first[0], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b00010 ), SV( 0b00010 ), SV( 0b00010 ), SV( 0b10110 ), SV( 0b10110 ), SV( 0b10110 ), SV( 0b10110 ), SV( 0b10110 ), SV( 0b10110 ), SV( 0b10110 ), SV( 0b10110 ), SV( 0b10110 ) }));

    // reverse sequences

    observedColours = unwindColours(contig.sequences.second[0], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b00010 ), SV( 0b00010 ), SV( 0b00010 ), SV( 0b11010 ), SV( 0b11010 ), SV( 0b11010 ), SV( 0b11010 ), SV( 0b11010 ), SV( 0b11010 ), SV( 0b11010 ), SV( 0b11010 ), SV( 0b11010 ) }));
}

TEST_F(SegmentIteratorTests, SeqBContig3Segments) {
    const Dna4Contig &contig = genomes[1].contigs[2];

    // forward sequences

    std::vector<_SegmentToken>
        observedTokens = unwindSegments(contig.sequences.first[0], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 12, 0b10110, false }, _SegmentToken{ 9, 12, 0b00010, true } }));

    // reverse sequences

    observedTokens = unwindSegments(contig.sequences.second[0], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 12, 0b11010, false }, _SegmentToken{ 9, 12, 0b00010, true } }));
}

/*

Contig 4
========

Annotations
-----------

0b0010
0b0110 -> contig_four	0	11	+
0b1010 -> contig_four	0	11	-

Forward Sequence
----------------

15 -->

A C G N A C N A C G T

0                   1
0 1 2 3 4 5 6 7 8 9 0

Reverse Sequence
----------------

17 -->

A C G T N G T N C G T

0                   1
0 1 2 3 4 5 6 7 8 9 0

*/

TEST_F(SegmentIteratorTests, SeqBContig4SegmentEdges) {
    const Dna4Contig &contig = genomes[1].contigs[3];

    // forward sequences

    std::vector<std::vector<size_t>>
        observedColours = unwindColours(contig.sequences.first[0], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b0010 ), SV( 0b0010 ), SV( 0b0010 ) }));

    observedColours = unwindColours(contig.sequences.first[1], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b0010 ), SV( 0b0010 ), SV( 0b0010 ), SV( 0b0110 ) }));

    // reverse sequences

    observedColours = unwindColours(contig.sequences.second[0], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b0010 ), SV( 0b0010 ), SV( 0b0010 ), SV( 0b1010 ) }));

    observedColours = unwindColours(contig.sequences.second[1], minContigSize);
    ASSERT_THAT(observedColours, ElementsAreArray({ SV( 0b0010 ), SV( 0b0010 ), SV( 0b0010 ) }));
}

TEST_F(SegmentIteratorTests, SeqBContig4Segments) {
    const Dna4Contig &contig = genomes[1].contigs[3];

    // forward sequences

    std::vector<_SegmentToken>
        observedTokens = unwindSegments(contig.sequences.first[0], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 3, 0b0010, true } }));

    observedTokens = unwindSegments(contig.sequences.first[1], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 4, 0b0110, false }, _SegmentToken{ 1, 4, 0b0010, true } }));

    // reverse sequences

    observedTokens = unwindSegments(contig.sequences.second[0], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 4, 0b1010, false }, _SegmentToken{ 1, 4, 0b0010, true } }));

    observedTokens = unwindSegments(contig.sequences.second[1], minContigSize);
    ASSERT_THAT(observedTokens, ElementsAreArray({ _SegmentToken{ 0, 3, 0b0010, true } }));
}
