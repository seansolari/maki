
#include "maki/build/graph/build_colours.hpp"
#include "maki/build/io/filter.hpp"
#include "maki/maki.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>

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
