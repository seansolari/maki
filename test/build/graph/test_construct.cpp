#include "maki/build/graph/build_colours.hpp"
#include "maki/build/graph/construct_cdbg.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/utils/tempfile.hpp"
#include "test_common.hpp"
#include "gmock/gmock.h"

/*

NNN A 1 0 1 *
NNN C 2 0 1 *
NNN G 3 0 1 *
NNN T 4 0 1 *
NNA G 3 0 1 *
ACA C 2 1 1
TCA N 0 0 1
NGA C 2 0 1 *
CGA G 3 2 1
TGA G 3 3 1
NTA C 2 0 1 *
GTA N 0 0 1
NNC G 3 0 1 *
CAC T 4 1 1
GAC T 4 3 1
TAC A 1 1 1
TAC T 4 2 1
CTC A 1 3 1
CTC G 3 2 1
GTC N 0 0 1
NNG A 1 0 1 *
NAG T 4 0 1 *
GAG T 4 2 1
GAG T 4 3 0
NCG A 1 0 1 *
TCG N 0 0 1
NTG A 1 0 1 *
GTG T 4 1 1
NNT A 1 0 1 *
NNT G 3 0 1 *
ACT C 2 2 1
ACT C 2 3 0
AGT A 1 2 1
AGT C 2 3 1
AGT G 3 1 1
TGT A 1 1 1

*/
class EgidiGraphTests : public testing::Test {
protected:
  EgidiGraphTests()
      : k(3), c(), genomes(3), out(tempio::create_temporary_directory()) {
    ParseFastaToGenome(genomes[0], STRING(CONTIG_ONE_PATH), c, k);
    ParseFastaToGenome(genomes[1], STRING(CONTIG_TWO_PATH), c, k);
    ParseFastaToGenome(genomes[2], STRING(CONTIG_THREE_PATH), c, k);
  }

  ~EgidiGraphTests() { std::filesystem::remove_all(out); }

  uint8_t k;
  Colours c;
  std::vector<Dna4Genome> genomes;
  std::filesystem::path out;
};

TEST_F(EgidiGraphTests, SuffixFill_1) {
  // Construct graph
  cdbg::construct(toView(genomes), MetaColours(std::move(c.ids)),
                  {.kmer_size = k, .suffix_size = 1, .out = out});

  // Load graph
  ColouredGraph g;
  ColouredGraph::FromDisk(g, out);

  // Structure
  EXPECT_EQ(g.k, k);
  EXPECT_THAT(g.C, ElementsAreArray({1, 8, 6, 7, 4}));
  EXPECT_THAT(g.F, ElementsAreArray({4, 8, 8, 7, 7}));

  // Edges
  wavelet_matrix XW = EdgesToWaveletMatrix(
      {0b1001, 0b1010, 0b1011, 0b1100, 0b1011, 0b1010, 0b1000, 0b1010, 0b1011,
       0b0011, 0b1010, 0b1000, 0b1011, 0b1100, 0b0100, 0b1001, 0b0100, 0b1001,
       0b1011, 0b1000, 0b1001, 0b1100, 0b0100, 0b1001, 0b1000, 0b1001, 0b1100,
       0b1001, 0b1011, 0b1010, 0b1001, 0b1010, 0b1011, 0b0001});
  sdsl::bit_vector Xl = {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1,
                         0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1};

  EXPECT_THAT(g.W, ContainerEq(XW));
  EXPECT_THAT(g.l, ContainerEq(Xl));

  // Colours
  auto craw = UnpackColours(g);
  EXPECT_THAT(
      craw.first,
      ElementsAreArray({0, 0, 0, 0, 0, 1, 0, 0, 2, 3, 0, 0, 0, 1, 3, 1, 2, 3,
                        2, 0, 0, 0, 2, 3, 0, 0, 0, 1, 0, 0, 2, 3, 2, 3, 1, 1}));
  EXPECT_THAT(
      craw.second,
      ElementsAreArray({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1}));
}

// larger k-mer and suffix size

class SmallFastaTests : public testing::Test {
protected:
  SmallFastaTests()
      : k(9), c(), genomes(1), out(tempio::create_temporary_directory()) {
    ParseFastaToGenome(genomes[0], STRING(SMALL_DUP_SEQ), c, k);
  }

  ~SmallFastaTests() { std::filesystem::remove_all(out); }

  uint8_t k;
  Colours c;
  std::vector<Dna4Genome> genomes;
  std::filesystem::path out;
};

TEST_F(SmallFastaTests, SuffixFill_3) {
  std::size_t s = 3;

  // Construct graph
  cdbg::construct(toView(genomes), MetaColours(std::move(c.ids)),
                  {.kmer_size = k, .suffix_size = s, .out = out});

  // Load graph
  ColouredGraph g;
  ColouredGraph::FromDisk(g, out);

  // Structure
  EXPECT_EQ(g.k, k);
  EXPECT_THAT(g.C, ElementsAreArray({1, 27, 33, 33, 27}));
  EXPECT_THAT(g.F, ElementsAreArray({2, 27, 33, 33, 27}));

  // Edges
  wavelet_matrix XW = EdgesToWaveletMatrix(
      {0b1010, 0b1011, 0b1011, 0b1011, 0b1010, 0b1100, 0b1010, 0b1011, 0b1010,
       0b1100, 0b1010, 0b1100, 0b1100, 0b1011, 0b1001, 0b1100, 0b1001, 0b1011,
       0b1010, 0b1010, 0b1100, 0b1011, 0b1001, 0b1011, 0b1011, 0b1001, 0b1001,
       0b1011, 0b1010, 0b1001, 0b1011, 0b1100, 0b1000, 0b1010, 0b1001, 0b1001,
       0b1100, 0b1011, 0b1100, 0b1001, 0b1011, 0b1001, 0b1001, 0b1001, 0b1100,
       0b1011, 0b1010, 0b1100, 0b1100, 0b1100, 0b1001, 0b1100, 0b1010, 0b1011,
       0b1010, 0b1001, 0b1011, 0b1010, 0b1011, 0b1011, 0b1100, 0b1010, 0b1100,
       0b1001, 0b1100, 0b1010, 0b1010, 0b1010, 0b1011, 0b1010, 0b1010, 0b1100,
       0b1001, 0b1010, 0b1001, 0b1011, 0b1001, 0b1100, 0b1001, 0b1011, 0b1010,
       0b1100, 0b1001, 0b1001, 0b1001, 0b1010, 0b1001, 0b1011, 0b1100, 0b1011,
       0b1010, 0b1000, 0b1100, 0b1011, 0b1010, 0b1011, 0b1011, 0b1010, 0b1011,
       0b1100, 0b1011, 0b1100, 0b1001, 0b1010, 0b1010, 0b1011, 0b1010, 0b1010,
       0b1011, 0b1100, 0b1011, 0b1100, 0b1001, 0b1100, 0b1011, 0b1010, 0b1010,
       0b1001, 0b1010, 0b1011, 0b1001, 0b1010});
  sdsl::bit_vector Xl = {
      0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  EXPECT_THAT(g.W, ContainerEq(XW));
  EXPECT_THAT(g.l, ContainerEq(Xl));

  // Colours
  auto craw = UnpackColours(g);
  EXPECT_THAT(
      craw.first,
      ElementsAreArray({0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
                        1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
                        1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1,
                        1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1}));
  EXPECT_THAT(craw.second, ElementsAreArray(std::vector<uint64_t>(122, 1)));
}
