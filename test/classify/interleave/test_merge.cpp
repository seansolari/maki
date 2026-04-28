
#include "maki/classify/interleave.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/utils/tempfile.hpp"
#include "test_common.hpp"
#include <functional>

using ClassifyFn = std::function<std::vector<int64_t>(
    const DeBruijnGraph *, const DeBruijnGraph *, size_t)>;

class BossClassificationTest : public ::testing::TestWithParam<ClassifyFn> {
protected:
  BossClassificationTest()
      : bufferBase(tempio::create_temporary_directory()),
        refBase(bufferBase / "reference"), qryBase(bufferBase / "query") {
    fs::create_directories(refBase);
    fs::create_directories(qryBase);
  }

  ~BossClassificationTest() { fs::remove_all(bufferBase); }

  ClassifyFn classify;
  fs::path bufferBase, refBase, qryBase;
  std::size_t k;
  ColouredGraph query, reference;

  void SetUp() override { classify = GetParam(); }

  void MakeAndLoadQuery(const std::vector<std::string> &kmers) {
    MakeGraph(kmers, k - 1, 1, qryBase);
    ColouredGraph::FromDisk(query, qryBase);
  }

  void MakeAndLoadReference(const std::vector<std::string> &kmers) {
    MakeGraph(kmers, k - 1, 1, refBase);
    ColouredGraph::FromDisk(reference, refBase);
  }

  void MakeAndLoadQueryFasta(const std::string &sequence) {
    MakeGraph(">header\n" + sequence + "\n", k - 1, 1, qryBase);
    ColouredGraph::FromDisk(query, qryBase);
  }

  void MakeAndLoadReferenceFasta(const std::string &sequence) {
    MakeGraph(">header\n" + sequence + "\n", k - 1, 1, refBase);
    ColouredGraph::FromDisk(reference, refBase);
  }

  void assert_basic_invariants(const std::vector<int64_t> &out) {
    ASSERT_EQ(out.size(), query.edges());

    for (size_t i = 0; i < out.size(); ++i) {
      if (out[i] != -1) {
        ASSERT_GE(out[i], 0);
        ASSERT_LT(static_cast<size_t>(out[i]), reference.edges());
        EXPECT_EQ(query.kmer(i), reference.kmer(out[i]));
      }
    }
  }
};

TEST_P(BossClassificationTest, MinimalBiologicalGraph) {
  k = 3;
  MakeAndLoadQuery({"ACG"});
  MakeAndLoadReference({"ACG"});

  auto out = classify(&query, &reference, /*grainsize=*/1);
  assert_basic_invariants(out);
}

TEST_P(BossClassificationTest, IdenticalGraphs) {
  k = 4;
  std::vector<std::string> kmers = {"ACGT", "CGTA", "GTAC", "TACG"};

  MakeAndLoadQuery(kmers);
  MakeAndLoadReference(kmers);

  auto out = classify(&query, &reference, 4);

  ASSERT_EQ(out.size(), query.edges());
  for (size_t i = 0; i < out.size(); ++i) {
    EXPECT_EQ(out[i], static_cast<int64_t>(i));
  }
}

TEST_P(BossClassificationTest, SameKmersDifferentConstructionOrder) {
  k = 4;
  MakeAndLoadQuery({"ACGT", "CGTA", "GTAC"});
  MakeAndLoadReference({"GTAC", "ACGT", "CGTA"});

  auto out = classify(&query, &reference, 2);
  assert_basic_invariants(out);
}

TEST_P(BossClassificationTest, PartialOverlap) {
  k = 4;
  MakeAndLoadQuery({"ACGT", "CGTA", "TTTT"});
  MakeAndLoadReference({"ACGT", "CGTA"});

  auto out = classify(&query, &reference, 1);
  assert_basic_invariants(out);

  // Expect at least one -1 because "TTTT" is absent biologically
  EXPECT_TRUE(
      std::any_of(out.begin(), out.end(), [](int64_t x) { return x == -1; }));
}

TEST_P(BossClassificationTest, DisjointGraphs) {
  k = 4;
  MakeAndLoadQuery({"AAAA", "TTTT"});
  MakeAndLoadReference({"CCCC", "GGGG"});

  auto out = classify(&query, &reference, 8);

  for (auto x : out) {
    EXPECT_EQ(x, -1);
  }
}

TEST_P(BossClassificationTest, KmerSizeVariation) {
  k = 6;
  MakeAndLoadQuery({"ACGTGA", "CGTGAC"});
  MakeAndLoadReference({"ACGTGA", "CGTGAC"});

  auto out_k6 = classify(&query, &reference, 1);
  assert_basic_invariants(out_k6);
}

TEST_P(BossClassificationTest, GrainSizeInvariance) {
  k = 4;
  MakeAndLoadQuery({"ACGT", "CGTA", "GTAC", "TACG"});
  MakeAndLoadReference({"ACGT", "CGTA", "GTAC", "TACG"});

  auto out1 = classify(&query, &reference, 1);
  auto out2 = classify(&query, &reference, 2);
  auto out3 = classify(&query, &reference, query.edges() * 4);

  EXPECT_EQ(out1, out2);
  EXPECT_EQ(out1, out3);
}

TEST_P(BossClassificationTest, LargeIdenticalGenomeStressTest) {
  const size_t genome_length = 200000;
  auto genome = random_dna_sequence(genome_length, /*seed=*/42);

  k = 31;
  MakeAndLoadQueryFasta(genome);
  MakeAndLoadReferenceFasta(genome);

  auto out = classify(&query, &reference, /*grainsize=*/1024);
  ASSERT_EQ(out.size(), query.edges());

  for (size_t i = 0; i < out.size(); ++i) {
    EXPECT_EQ(out[i], static_cast<int64_t>(i));
  }
}

TEST_P(BossClassificationTest, LargeGenomeSingleSNPSensitivity) {
  const size_t genome_length = 200000;

  auto ref_genome = random_dna_sequence(genome_length, 1337);
  auto mut_genome = apply_snp(ref_genome, genome_length / 2);

  k = 31;
  MakeAndLoadReferenceFasta(ref_genome);
  MakeAndLoadQueryFasta(mut_genome);

  auto out = classify(&query, &reference, 2048);
  assert_basic_invariants(out);

  // Count unmatched edges
  size_t missing = std::count(out.begin(), out.end(), -1);

  // A single SNP disrupts at most k kmers biologically,
  // but BOSS padding adds a few more — we expect O(k) failures
  EXPECT_GT(missing, 0u);
  EXPECT_LT(missing, 4 * k);
}

TEST_P(BossClassificationTest, LargeGenomeClusteredSNPSensitivity) {
  const size_t genome_length = 200000;

  auto genome = random_dna_sequence(genome_length, 9001);

  auto mutated = genome;
  for (size_t i = 0; i < 10; ++i) {
    mutated = apply_snp(mutated, genome_length / 2 + i);
  }

  k = 31;
  MakeAndLoadReferenceFasta(genome);
  MakeAndLoadQueryFasta(mutated);

  auto out = classify(&query, &reference, 4096);
  assert_basic_invariants(out);

  size_t missing = std::count(out.begin(), out.end(), -1);
  EXPECT_GT(missing, k);
  EXPECT_LT(missing, 20 * k);
}

TEST_P(BossClassificationTest, LargeGraphGrainSizeInvariance) {
  const size_t genome_length = 150000;
  auto genome = random_dna_sequence(genome_length, 31415);

  k = 31;
  MakeAndLoadQueryFasta(genome);
  MakeAndLoadReferenceFasta(genome);

  auto out1 = classify(&query, &reference, 1);
  auto out2 = classify(&query, &reference, 512);
  auto out3 = classify(&query, &reference, query.edges());

  EXPECT_EQ(out1, out2);
  EXPECT_EQ(out1, out3);
}

INSTANTIATE_TEST_SUITE_P(BossClassificationImplementations,
                         BossClassificationTest,
                         ::testing::Values(ClassifySmall, ClassifyLarge));
