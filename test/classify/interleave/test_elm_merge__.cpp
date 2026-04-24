#include <gtest/gtest.h>

#include "dbg_api.hpp"
#include "elm_api.hpp"
#include "test_common.hpp"

// ------------------------------------------------------------
// Register ELM implementations here
// ------------------------------------------------------------

#include "elm_impl_a.hpp"
#include "elm_impl_b.hpp"

using Implementations = ::testing::Types<ELMImplA, ELMImplB>;

template <typename Impl> class ELMMergeTest : public ::testing::Test {
protected:
  static constexpr size_t k = 3;
};

TYPED_TEST_SUITE(ELMMergeTest, Implementations);

// ------------------------------------------------------------
// Refinement tests (run first!)
// ------------------------------------------------------------

TYPED_TEST(ELMMergeTest, RefinementMonotonicity) {
  if (!TypeParam::HasRefinementTrace())
    return;

  auto g1 = TypeParam::Build(this->k, {"ATG", "TGA"});
  auto g2 = TypeParam::Build(this->k, {"TGA", "GAT"});
  auto trace = TypeParam::Trace(g1, g2);

  size_t rounds = trace.NumRounds();
  for (size_t r = 1; r < rounds; ++r) {
    EXPECT_GE(trace.BlockCount(r), trace.BlockCount(r - 1));
  }
}

TYPED_TEST(ELMMergeTest, RefinementStabilisesByK) {
  if (!TypeParam::HasRefinementTrace())
    return;

  auto g1 = TypeParam::Build(this->k, {"AAA", "AAT"});
  auto g2 = TypeParam::Build(this->k, {"AAC", "AAG"});
  auto trace = TypeParam::Trace(g1, g2);

  EXPECT_LE(trace.NumRounds(), this->k);
  EXPECT_EQ(trace.BlockDigest(trace.NumRounds() - 1),
            trace.BlockDigest(trace.NumRounds()));
}

// ------------------------------------------------------------
// Correctness of final merge (oracle-based)
// ------------------------------------------------------------

TYPED_TEST(ELMMergeTest, MergeEqualsNaiveUnion) {
  auto g1 = TypeParam::Build(this->k, {"ATG", "TGA"});
  auto g2 = TypeParam::Build(this->k, {"TGA", "GAT"});
  auto merged = TypeParam::Merge(g1, g2);

  auto view = TypeParam::View(merged);
  NaiveDBG expected =
      NaiveDBG::Union(NaiveDBG::FromKmers(this->k, {"ATG", "TGA"}),
                      NaiveDBG::FromKmers(this->k, {"TGA", "GAT"}));

  EXPECT_TRUE(ToNaiveDBG(*view).Equals(expected));
}

TYPED_TEST(ELMMergeTest, MergeIsCommutative) {
  auto g1 = TypeParam::Build(this->k, {"AAA", "AAT"});
  auto g2 = TypeParam::Build(this->k, {"AAT", "ATG"});

  auto m12 = TypeParam::Merge(g1, g2);
  auto m21 = TypeParam::Merge(g2, g1);

  EXPECT_EQ(CanonicalDigest(*TypeParam::View(m12)),
            CanonicalDigest(*TypeParam::View(m21)));
}

TYPED_TEST(ELMMergeTest, MergeIdempotence) {
  auto g = TypeParam::Build(this->k, {"ATG", "TGA"});
  auto merged = TypeParam::Merge(g, g);

  EXPECT_EQ(CanonicalDigest(*TypeParam::View(g)),
            CanonicalDigest(*TypeParam::View(merged)));
}

TYPED_TEST(ELMMergeTest, RandomRefinementMatchesNaiveSuffixClasses) {
  if (!TypeParam::HasRefinementTrace())
    return;

  constexpr size_t trials = 50;
  constexpr size_t kmers_per_graph = 40;

  for (size_t t = 0; t < trials; ++t) {
    auto km1 = GenerateRandomKmers(this->k, kmers_per_graph, 1000 + t);
    auto km2 = GenerateRandomKmers(this->k, kmers_per_graph, 2000 + t);

    auto g1 = TypeParam::Build(this->k, km1);
    auto g2 = TypeParam::Build(this->k, km2);
    auto trace = TypeParam::Trace(g1, g2);

    // Build naive edge list
    NaiveDBG naive = NaiveDBG::Union(NaiveDBG::FromKmers(this->k, km1),
                                     NaiveDBG::FromKmers(this->k, km2));

    std::vector<std::string> edges;
    for (const auto &e : naive.Edges()) {
      edges.push_back(std::get<0>(e) + std::string(1, std::get<1>(e)));
    }

    for (size_t r = 0; r <= trace.NumRounds(); ++r) {
      auto blocks = trace.BlockAssignment(r);

      // Check: edges with same r-suffix must be in same block
      for (size_t i = 0; i < edges.size(); ++i) {
        for (size_t j = 0; j < edges.size(); ++j) {
          auto suf_i = edges[i].substr(this->k - r);
          auto suf_j = edges[j].substr(this->k - r);

          if (suf_i == suf_j) {
            EXPECT_EQ(blocks[i], blocks[j]) << "Mismatch at round " << r;
          }
        }
      }
    }
  }
}

TYPED_TEST(ELMMergeTest, RandomRefinementIsMonotone) {
  if (!TypeParam::HasRefinementTrace())
    return;

  constexpr size_t trials = 100;

  for (size_t t = 0; t < trials; ++t) {
    auto km1 = GenerateRandomKmers(this->k, 30, 3000 + t);
    auto km2 = GenerateRandomKmers(this->k, 30, 4000 + t);

    auto trace = TypeParam::Trace(TypeParam::Build(this->k, km1),
                                  TypeParam::Build(this->k, km2));

    for (size_t r = 1; r < trace.NumRounds(); ++r) {
      EXPECT_GE(trace.BlockCount(r), trace.BlockCount(r - 1));
    }
  }
}

TYPED_TEST(ELMMergeTest, LargeMergeDigestStability) {
  constexpr size_t large = 10000;

  std::vector<std::string> km1, km2;
  GenerateOverlappingKmers(this->k, large, 0.5, 9000, km1, km2);

  auto g1 = TypeParam::Build(this->k, km1);
  auto g2 = TypeParam::Build(this->k, km2);

  auto merged = TypeParam::Merge(g1, g2);
  auto merged2 = TypeParam::Merge(g2, g1);

  EXPECT_EQ(CanonicalDigest(*TypeParam::View(merged)),
            CanonicalDigest(*TypeParam::View(merged2)));
}

TYPED_TEST(ELMMergeTest, LargeMergeAssociativity) {
  constexpr size_t block = 5000;

  auto kmA = GenerateRandomKmers(this->k, block, 1111);
  auto kmB = GenerateRandomKmers(this->k, block, 2222);
  auto kmC = GenerateRandomKmers(this->k, block, 3333);

  auto gA = TypeParam::Build(this->k, kmA);
  auto gB = TypeParam::Build(this->k, kmB);
  auto gC = TypeParam::Build(this->k, kmC);

  auto m1 = TypeParam::Merge(TypeParam::Merge(gA, gB), gC);
  auto m2 = TypeParam::Merge(gA, TypeParam::Merge(gB, gC));

  EXPECT_EQ(CanonicalDigest(*TypeParam::View(m1)),
            CanonicalDigest(*TypeParam::View(m2)));
}

TYPED_TEST(ELMMergeTest, LargeMergeSanityChecks) {
  constexpr size_t large = 20000;

  auto km = GenerateRandomKmers(this->k, large, 424242);

  auto g = TypeParam::Build(this->k, km);
  auto merged = TypeParam::Merge(g, g);

  auto view = TypeParam::View(merged);

  // Idempotence on large data
  EXPECT_EQ(CanonicalDigest(*view), CanonicalDigest(*TypeParam::View(g)));

  // No node explosion
  EXPECT_LE(view->EnumerateNodes().size(), large * this->k);
}
