
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "test_common.hpp"
#include <gtest/gtest.h>

TEST(LongSuffixGateSerial, InitiallyEmpty) {
  LongSuffixGate gate(/*length=*/4);
  EXPECT_EQ(gate.count(), 0u);
}

TEST(LongSuffixGateSerial, SingleSuffixAllEdgesOnce) {
  LongSuffixGate gate(4);
  auto sfx = makeLongSuffixFromDNA("AC");

  for (uint8_t edge : AllEdges) {
    auto result = gate.trySet(sfx, edge);
    EXPECT_TRUE(result.has_value());
  }

  EXPECT_EQ(gate.count(), AllEdges.size());
}

TEST(LongSuffixGateSerial, DuplicateInsertRejected) {
  LongSuffixGate gate(4);
  auto sfx = makeLongSuffixFromDNA("GTA");

  auto first = gate.trySet(sfx, A);
  EXPECT_TRUE(first.has_value());

  auto second = gate.trySet(sfx, A);
  EXPECT_FALSE(second.has_value());

  EXPECT_EQ(gate.count(), 1u);
}

TEST(LongSuffixGateSerial, SameEdgeDifferentSuffixes) {
  LongSuffixGate gate(4);

  auto s1 = makeLongSuffixFromDNA("A");
  auto s2 = makeLongSuffixFromDNA("AA");
  auto s3 = makeLongSuffixFromDNA("AAA");
  auto s4 = makeLongSuffixFromDNA("AAAA");

  EXPECT_TRUE(gate.trySet(s1, T).has_value());
  EXPECT_TRUE(gate.trySet(s2, T).has_value());
  EXPECT_TRUE(gate.trySet(s3, T).has_value());
  EXPECT_TRUE(gate.trySet(s4, T).has_value());

  EXPECT_EQ(gate.count(), 4u);
}

TEST(LongSuffixGateSerial, TerminalEdgeDistinctFromDNAEdges) {
  LongSuffixGate gate(4);
  auto sfx = makeLongSuffixFromDNA("CG");

  EXPECT_TRUE(gate.trySet(sfx, TerminalEdge).has_value());
  EXPECT_TRUE(gate.trySet(sfx, A).has_value());

  EXPECT_EQ(gate.count(), 2u);
}

TEST(LongSuffixGateSerial, AllSuffixLengthsUpToFour) {
  LongSuffixGate gate(4);

  std::vector<std::string> suffixes{"", "A", "AC", "ACG", "ACGT"};

  size_t expected = 0;

  for (auto &dna : suffixes) {
    auto sfx = makeLongSuffixFromDNA(dna);
    for (uint8_t edge : AllEdges) {
      if (gate.trySet(sfx, edge))
        ++expected;
    }
  }

  EXPECT_EQ(gate.count(), expected);
}

TEST(LongSuffixGateParallel, SameSuffixSameEdgeContention) {
  LongSuffixGate gate(4);
  auto sfx = makeLongSuffixFromDNA("ACG");

  constexpr int threads = 16;
  std::vector<std::thread> workers;
  std::atomic<int> successes{0};

  for (int i = 0; i < threads; ++i) {
    workers.emplace_back([&] {
      if (gate.trySet(sfx, G).has_value())
        ++successes;
    });
  }

  for (auto &t : workers)
    t.join();

  EXPECT_EQ(successes.load(), 1);
  EXPECT_EQ(gate.count(), 1u);
}

TEST(LongSuffixGateParallel, SameSuffixDifferentEdgesParallel) {
  LongSuffixGate gate(4);
  auto sfx = makeLongSuffixFromDNA("TT");

  std::vector<std::thread> workers;

  for (uint8_t edge : AllEdges) {
    workers.emplace_back([&, edge] { gate.trySet(sfx, edge); });
  }

  for (auto &t : workers)
    t.join();

  EXPECT_EQ(gate.count(), AllEdges.size());
}

TEST(LongSuffixGateParallel, MixedSuffixesAndEdges) {
  LongSuffixGate gate(4);
  std::vector<std::string> suffixes{"A", "C", "G", "T", "AA", "AC", "GT", "TT"};

  std::vector<std::thread> workers;

  for (auto &dna : suffixes) {
    for (uint8_t edge : AllEdges) {
      workers.emplace_back([&, dna, edge] {
        auto sfx = makeLongSuffixFromDNA(dna);
        gate.trySet(sfx, edge);
      });
    }
  }

  for (auto &t : workers)
    t.join();

  EXPECT_EQ(gate.count(), suffixes.size() * AllEdges.size());
}

TEST(LongSuffixGateParallel, RepeatedParallelReentry) {
  for (int round = 0; round < 100; ++round) {
    LongSuffixGate gate(4);
    auto sfx = makeLongSuffixFromDNA("ACGT");

    std::vector<std::thread> threads;
    for (int i = 0; i < 8; ++i) {
      threads.emplace_back([&] {
        for (uint8_t e : AllEdges)
          gate.trySet(sfx, e);
      });
    }

    for (auto &t : threads)
      t.join();

    EXPECT_EQ(gate.count(), AllEdges.size());
  }
}

TEST(LongSuffixGateSerial, ExhaustiveAllSuffixEdgeCombinations) {
  constexpr size_t MAX_LEN = 4; // changeable parameter

  LongSuffixGate gate(MAX_LEN);

  // Enumerate all DNA suffixes of length <= MAX_LEN
  std::vector<std::string> suffixes;
  enumerateDNA(MAX_LEN, suffixes);

  size_t expectedInserted = 0;

  for (const auto &dna : suffixes) {
    auto sfx = makeLongSuffixFromDNA(dna);

    for (uint8_t edge : AllEdges) {
      auto result = gate.trySet(sfx, edge);

      // Every combination must succeed exactly once
      ASSERT_TRUE(result.has_value()) << "Failed insertion: suffix=\"" << dna
                                      << "\" edge=" << static_cast<int>(edge);

      ++expectedInserted;
    }
  }

  // Final count must match exact combinatorial expectation
  EXPECT_EQ(expectedInserted, ShortSuffix::numSuffixes(MAX_LEN) * 5);
  EXPECT_EQ(gate.count(), expectedInserted);
}

TEST(LongSuffixGateParallel, ExhaustiveAllSuffixEdgeCombinations) {
  constexpr size_t MAX_LEN = 4; // changeable parameter

  LongSuffixGate gate(MAX_LEN);

  // Enumerate all DNA suffixes of length <= MAX_LEN
  std::vector<std::string> suffixes;
  enumerateDNA(MAX_LEN, suffixes);

  // Total expected combinations:
  // (sum_{k=0..MAX_LEN} 4^k) * 5 edges
  const size_t expectedTotal = suffixes.size() * AllEdges.size();
  EXPECT_EQ(expectedTotal, ShortSuffix::numSuffixes(MAX_LEN) * 5);

  // We parallelize over suffixes; each task handles all edges
  std::vector<std::thread> workers;
  workers.reserve(suffixes.size());

  for (const auto &dna : suffixes) {
    workers.emplace_back([&, dna] {
      auto sfx = makeLongSuffixFromDNA(dna);

      for (uint8_t edge : AllEdges) {
        gate.trySet(sfx, edge);
      }
    });
  }

  for (auto &t : workers)
    t.join();

  EXPECT_EQ(gate.count(), expectedTotal);
}
