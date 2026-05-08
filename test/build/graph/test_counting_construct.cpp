
#include "maki/core/utils/tempfile.hpp"
#include "test_common.hpp"

class WeightedGraphTest : public testing::Test {
protected:
  WeightedGraphTest() : bufferBase(tempio::create_temporary_directory()) {}

  ~WeightedGraphTest() { fs::remove_all(bufferBase); }

  inline void verify_graph(const std::vector<std::string> &seqs,
                           std::size_t k) {
    VerifyGraph(seqs, k, bufferBase);
  }

  fs::path bufferBase;
};

TEST_F(WeightedGraphTest, SingleSequenceNoRepeat) {
  std::vector<std::string> seqs = {"ACGTACGT"};
  verify_graph(seqs, 4);
}

TEST_F(WeightedGraphTest, SimpleRepeat) {
  std::vector<std::string> seqs = {"AAAAAA"};
  verify_graph(seqs, 4);
}

TEST_F(WeightedGraphTest, MotifRepeat) {
  std::vector<std::string> seqs = {motif_repeat("ATGC", 100)};
  verify_graph(seqs, 4);
}

TEST_F(WeightedGraphTest, RandomLowRepeat) {
  std::mt19937 rng(42);
  auto seqs = generate_sequences(5, 200, 0.1, rng);
  verify_graph(seqs, 5);
}

TEST_F(WeightedGraphTest, RandomHighRepeat) {
  std::mt19937 rng(42);
  auto seqs = generate_sequences(5, 200, 0.8, rng);
  verify_graph(seqs, 5);
}

TEST_F(WeightedGraphTest, VaryLengths) {
  std::mt19937 rng(42);

  for (size_t len : {50, 100, 500}) {
    auto seqs = generate_sequences(3, len, 0.3, rng);
    verify_graph(seqs, 4);
  }
}

TEST_F(WeightedGraphTest, VaryK) {
  std::mt19937 rng(42);
  auto seqs = generate_sequences(4, 150, 0.2, rng);

  for (size_t k : {4, 5, 7}) {
    verify_graph(seqs, k);
  }
}

TEST_F(WeightedGraphTest, EmptySequence) {
  std::vector<std::string> seqs = {""};
  verify_graph(seqs, 4);
}

TEST_F(WeightedGraphTest, KGreaterThanLength) {
  std::vector<std::string> seqs = {"ACG"};
  verify_graph(seqs, 10);
}

TEST_F(WeightedGraphTest, SingleCharAlphabet) {
  std::vector<std::string> seqs = {"CCCCCCCC"};
  verify_graph(seqs, 4);
}
