#include "maki/core/seq/random.hpp"

Dna4Sequence RandomSequence(size_t n, RandomDna4 &g) {
  Dna4Sequence seq;
  seq.reserve(n);

  for (size_t i = 0; i < n; ++i) {
    seq.push_back(g.next());
  }

  return seq;
}

Dna4Sequence RandomSequence(size_t n, size_t seed) {
  RandomDna4 g(seed);
  return RandomSequence(n, g);
}

Dna4Sequence RandomSequence(size_t low, size_t high, size_t seed) {
  // sequence length generator

  std::mt19937 gen(seed);
  std::uniform_int_distribution<size_t> len(low, high);

  // random sequence generator

  RandomDna4 dna(seed);

  // generate sequence

  return RandomSequence(len(gen), dna);
}

std::vector<Dna4Genome> RandomGenomes(size_t n, size_t seed, size_t low,
                                      size_t high) {
  std::vector<Dna4Genome> colln(n);

  // sequence length generator

  std::mt19937 gen(seed);
  std::uniform_int_distribution<size_t> len(low, high);

  // random sequence generator

  RandomDna4 dna(seed);

  // generate sequences

  for (size_t i = 0; i < n; ++i) {
    Dna4Contig &ctg = colln[i].contigs.emplace_back();
    AnnotatedSequence const &rec =
        ctg.sequences[0].emplace_back(RandomSequence(len(gen), dna), i + 1);
    ctg.totalLength += rec.sequence.size();
  }

  return colln;
}

std::vector<Dna4Genome> RandomMutations(Dna4Sequence const &seq, size_t n,
                                        float rateLimit, size_t seed) {
  std::vector<Dna4Genome> colln(n);

  // random generators

  std::mt19937 gen(seed);

  RandomDna4 dna(seed);

  size_t mutLimit =
      std::max((size_t)(rateLimit * (float)seq.size()), (size_t)1u);
  std::uniform_int_distribution<size_t> numMuts(1, mutLimit);
  std::uniform_int_distribution<size_t> mutPos(0, seq.size() - 1);

  // generate sequences

  for (size_t i = 0; i < n; ++i) {
    Dna4Contig &ctg = colln[i].contigs.emplace_back();
    AnnotatedSequence &rec = ctg.sequences[0].emplace_back(seq, i + 1);
    ctg.totalLength += rec.sequence.size();

    size_t z = 0, M = numMuts(gen);

    for (; z < M; ++z)
      rec.sequence[mutPos(gen)] = dna.next();
  }

  return colln;
}
