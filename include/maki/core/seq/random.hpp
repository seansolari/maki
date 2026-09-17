
#include "maki/build/io/gff.hpp"
#include "maki/core/seq/io.hpp"

using namespace seqan3::literals;

struct RandomDna4 {
protected:
  std::array<seqan3::dna4, 4> _ALPHA = {'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};
  std::mt19937 gen;
  std::uniform_int_distribution<> distr;

public:
  template <class SeedSeq>
  RandomDna4(SeedSeq &&seed) : gen(std::forward<SeedSeq>(seed)), distr(0, 3) {}
  seqan3::dna4 next() { return _ALPHA[distr(gen)]; }
};

Dna4Sequence RandomSequence(std::size_t n, RandomDna4 &g);
Dna4Sequence RandomSequence(std::size_t n, std::size_t seed);
Dna4Sequence RandomSequence(std::size_t low, std::size_t high,
                            std::size_t seed);
std::vector<Dna4Genome> RandomGenomes(std::size_t n, std::size_t seed,
                                      std::size_t low, std::size_t high);
std::vector<Dna4Genome> RandomMutations(Dna4Sequence const &seq, std::size_t n,
                                        float rateLimit, std::size_t seed);
