#pragma once

#include "maki/build/io/gff.hpp"
#include "maki/core/seq/concepts.hpp"

// ---------------------------------------------------------------------------
// Containers
// ---------------------------------------------------------------------------

struct FastaFragment {
  FastaFragment(Dna4SequenceConstIter it, std::size_t begin, std::size_t end,
                uint64_t id)
      : _it(it), _begin(begin), _end(end), _id(id) {}

  inline std::size_t size() const noexcept { return _end - _begin; }
  inline uint64_t id() const noexcept { return _id; }

  inline void setEndToTerminal() { _terminal = true; }
  inline bool endIsTerminal() const noexcept { return _terminal; }

  std::size_t numTerminals(std::size_t k) const;
  std::size_t numKmers(std::size_t k) const;

  using view_type =
      SequenceFragment<std::ranges::subrange<Dna4SequenceConstIter>>;

  view_type terminals() const;
  view_type fragments(std::size_t k) const;

protected:
  Dna4SequenceConstIter _it;
  std::size_t _begin;
  std::size_t _end;
  uint64_t _id;
  bool _terminal = false;
};

static_assert(sequence_like<FastaFragment>);

class ChunkedDna4Genome {
public:
  ChunkedDna4Genome(Dna4Genome &&genome, std::size_t ntChunksize,
                    std::size_t overlap);

  static std::size_t calculateChunksize(const Dna4Genome &genome, std::size_t granularity);  
  constexpr std::size_t length() const noexcept { return _genome.length(); }
  constexpr std::size_t size() const noexcept { return _chunks.size(); }
  constexpr const FastaFragment &operator[](std::size_t i) const noexcept { return _chunks[i]; }

private:
  Dna4Genome _genome;
  std::vector<FastaFragment> _chunks;

protected:
  void _chunk(std::size_t granularity, std::size_t overlap);
};

static_assert(container_span<ChunkedDna4Genome>);

// ---------------------------------------------------------------------------
// API
// ---------------------------------------------------------------------------

void parseFastaStream(Dna4Genome &genome, std::istream &fastaStream,
                      Colours &colours, std::size_t minContigSize);

Dna4Genome parseGenome(const std::string &fastaFile, Colours &colours,
                       std::size_t k);

ChunkedDna4Genome chunkFNA(Dna4Genome &&genome, std::size_t ntChunksize,
                           std::size_t k);
