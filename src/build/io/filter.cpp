#include "maki/build/io/filter.hpp"
#include "maki/core/utils/logging.hpp"

#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/io/sequence_file/all.hpp>

// ---------------------------------------------------------------------------
// Containers
// ---------------------------------------------------------------------------

std::size_t FastaFragment::numTerminals(std::size_t k) const {
  if (_begin == 0)
    return k;
  else
    return 0;
}

std::size_t FastaFragment::numKmers(std::size_t k) const {
  return size() - k + 1;
}

FastaFragment::view_type FastaFragment::terminals() const {
  if (_begin == 0)
    return SequenceFragment(std::ranges::subrange(_it + _begin, _it + _end),
                            _id, _terminal);
  else
    return SequenceFragment(std::ranges::subrange<Dna4SequenceConstIter>(), _id,
                            _terminal);
}

FastaFragment::view_type
FastaFragment::fragments([[maybe_unused]] std::size_t k) const {
  return SequenceFragment(std::ranges::subrange(_it + _begin, _it + _end), _id,
                          _terminal);
}

ChunkedDna4Genome::ChunkedDna4Genome(Dna4Genome &&genome,
                                     std::size_t ntChunksize,
                                     std::size_t overlap)
    : _genome(std::move(genome)), _chunks() {
  for (auto &&fragment : _genome.raw_contigs()) {
    if (fragment.size() < overlap)
      LOG_WARN() << "fragment " << fragment.id() << " size " << fragment.size()
                 << " less than required overlap " << overlap
                 << ", disregarding.";

    else {
      auto ptr = fragment.data().base().cbegin();

      std::size_t i = overlap;
      std::size_t j;

      while (i < fragment.size()) {
        j = std::min(fragment.size(), i + ntChunksize);
        _chunks.emplace_back(ptr, i - overlap, j, fragment.id());
        i = j;
      }

      _chunks.back().setEndToTerminal();
    }
  }
}

std::size_t ChunkedDna4Genome::calculateChunksize(const Dna4Genome &genome,
                                                  std::size_t granularity) {
  return std::max(genome.medianContigSize() / granularity, std::size_t{10'000});
}

// ---------------------------------------------------------------------------
// API
// ---------------------------------------------------------------------------

/**
 * Parse FASTA stream, assigning colours by contig accession.
 */
void parseFastaStream(Dna4Genome &genome, std::istream &fastaStream,
                      Colours &colours, std::size_t minContigSize) {
  seqan3::sequence_file_input seqInput(fastaStream, seqan3::format_fasta{});

  for (auto &&rec : seqInput) {
    Dna4Contig &out = genome.contigs.emplace_back(
        rec.id().substr(0, rec.id().find(' ')), minContigSize);

    uint64_t seqId = colours.getOrAssign(std::string(out.accn));

    // insert forward and reverse sequences

    SequenceVector &fwdSeqs = out.sequences[0], &revSeqs = out.sequences[1];

    for (auto &&subSeq : rec.sequence() | std::views::split('N'_dna5) |
                             parsing::filterByLength(minContigSize)) {
      out.totalLength += 2 * subSeq.size();

      // insert forward

      auto &fwd = fwdSeqs.emplace_back(seqId, subSeq.size());

      for (auto &&c : subSeq | parsing::transformToDna4) {
        fwd.sequence.push_back(c);
      }

      // insert reverse

      auto &rev = revSeqs.emplace_back(seqId, subSeq.size());

      for (auto &&c : subSeq | std::views::reverse | seqan3::views::complement |
                          parsing::transformToDna4) {
        rev.sequence.push_back(c);
      }
    }
  }
}

ChunkedDna4Genome parseFilterFNA(const std::string &fastaFile, Colours &colours,
                                 std::size_t granularity, std::size_t k) {
  // base input file stream that reads bytes
  zstr::ifstream zis(fastaFile);
  std::string fastaData(std::istreambuf_iterator<char>(zis), {});

  // process stream data
  Dna4Genome genome{};

  std::istringstream fs(fastaData);
  parseFastaStream(genome, fs, colours, k);

  // Chunk data
  std::size_t chunksize =
      ChunkedDna4Genome::calculateChunksize(genome, granularity);
  return ChunkedDna4Genome(std::move(genome), chunksize, k);
}
