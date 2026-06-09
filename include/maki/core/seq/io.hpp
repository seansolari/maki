
#pragma once
#include <string>
#include <string_view>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

namespace parsing {
template <typename O> inline O dna4ToRank(seqan3::dna4 &&nt) {
  return static_cast<O>(seqan3::to_rank(std::forward<seqan3::dna4>(nt)));
}

inline uint8_t dna4ToShort(seqan3::dna4 &&nt) {
  return dna4ToRank<uint8_t>(std::forward<seqan3::dna4>(nt));
}

inline uint64_t dna4ToLong(seqan3::dna4 &&nt) {
  return dna4ToRank<uint64_t>(std::forward<seqan3::dna4>(nt));
}

inline uint64_t dna4ToDna5(uint64_t &&nt) { return (nt + 1) & 0b00000111; }

inline seqan3::dna4 dna5ToDna4(uint64_t &&rnk) {
  return seqan3::assign_rank_to((rnk - 1) & 0b11, seqan3::dna4{});
}
} // namespace parsing

using Dna4Sequence = seqan3::bitpacked_sequence<seqan3::dna4>;
using Dna4SequenceConstIter = Dna4Sequence::const_iterator;

/**
 * writeMerTo
 * ------------
 *
 * Encode base-pairs from input iterator into values pointed
 * at by output iterator. E.g. encodes 10-mer into bytes of
 * the following orientation:
 *
 *      `| 4 3 2 1 | 8 7 6 5 | _ _ 10 9 |`
 *
 */
template <typename OutputPtr, typename O = std::remove_pointer_t<OutputPtr>>
inline OutputPtr writeMerTo(Dna4SequenceConstIter it, OutputPtr out,
                            uint8_t k) {
  constexpr size_t bpPerRecord = sizeof(O) * 4; // base-pairs per output record

  // write complete cells
  for (size_t i = 0; i < k / bpPerRecord; ++i) {
    O &val = *out++;
    for (size_t j = 0; j < bpPerRecord; ++j) {
      val |= parsing::dna4ToRank<O>(*it++) << (2 * j);
    }
  }

  // write partial cell
  O &val = *out;
  for (size_t j = 0; j < k % bpPerRecord; ++j) {
    val |= parsing::dna4ToRank<O>(*it++) << (2 * j);
  }
  return out;
}

template <typename OutputPtr, typename O = std::remove_pointer_t<OutputPtr>>
  requires std::random_access_iterator<OutputPtr>
inline void reverseWriteMerTo(Dna4SequenceConstIter it, OutputPtr out,
                              uint8_t k, uint8_t n) {
  constexpr size_t bpPerRecord = sizeof(O) * 4; // base-pairs per output record
  for (long i = k - 1, end = k - 1 - n; i > end; --i)
    *(out + (i / bpPerRecord)) |= parsing::dna4ToRank<O>(*it--)
                                  << (2 * (i % bpPerRecord));
}

std::string toString(const Dna4Sequence &);
std::string toString(std::vector<Dna4Sequence> const &, const char *);

enum InputFileType {
  FastaFileType,
  Gff3FileType,
  FastQFileType,
  UnknownFileType
};

constexpr std::string_view COMPRESSED_EXTENSIONS[2] = {".tar.gz", ".gz"};
constexpr std::string_view FNA_EXTENSIONS[3] = {".fna", ".fa", ".fasta"};
constexpr std::string_view GFF_EXTENSIONS[2] = {".gff3", ".gff"};
constexpr std::string_view FQ_EXTENSIONS[2] = {".fq", ".fastq"};

std::string_view removeCompressedExtensions(std::string_view filePath);

InputFileType detectFileType(std::string_view inputFile);
std::string_view extractSequenceName(std::string_view path);

struct GenomeManifest {
  std::vector<std::string> files;
  InputFileType type;
};

GenomeManifest readFilePaths(const char *manifest_file, std::size_t col,
                             const char sep, InputFileType filter);

GenomeManifest readManifest(const char *manifest_file, InputFileType filter);

struct DataFilePair {
  std::string forwardFile, reverseFile;

  DataFilePair() = default;
  DataFilePair(std::string &&fwd, std::string &&rev)
      : forwardFile(std::move(fwd)), reverseFile(std::move(rev)) {}

  friend std::ostream &operator<<(std::ostream &, const DataFilePair &obj);
  inline bool operator==(DataFilePair const &) const = default;
  bool gzCompressed() const;
};
