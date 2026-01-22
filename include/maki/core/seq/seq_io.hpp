
#pragma once
#include <string>
#include <string_view>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>

namespace parsing
{
  template <typename O>
  inline O dna4ToRank(seqan3::dna4 &&nt)
  {
    return static_cast<O>(seqan3::to_rank(std::forward<seqan3::dna4>(nt)));
  }

  inline uint8_t dna4ToShort(seqan3::dna4 &&nt)
  {
    return dna4ToRank<uint8_t>(std::forward<seqan3::dna4>(nt));
  }

  inline uint64_t dna4ToLong(seqan3::dna4 &&nt)
  {
    return dna4ToRank<uint64_t>(std::forward<seqan3::dna4>(nt));
  }

  inline uint64_t dna4ToDna5(uint64_t &&nt)
  {
    return (nt + 1) & 0b00000111;
  }

  inline seqan3::dna4 dna5ToDna4(uint64_t &&rnk)
  {
    return seqan3::assign_rank_to((rnk - 1) & 0b11, seqan3::dna4{});
  }

} // namespace parsing

using Dna4Sequence = seqan3::bitpacked_sequence<seqan3::dna4>;
using Dna4SequenceConstIter = Dna4Sequence::const_iterator;

std::string toString(const Dna4Sequence &);
std::string toString(std::vector<Dna4Sequence> const &, const char *);

enum InputFileType
{
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
std::vector<std::string> readFilePaths(const char *manifest_file, InputFileType filter);

std::string_view extractSequenceName(std::string_view path);
