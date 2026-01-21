#pragma once
#include <array>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <random>
#include <ranges>
#include <unordered_set>
#include <span>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <glog/logging.h>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/all.hpp>
#include <zstr.hpp>
#include "maki/colour/encoder.hpp"

namespace fs = std::filesystem;

using namespace seqan3::literals;

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

std::string toString(const Dna4Sequence &);
std::string toString(std::vector<Dna4Sequence> const &, const char *);

// ---------------------------------------------------------------------------
// Parsing GFF files
// ---------------------------------------------------------------------------

struct gffToken
{
  uint64_t featureId;
  size_t begin, end;

  bool operator==(const gffToken &) const = default;
  bool operator!=(const gffToken &) const = default;

  size_t size() const noexcept { return end - begin; }
};

struct gffRecord
{
  std::string accn;        // sequence accession
  int64_t begin = 0;       // start (field 4), GFF: 1-based inclusive
  int64_t end = 0;         // end   (field 5), GFF: 1-based inclusive
  char strand = '.';       // '+', '-', '.', '?' (field 7)
  uint64_t featureId = 0;  // (internal) feature ID

  friend bool operator==(const gffRecord &, const gffRecord &);
  friend bool operator!=(const gffRecord &, const gffRecord &);

  inline gffToken asToken() const noexcept { return {featureId, begin, end}; }
  operator gffToken() const noexcept { return asToken(); }
};

// Parse GFF line, extracting seed name if it has an ID attribute
std::pair<gffRecord,std::optional<std::string>> parseGffLine(std::string_view line);

using RecordVector = std::vector<gffRecord>;
using RecordVectorIter = RecordVector::iterator;
using RecordVectorConstIter = RecordVector::const_iterator;

using AnnotationTokenVector = std::vector<gffToken>;
using AnnotationTokenVectorIter = AnnotationTokenVector::iterator;
using AnnotationTokenVectorConstIter = AnnotationTokenVector::const_iterator;

struct AnnotRange
{
  RecordVectorIter begin, end;
  void sortByStart();
  AnnotRange getStrand(char c) const;
  RecordVectorIter encloseStart(size_t z) const;
  size_t size() const noexcept { return end - begin; }
};

struct GenomeAnnotationList
{
  RecordVector records;
  size_t size() const noexcept { return records.size(); }
  void sortContigStrand();
  void sort();
  AnnotRange findContig(const std::string &id);

protected:
  static void ensureDisjointRegion(RecordVectorConstIter it, RecordVectorConstIter end);
};

struct ContigIdLess
{
  bool operator()(const gffRecord &lhs, const gffRecord &rhs) const { return lhs.accn < rhs.accn; }
  bool operator()(const gffRecord &lhs, const std::string &rhs) const { return lhs.accn < rhs; }
  bool operator()(const std::string &lhs, const gffRecord &rhs) const { return lhs < rhs.accn; }
};

struct StrandLess
{
  bool operator()(const gffRecord &lhs, const gffRecord &rhs) const { return lhs.strand < rhs.strand; }
  bool operator()(const gffRecord &lhs, char rhs) const { return lhs.strand < rhs; }
  bool operator()(char lhs, const gffRecord &rhs) const { return lhs < rhs.strand; }
};

struct BeginPosLess
{
  bool operator()(const gffRecord &lhs, const gffRecord &rhs) const { return lhs.begin < rhs.begin; }
  bool operator()(const gffRecord &lhs, size_t rhs) const { return lhs.begin < rhs; }
  bool operator()(size_t lhs, const gffRecord &rhs) const { return lhs < rhs.begin; }
};

struct ContigStrandLess
{
  bool operator()(const gffRecord &lhs, const gffRecord &rhs) const;
};

struct AnnotLess
{
  bool operator()(const gffRecord &lhs, const gffRecord &rhs) const;
};

using Dna4SequenceConstIter = Dna4Sequence::const_iterator;
using AnnotationTokenList = std::vector<gffToken>;

// For an edge to be annotated with a colour, that annotation must have started at least `k` positions ago.
class AnnotatedSequenceEdgeIterator
{
public:
  AnnotatedSequenceEdgeIterator(Dna4SequenceConstIter seqIter, size_t startPos, uint8_t windowSize, uint64_t nullFeatureId, AnnotationTokenVectorConstIter annotIter, AnnotationTokenVectorConstIter annotEnd);
  AnnotatedSequenceEdgeIterator(Dna4SequenceConstIter seqIter, uint8_t windowSize, uint64_t nullFeatureId, AnnotationTokenVectorConstIter annotIter, AnnotationTokenVectorConstIter annotEnd) : AnnotatedSequenceEdgeIterator(seqIter, 0, windowSize, nullFeatureId, annotIter, annotEnd) {}
  AnnotatedSequenceEdgeIterator(Dna4SequenceConstIter seqIter) : inputSeqIter(seqIter) {}

protected:
  Dna4SequenceConstIter inputSeqIter;
  size_t inputSeqPos;
  uint8_t k;
  uint64_t nullFeatureId;
  AnnotationTokenVectorConstIter inputAnnotIter, inputAnnotEnd;
  AnnotationTokenList currentAnnots;
  size_t numCurrentAnnots;

public:
  bool operator==(const AnnotatedSequenceEdgeIterator &other) const noexcept { return inputSeqIter == other.inputSeqIter; }
  bool operator!=(const AnnotatedSequenceEdgeIterator &other) const noexcept { return inputSeqIter != other.inputSeqIter; }
  
  void operator++();
  std::iter_value_t<Dna4SequenceConstIter> operator*() const { return *inputSeqIter; }
  
  uint64_t getNullFeatureId() const { return nullFeatureId; }
  bool positionIsAnnotated() const { return numCurrentAnnots != 0; }
  size_t getNumCurrentAnnots() const { return std::max(numCurrentAnnots, (size_t)1); }
  const AnnotationTokenList &viewCurrentAnnots() const { return currentAnnots; }
  std::string currentAnnotsToString() const;
};

struct AnnotatedSequenceSegment
{
  Dna4SequenceConstIter seqBegin, seqEnd;
  uint64_t featureId;
  bool endIsTerminal;
  size_t size() const noexcept { return seqEnd - seqBegin; }
};

class AnnotatedSequenceSegmentSentinel
{
};

class AnnotatedSequenceSegmentIterator
{
public:
    AnnotatedSequenceSegmentIterator(Dna4SequenceConstIter seqIter, size_t seqLength, uint8_t k_, uint64_t nullFeatureId, AnnotationTokenVectorConstIter annotIter, AnnotationTokenVectorConstIter annotEnd)
        : inputSeqIter(seqIter)
        , inputEdgePos(k_)
        , inputSeqLength(seqLength)
        , k(k_)
        , nullFeatureId(nullFeatureId)
        , inputAnnotIter(annotIter)
        , inputAnnotEnd(annotEnd)
    {
        operator++();
    }

protected:
  Dna4SequenceConstIter inputSeqIter;
  size_t inputEdgePos, inputSeqLength;
  uint8_t k;
  uint64_t nullFeatureId;
  AnnotationTokenVectorConstIter inputAnnotIter, inputAnnotEnd;
  size_t currentSegmentBegin, currentSegmentEnd;
  uint64_t currentSegmentFeatureId;

public:
  inline bool operator==([[maybe_unused]] const AnnotatedSequenceSegmentSentinel &) const noexcept
  {
    return (currentSegmentBegin == inputSeqLength) && (currentSegmentEnd == inputSeqLength);
  }

  inline bool operator!=([[maybe_unused]] const AnnotatedSequenceSegmentSentinel &) const noexcept
  {
    return (currentSegmentBegin != inputSeqLength) || (currentSegmentEnd != inputSeqLength);
  }

  AnnotatedSequenceSegment operator*() const;

  // does not count terminal edge
  inline size_t segmentNumInternalEdges() const noexcept { return currentSegmentEnd - currentSegmentBegin - k; }

  void operator++();
};

struct AnnotatedSequence
{
  Dna4Sequence sequence;
  AnnotationTokenVector annotations;
  uint64_t nullFeatureId;

  AnnotatedSequence(Dna4Sequence const &seq) : sequence(seq), annotations(), nullFeatureId(0) {}
  AnnotatedSequence(Dna4Sequence const &seq, uint64_t id) : sequence(seq), annotations(), nullFeatureId(id) {}
  AnnotatedSequence(Dna4Sequence &&seq) : sequence(std::move(seq)), annotations(), nullFeatureId(0) {}
  AnnotatedSequence(Dna4Sequence &&seq, uint64_t id) : sequence(std::move(seq)), annotations(), nullFeatureId(id) {}
  AnnotatedSequence(uint64_t id, size_t length) : sequence(), annotations(), nullFeatureId(id){ sequence.reserve(length); }
  
  inline AnnotatedSequenceEdgeIterator colouredSegmentsEdgeBegin(uint8_t windowSize) const { return AnnotatedSequenceEdgeIterator(sequence.cbegin(), windowSize, annotations.cbegin(), annotations.cend()); }
  inline AnnotatedSequenceEdgeIterator colouredSegmentsEdgeBegin(size_t startPos, uint8_t windowSize) const { return AnnotatedSequenceEdgeIterator(sequence.cbegin(), static_cast<int32_t>(startPos), windowSize, annotations.cbegin(), annotations.cend()); }
  inline AnnotatedSequenceEdgeIterator colouredSegmentsEdgeEnd() const { return AnnotatedSequenceEdgeIterator(sequence.cend()); }
  inline AnnotatedSequenceSegmentIterator colouredSegmentsBegin(uint8_t windowSize) const { return AnnotatedSequenceSegmentIterator(sequence.cbegin(), sequence.size(), windowSize, annotations.cbegin(), annotations.cend()); }
  inline AnnotatedSequenceSegmentSentinel colouredSegmentsEnd() const { return AnnotatedSequenceSegmentSentinel(); }
  
  size_t numKmers(uint8_t k) const;
};

using SequenceVector = std::vector<AnnotatedSequence>;

struct Dna4Contig
{
  std::pair</*forward*/SequenceVector,/*reverse*/SequenceVector> sequences;
  std::string accn;
  size_t minFragmentSize = 0;
  size_t totalLength = 0;
  size_t numAnnotations = 0;

  Dna4Contig() = default;
  Dna4Contig(std::string accn_, size_t minContigSize_)
    : accn(accn_), minFragmentSize(minContigSize_) {}

  void insert(const seqan3::dna5_vector &, AnnotRange annots, size_t minAnnotSize);
  size_t numKmers(uint8_t k) const;
  inline size_t numFragments() const noexcept { return sequences.first.size() + sequences.second.size(); }

protected:
  std::vector<size_t> _findNs(const seqan3::dna5_vector &);
  void _insertOrientation(SequenceVector &, R &&seqView, uint64_t seqFeatureId, const std::vector<size_t> &nPos, AnnotRange annots, size_t minAnnotSize);
  void _insertFragment(SequenceVector &vec, Iter seqIter, size_t startPos, size_t endPos, uint64_t seqFeatureId, AnnotRange &annots, AnnotationTokenList &currentAnnots, size_t minAnnotSize);
};

struct Dna4Genome
{
  std::vector<Dna4Contig> contigs;
  Dna4Genome() = default;
  size_t length() const;
  // includes concluding terminal edge, but not starting terminals
  size_t numKmers(uint8_t k) const;
  inline size_t numContigs() const noexcept { return contigs.size(); }
  size_t numFragments() const;
  // number of starting terminals, does not include trailing terminal edge
  size_t numTerminals(uint8_t k) const;
  size_t medianContigSize() const;
  // number of bytes used
  size_t rss() const;
  void reset_memory();
};

void parseAnnotationStream(GenomeAnnotationList &annots, Colours &colours, zstr::ifstream &gffStream);
void parseFastaStream(Dna4Genome &genome, std::istream &fastaStream, size_t minContigSize);
Dna4Genome parseGenome(fs::path fastaFile, size_t minContigSize);
Dna4Genome parseAnnotatedGenome(fs::path gffFile, Colours &colours, size_t minContigSize);
std::vector<fs::path> readFilePaths(const char *manifest_file);

class Dna4GenomeVector
{
public:
  using BaseGenomeVector = std::vector<Dna4Genome>;
  using BaseGenomeVectorIter = BaseGenomeVector::iterator;
  using BaseGenomeVectorConstIter = BaseGenomeVector::const_iterator;

protected:
  BaseGenomeVector genomes;

public:
  Dna4GenomeVector(size_t n) : genomes(n) {}
  operator BaseGenomeVector const &() const noexcept { return genomes; }
  BaseGenomeVector const &buffer() const noexcept { return genomes; }
  // sends current data to temporary object, whose deconstructor will then be called
  void reset_memory() { BaseGenomeVector().swap(genomes); }
  size_t numGenomes() const { return genomes.size(); }
  size_t totalLength() const;
  size_t numKmers(uint8_t windowSize) const;
  Dna4Genome &operator[](size_t i) noexcept { return genomes[i]; }
  const Dna4Genome &operator[](size_t i) const noexcept { return genomes[i]; }
  Dna4Genome *data() noexcept { return genomes.data(); }
  const Dna4Genome *data() const noexcept { return genomes.data(); }
  BaseGenomeVectorIter begin() { return genomes.begin(); }
  BaseGenomeVectorIter end() { return genomes.end(); }
  BaseGenomeVectorConstIter begin() const { return genomes.cbegin(); }
  BaseGenomeVectorConstIter end() const { return genomes.cend(); }
  BaseGenomeVectorConstIter cbegin() const { return genomes.cbegin(); }
  BaseGenomeVectorConstIter cend() const { return genomes.cend(); }
};

Dna4GenomeVector loadGenomes(const std::vector<fs::path> &fastaFiles, Colours &colours, size_t minContigSize);
Dna4GenomeVector loadGenomes(const std::vector<fs::path> &fastaFiles, size_t minContigSize);

// ---------------------------------------------------------------------------
// Parsing large FNA files (filters)
// ---------------------------------------------------------------------------

class ChunkedDna4Genome
{
public:
  typedef std::vector<SequenceRange>::const_iterator iterator;
  ChunkedDna4Genome(Dna4Genome &&genome, size_t granularity, size_t overlap);
  auto begin() const { return _chunks.begin(); }
  auto end() const { return _chunks.end(); }
  operator Dna4Genome const &() const noexcept { return _genome; }
  operator std::vector<SequenceRange> &() noexcept { return _chunks; }
  operator const std::vector<SequenceRange> &() const noexcept { return _chunks; }
  SequenceRange const &operator[](size_t i) const { return _chunks[i]; }
  size_t chunks() const noexcept { return _chunks.size(); }
  size_t numEdges() const noexcept { return _genome.length() + _genome.numFragments(); }
  size_t sequenceRSS() const { return _genome.rss(); }
  size_t chunksRSS() const { return sizeof(SequenceRange) * _chunks.capacity(); }
  void reset_memory();

protected:
  Dna4Genome _genome;
  std::vector<AnnotatedSequenceSegment> _chunks;
};

// ---------------------------------------------------------------------------
// API
// ---------------------------------------------------------------------------

class SequenceContainer
{
  virtual std::size_t numTerminals(std::size_t k_) const =0;

  virtual Dna4SequenceConstIter begin() const =0;
  virtual Dna4SequenceConstIter end() const =0;
};

std::vector<SequenceContainer*> constructView(const Dna4GenomeVector&, const ChunkedDna4Genome&);
