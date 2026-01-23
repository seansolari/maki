#pragma once
#include <cassert>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <optional>
#include <string>
#include <vector>

#include <indicators/progress_bar.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <zstr.hpp>

#include "maki/build/graph/build_colours.hpp"
#include "maki/core/seq/seq_concepts.hpp"
#include "maki/core/seq/seq_io.hpp"

using namespace seqan3::literals;

// ---------------------------------------------------------------------------
// Parsing GFF files
// ---------------------------------------------------------------------------

struct gffToken {
  uint64_t featureId;
  int64_t begin;
  int64_t end;

  bool operator==(const gffToken &) const;
  bool operator!=(const gffToken &) const;

  std::size_t size() const noexcept { return end - begin; }
};

std::ostream &operator<<(std::ostream &, const gffToken &);

struct gffRecord {
  std::string accn;       // sequence accession
  int64_t begin = 0;      // start (field 4), GFF: 1-based inclusive
  int64_t end = 0;        // end   (field 5), GFF: 1-based inclusive
  char strand = '.';      // '+', '-', '.', '?' (field 7)
  uint64_t featureId = 0; // (internal) feature ID

  friend bool operator==(const gffRecord &, const gffRecord &);
  friend bool operator!=(const gffRecord &, const gffRecord &);

  inline gffToken asToken() const noexcept { return {featureId, begin, end}; }
  operator gffToken() const noexcept { return asToken(); }
};

// Parse GFF line, extracting seed name if it has an ID attribute
std::pair<gffRecord, std::optional<std::string>>
parseGffLine(std::string_view line);

using RecordVector = std::vector<gffRecord>;
using RecordVectorIter = RecordVector::iterator;
using RecordVectorConstIter = RecordVector::const_iterator;

using AnnotationTokenVector = std::vector<gffToken>;
using AnnotationTokenVectorIter = AnnotationTokenVector::iterator;
using AnnotationTokenVectorConstIter = AnnotationTokenVector::const_iterator;

struct AnnotRange;

struct GenomeAnnotationList {
  RecordVector records;

protected:
  struct AnnotLess {
    bool operator()(const gffRecord &lhs, const gffRecord &rhs) const;
  };
  struct ContigStrandLess {
    bool operator()(const gffRecord &lhs, const gffRecord &rhs) const;
  };
  struct ContigIdLess {
    bool operator()(const gffRecord &lhs, const gffRecord &rhs) const {
      return lhs.accn < rhs.accn;
    }
    bool operator()(const gffRecord &lhs, const std::string &rhs) const {
      return lhs.accn < rhs;
    }
    bool operator()(const std::string &lhs, const gffRecord &rhs) const {
      return lhs < rhs.accn;
    }
  };

public:
  void sort();
  void sortContigStrand();
  AnnotRange findContig(const std::string &id);
};

/**
 * Represents a range of annotations, usually from a specific contig, sorted by
 * strand.
 */
struct AnnotRange {
  RecordVectorIter begin;
  RecordVectorIter end;

  std::size_t size() const noexcept { return end - begin; }

protected:
  struct StrandLess {
    bool operator()(const gffRecord &lhs, const gffRecord &rhs) const {
      return lhs.strand < rhs.strand;
    }
    bool operator()(const gffRecord &lhs, char rhs) const {
      return lhs.strand < rhs;
    }
    bool operator()(char lhs, const gffRecord &rhs) const {
      return lhs < rhs.strand;
    }
  };

  struct BeginPosLess {
    bool operator()(const gffRecord &lhs, const gffRecord &rhs) const {
      return lhs.begin < rhs.begin;
    }
    bool operator()(const gffRecord &lhs, int64_t rhs) const {
      return lhs.begin < rhs;
    }
    bool operator()(int64_t lhs, const gffRecord &rhs) const {
      return lhs < rhs.begin;
    }
  };

public:
  void sortByStart();
  AnnotRange getStrand(char c) const;
  RecordVectorIter encloseStart(int64_t z) const;
};

using AnnotationTokenList = std::vector<gffToken>;

// For an edge to be annotated with a colour, that annotation must have started
// at least `k` positions ago.
class AnnotatedSequenceEdgeIterator {
public:
  AnnotatedSequenceEdgeIterator(Dna4SequenceConstIter seqIter, int64_t startPos,
                                std::size_t windowSize, uint64_t nullFeatureId,
                                AnnotationTokenVectorConstIter annotIter,
                                AnnotationTokenVectorConstIter annotEnd);
  AnnotatedSequenceEdgeIterator(Dna4SequenceConstIter seqIter,
                                std::size_t windowSize, uint64_t nullFeatureId,
                                AnnotationTokenVectorConstIter annotIter,
                                AnnotationTokenVectorConstIter annotEnd)
      : AnnotatedSequenceEdgeIterator(seqIter, 0, windowSize, nullFeatureId,
                                      annotIter, annotEnd) {}
  AnnotatedSequenceEdgeIterator(Dna4SequenceConstIter seqIter)
      : inputSeqIter(seqIter) {}

protected:
  Dna4SequenceConstIter inputSeqIter;
  int64_t inputSeqPos;
  int64_t k;
  uint64_t nullFeatureId;
  AnnotationTokenVectorConstIter inputAnnotIter, inputAnnotEnd;
  AnnotationTokenList currentAnnots;
  std::size_t numCurrentAnnots;

public:
  bool operator==(const AnnotatedSequenceEdgeIterator &other) const noexcept {
    return inputSeqIter == other.inputSeqIter;
  }
  bool operator!=(const AnnotatedSequenceEdgeIterator &other) const noexcept {
    return inputSeqIter != other.inputSeqIter;
  }

  void operator++();
  std::iter_value_t<Dna4SequenceConstIter> operator*() const {
    return *inputSeqIter;
  }

  uint64_t getNullFeatureId() const { return nullFeatureId; }
  bool positionIsAnnotated() const { return numCurrentAnnots != 0; }
  std::size_t getNumCurrentAnnots() const {
    return std::max(numCurrentAnnots, (std::size_t)1);
  }
  const AnnotationTokenList &viewCurrentAnnots() const { return currentAnnots; }
  std::string currentAnnotsToString() const;
};

// class FeatureSegment final : public SequenceContainer {
// public:
//   FeatureSegment(Dna4SequenceConstIter seq_, size_t begin_, size_t end_,
//                  uint64_t id_, bool terminal = false);
//   Dna4SequenceConstIter begin() const { return _it + _begin; }
//   Dna4SequenceConstIter end() const { return _it + _end; }
//   std::size_t size() const noexcept { return _end - _begin; }
//   std::size_t numInternalKmers(std::size_t k) const noexcept {
//     assert(size() >= k);
//     return size() - k;
//   }
//   std::size_t numKmers(std::size_t k) const noexcept;
//   std::size_t numTerminals(std::size_t k) const noexcept;
//   std::size_t numEdges(std::size_t k) const noexcept;
//   bool endIsTerminal() const noexcept { return _endIsTerminal; }
//   void setEndToTerminal() { _endIsTerminal = true; }
//
// protected:
//   Dna4SequenceConstIter _it;
//   std::size_t _begin;
//   std::size_t _end;
//   uint64_t _id;
//   bool _endIsTerminal;
// };

class AnnotatedSequenceSegmentIterator {
public:
  AnnotatedSequenceSegmentIterator() = default;
  AnnotatedSequenceSegmentIterator(Dna4SequenceConstIter seqIter,
                                   std::size_t seqLength, std::size_t k_,
                                   uint64_t nullFeatureId,
                                   AnnotationTokenVectorConstIter annotIter,
                                   AnnotationTokenVectorConstIter annotEnd);
  AnnotatedSequenceSegmentIterator(AnnotatedSequenceSegmentIterator &&) =
      default;
  AnnotatedSequenceSegmentIterator &
  operator=(AnnotatedSequenceSegmentIterator &&) = default;

public:
  using value_type = SequenceFragment;
  using difference_type = std::ptrdiff_t;
  using iterator_category = std::input_iterator_tag;
  using reference = SequenceFragment;

protected:
  Dna4SequenceConstIter inputSeqIter;
  int64_t inputEdgePos;
  int64_t inputSeqLength;
  int64_t k;

  uint64_t nullFeatureId;
  AnnotationTokenVectorConstIter inputAnnotIter;
  AnnotationTokenVectorConstIter inputAnnotEnd;

  int64_t currentSegmentBegin = 0;
  int64_t currentSegmentEnd = 0;
  uint64_t currentSegmentFeatureId = 0;

public:
  SequenceFragment operator*() const;

  AnnotatedSequenceSegmentIterator &operator++();
  AnnotatedSequenceSegmentIterator operator++(int);

  friend bool operator==(const AnnotatedSequenceSegmentIterator &,
                         const std::default_sentinel_t &);
  friend bool operator==(const std::default_sentinel_t &,
                         const AnnotatedSequenceSegmentIterator &);
  friend bool operator!=(const AnnotatedSequenceSegmentIterator &,
                         const std::default_sentinel_t &);
  friend bool operator!=(const std::default_sentinel_t &,
                         const AnnotatedSequenceSegmentIterator &);

  // does not count terminal edge
  inline std::size_t segmentNumInternalEdges() const noexcept {
    return static_cast<std::size_t>(currentSegmentEnd - currentSegmentBegin) -
           k;
  }
};

struct AnnotatedSequence {
  Dna4Sequence sequence;
  AnnotationTokenVector annotations;
  uint64_t nullFeatureId;

  AnnotatedSequence(const Dna4Sequence &seq)
      : sequence(seq), annotations(), nullFeatureId(0) {}
  AnnotatedSequence(const Dna4Sequence &seq, uint64_t id)
      : sequence(seq), annotations(), nullFeatureId(id) {}
  AnnotatedSequence(Dna4Sequence &&seq)
      : sequence(std::move(seq)), annotations(), nullFeatureId(0) {}
  AnnotatedSequence(Dna4Sequence &&seq, uint64_t id)
      : sequence(std::move(seq)), annotations(), nullFeatureId(id) {}
  AnnotatedSequence(uint64_t id, std::size_t length)
      : sequence(), annotations(), nullFeatureId(id) {
    sequence.reserve(length);
  }

  inline AnnotatedSequenceEdgeIterator
  colouredSegmentsEdgeBegin(std::size_t windowSize) const {
    return AnnotatedSequenceEdgeIterator(sequence.cbegin(), windowSize,
                                         nullFeatureId, annotations.cbegin(),
                                         annotations.cend());
  }
  inline AnnotatedSequenceEdgeIterator
  colouredSegmentsEdgeBegin(int64_t startPos, std::size_t windowSize) const {
    return AnnotatedSequenceEdgeIterator(
        sequence.cbegin(), startPos, windowSize, nullFeatureId,
        annotations.cbegin(), annotations.cend());
  }
  inline AnnotatedSequenceEdgeIterator colouredSegmentsEdgeEnd() const {
    return AnnotatedSequenceEdgeIterator(sequence.cend());
  }
  inline AnnotatedSequenceSegmentIterator
  colouredSegmentsBegin(std::size_t windowSize) const {
    return AnnotatedSequenceSegmentIterator(
        sequence.cbegin(), sequence.size(), windowSize, nullFeatureId,
        annotations.cbegin(), annotations.cend());
  }
  inline std::default_sentinel_t colouredSegmentsEnd() const { return {}; }

  struct FragmentView {
    const AnnotatedSequence &seq;
    std::size_t k;

    using iterator = AnnotatedSequenceSegmentIterator;
    using sentinel = std::default_sentinel_t;

    iterator begin() const { return seq.colouredSegmentsBegin(k); }
    sentinel end() const { return {}; }
  };
  static_assert(std::input_iterator<FragmentView::iterator>);
  static_assert(
      std::sentinel_for<FragmentView::sentinel, FragmentView::iterator>);
  static_assert(std::ranges::input_range<FragmentView>);

  FragmentView fragments(std::size_t k_) const;
  std::size_t numKmers(std::size_t k) const;
  SequenceFragment view() const;
};

using SequenceVector = std::vector<AnnotatedSequence>;

struct Dna4Contig {
  std::pair</*forward*/ SequenceVector, /*reverse*/ SequenceVector> sequences;
  std::string accn;
  std::size_t minFragmentSize = 0;
  std::size_t totalLength = 0;
  std::size_t numAnnotations = 0;

  Dna4Contig() = default;
  Dna4Contig(std::string &&accn_, std::size_t minContigSize_)
      : accn(accn_), minFragmentSize(minContigSize_) {}
  Dna4Contig(std::string_view accn_, std::size_t minContigSize_)
      : accn(accn_), minFragmentSize(minContigSize_) {}

  void insert(const seqan3::dna5_vector &, uint64_t seqFeatureId,
              AnnotRange annots, std::size_t minAnnotSize);
  inline std::size_t numFragments() const noexcept {
    return sequences.first.size() + sequences.second.size();
  }

  auto view() const;
  std::size_t numKmers(std::size_t) const;

protected:
  // Insert view of annotated fragment.
  void _insertFragment(SequenceVector &vec,
                       std::ranges::random_access_range auto &&sequence,
                       int64_t startPos, int64_t endPos, uint64_t seqFeatureId,
                       AnnotRange &annots, AnnotationTokenList &currentAnnots,
                       std::size_t minAnnotSize);

  // Insert annotated contig.
  void _insertOrientation(SequenceVector &vec,
                          std::ranges::random_access_range auto &&sequence,
                          uint64_t seqFeatureId,
                          const std::vector<int64_t> &nPositions,
                          AnnotRange annots, std::size_t minAnnotSize);
};

struct Dna4Genome final : public SequenceContainer {
  std::vector<Dna4Contig> contigs;

  virtual std::size_t numTerminals(std::size_t k) const override final;
  virtual poly_input_range<SequenceFragment> terminals() const override final;
  virtual std::size_t numKmers(std::size_t k) const override final;
  virtual poly_input_range<SequenceFragment>
  fragments(std::size_t k) const override final;

  std::size_t length() const;
  std::size_t numFragments() const;
  std::size_t medianContigSize() const;
};

// ---------------------------------------------------------------------------
// Parsing large FNA files (filters)
// ---------------------------------------------------------------------------

struct RestrictedSequenceFragment final : public SequenceContainer {
  RestrictedSequenceFragment(Dna4SequenceConstIter it, std::size_t begin,
                             std::size_t end, uint64_t id)
      : _it(it), _begin(begin), _end(end), _id(id) {}

  inline std::size_t size() const noexcept { return _end - _begin; }
  inline void setEndToTerminal() { _terminal = true; }

  /**
   * Filter sequences do not contribute new k-mers, so therefore contribute no
   * terminals.
   */
  virtual constexpr std::size_t
  numTerminals([[maybe_unused]] std::size_t k) const override final {
    return 0u;
  }
  virtual constexpr poly_input_range<SequenceFragment>
  terminals() const override final {
    return {};
  }
  virtual std::size_t numKmers(std::size_t k) const override final;
  virtual poly_input_range<SequenceFragment>
  fragments(std::size_t k) const override final;

protected:
  Dna4SequenceConstIter _it;
  std::size_t _begin;
  std::size_t _end;
  uint64_t _id;
  bool _terminal = false;
};

struct ChunkedDna4Genome {
  Dna4Genome genome;
  std::vector<RestrictedSequenceFragment> chunks;

  void chunk(std::size_t granularity, std::size_t overlap);
};

std::size_t chunks(const std::vector<ChunkedDna4Genome> &genomes);

// ---------------------------------------------------------------------------
// API
// ---------------------------------------------------------------------------

void parseAnnotationStream(GenomeAnnotationList &annots, Colours &colours,
                           zstr::ifstream &gffStream);
void parseFastaStream(Dna4Genome &genome, std::istream &fastaStream,
                      std::size_t minContigSize);
void parseFastaStream(Dna4Genome &genome, std::istream &fastaStream,
                      Colours &colours, std::size_t minContigSize);

Dna4Genome parseGFF(const std::string &gffFile, Colours &colours,
                    std::size_t k);
ChunkedDna4Genome parseFilterFNA(const std::string &fastaFile, Colours &colours,
                                 std::size_t granularity, std::size_t k);
std::vector<const SequenceContainer *>
combineViews(const std::vector<Dna4Genome> &,
             const std::vector<ChunkedDna4Genome> &);

template <typename ReturnType, typename... Args>
std::vector<ReturnType> parse(const std::vector<std::string> &files,
                              ReturnType (*fn)(const std::string &, Args...),
                              Args &&...args) {
  std::size_t N = files.size();
  std::vector<ReturnType> result(N);

  indicators::ProgressBar pbar{
      indicators::option::BarWidth{50},
      indicators::option::Start{"["},
      indicators::option::Fill{"="},
      indicators::option::Lead{">"},
      indicators::option::Remainder{" "},
      indicators::option::End{"]"},
      indicators::option::PrefixText{" Parsing genomes "},
      indicators::option::ForegroundColor{indicators::Color::green},
      indicators::option::ShowElapsedTime{true},
      indicators::option::ShowRemainingTime{true},
      indicators::option::FontStyles{
          std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
      indicators::option::MaxProgress{N}};

  tbb::parallel_for((std::size_t)0, N, (std::size_t)1, [&](std::size_t i) {
    result[i] = fn(files[i], args...);
    pbar.tick();
  });

  return result;
}
