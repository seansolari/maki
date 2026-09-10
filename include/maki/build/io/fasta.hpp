#pragma once
#include <cassert>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <optional>
#include <ranges>
#include <string>
#include <type_traits>
#include <vector>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <zstr.hpp>

#include "maki/build/graph/build_colours.hpp"
#include "maki/core/seq/concepts.hpp"
#include "maki/core/seq/io.hpp"

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
class AnnotatedEdgeIterator {
public:
  AnnotatedEdgeIterator(Dna4SequenceConstIter seqIter, int64_t startPos,
                        std::size_t windowSize, uint64_t nullFeatureId,
                        AnnotationTokenVectorConstIter annotIter,
                        AnnotationTokenVectorConstIter annotEnd);

  AnnotatedEdgeIterator(Dna4SequenceConstIter seqIter, std::size_t windowSize,
                        uint64_t nullFeatureId,
                        AnnotationTokenVectorConstIter annotIter,
                        AnnotationTokenVectorConstIter annotEnd)
      : AnnotatedEdgeIterator(seqIter, 0, windowSize, nullFeatureId, annotIter,
                              annotEnd) {}

  AnnotatedEdgeIterator(Dna4SequenceConstIter seqIter)
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
  bool operator==(const AnnotatedEdgeIterator &other) const noexcept {
    return inputSeqIter == other.inputSeqIter;
  }
  bool operator!=(const AnnotatedEdgeIterator &other) const noexcept {
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

class AnnotatedSegmentIterator {
public:
  AnnotatedSegmentIterator() = default;

  AnnotatedSegmentIterator(Dna4SequenceConstIter seqIter, std::size_t seqLength,
                           std::size_t k_, uint64_t nullFeatureId,
                           AnnotationTokenVectorConstIter annotIter,
                           AnnotationTokenVectorConstIter annotEnd);

  AnnotatedSegmentIterator(AnnotatedSegmentIterator &&) = default;

  AnnotatedSegmentIterator &operator=(AnnotatedSegmentIterator &&) = default;

public:
  using reference = SequenceFragment<
      std::ranges::subrange<seqan3::detail::random_access_iterator<
          const seqan3::bitpacked_sequence<seqan3::dna4>>>>;
  using value_type = std::remove_cvref_t<reference>;
  using difference_type = std::ptrdiff_t;
  using iterator_category = std::input_iterator_tag;

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
  value_type operator*() const;

  AnnotatedSegmentIterator &operator++();
  AnnotatedSegmentIterator operator++(int);

  friend bool operator==(const AnnotatedSegmentIterator &,
                         const std::default_sentinel_t &);
  friend bool operator==(const std::default_sentinel_t &,
                         const AnnotatedSegmentIterator &);
  friend bool operator!=(const AnnotatedSegmentIterator &,
                         const std::default_sentinel_t &);
  friend bool operator!=(const std::default_sentinel_t &,
                         const AnnotatedSegmentIterator &);

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

  inline AnnotatedEdgeIterator
  colouredSegmentsEdgeBegin(std::size_t windowSize) const {
    return AnnotatedEdgeIterator(sequence.cbegin(), windowSize, nullFeatureId,
                                 annotations.cbegin(), annotations.cend());
  }
  inline AnnotatedEdgeIterator
  colouredSegmentsEdgeBegin(int64_t startPos, std::size_t windowSize) const {
    return AnnotatedEdgeIterator(sequence.cbegin(), startPos, windowSize,
                                 nullFeatureId, annotations.cbegin(),
                                 annotations.cend());
  }
  inline AnnotatedEdgeIterator colouredSegmentsEdgeEnd() const {
    return AnnotatedEdgeIterator(sequence.cend());
  }
  inline AnnotatedSegmentIterator
  colouredSegmentsBegin(std::size_t windowSize) const {
    return AnnotatedSegmentIterator(sequence.cbegin(), sequence.size(),
                                    windowSize, nullFeatureId,
                                    annotations.cbegin(), annotations.cend());
  }
  inline std::default_sentinel_t colouredSegmentsEnd() const { return {}; }

  struct FragmentView {
    const AnnotatedSequence &seq;
    std::size_t k;

    using iterator = AnnotatedSegmentIterator;
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

  using view_type = SequenceFragment<
      std::ranges::ref_view<const seqan3::bitpacked_sequence<seqan3::dna4>>>;

  view_type view() const;
};

using SequenceVector = std::vector<AnnotatedSequence>;

struct Dna4Contig {
  std::array<SequenceVector, 2> sequences;
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
    return sequences[0].size() + sequences[1].size();
  }

  auto view() const { return sequences | std::views::join; }

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

namespace detail {

constexpr auto flatten_contig_strands =
    std::views::transform([](const Dna4Contig &ctg) { return ctg.view(); }) |
    std::views::join;

constexpr auto discard_annotations = std::views::transform(
    [](const AnnotatedSequence &seq) { return seq.view(); });

constexpr auto join_annotated_segments(std::size_t k) {
  return std::views::transform(
             [k](const AnnotatedSequence &seq) { return seq.fragments(k); }) |
         std::views::join;
}

template <std::ranges::random_access_range R>
using contig_terminal_view_type =
    decltype(std::declval<R>() | flatten_contig_strands | discard_annotations);

template <std::ranges::random_access_range R>
using contig_fragment_view_type =
    decltype(std::declval<R>() | flatten_contig_strands |
             join_annotated_segments(0));

} // namespace detail

struct Dna4Genome {
  std::vector<Dna4Contig> contigs;

  Dna4Genome() = default;
  Dna4Genome(Dna4Genome &&) = default;
  Dna4Genome &operator=(Dna4Genome &&) = default;
  Dna4Genome(const Dna4Genome &) = delete;
  Dna4Genome &operator=(const Dna4Genome &) = delete;

  std::size_t numTerminals(std::size_t k) const;
  std::size_t numKmers(std::size_t k) const;

public:
  using terminal_view_type =
      detail::contig_terminal_view_type<const decltype(contigs) &>;

  terminal_view_type terminals() const {
    return contigs | detail::flatten_contig_strands |
           detail::discard_annotations;
  }

  using fragment_view_type =
      detail::contig_fragment_view_type<const decltype(contigs) &>;

  fragment_view_type fragments(std::size_t k) const {
    return contigs | detail::flatten_contig_strands |
           detail::join_annotated_segments(k);
  }

  std::size_t length() const;
  std::size_t numFragments() const;
  std::size_t medianContigSize() const;
};

static_assert(sequence_container_like<Dna4Genome>);

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

template <typename ReturnType, typename... Args>
std::vector<ReturnType> parse(const std::vector<std::string> &files,
                              ReturnType (*fn)(const std::string &, Args...),
                              Args &&...args) {
  std::size_t N = files.size();
  std::vector<ReturnType> result(N);

  tbb::parallel_for((std::size_t)0, N, (std::size_t)1,
                    [&](std::size_t i) { result[i] = fn(files[i], args...); });

  return result;
}
