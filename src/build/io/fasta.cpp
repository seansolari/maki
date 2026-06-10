#include "maki/build/io/fasta.hpp"
#include "maki/core/seq/io.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <limits>
#include <numeric>
#include <optional>
#include <ranges>
#include <string_view>

#include <oneapi/tbb.h>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/io/sequence_file/all.hpp>

bool gffToken::operator==(const gffToken &rhs) const {
  return (featureId == rhs.featureId) && (begin == rhs.begin) &&
         (end == rhs.end);
}

bool gffToken::operator!=(const gffToken &rhs) const {
  return (featureId != rhs.featureId) || (begin != rhs.begin) ||
         (end != rhs.end);
}

std::ostream &operator<<(std::ostream &os, const gffToken &tkn) {
  return os << "GFF-token[" << tkn.featureId << ", " << tkn.begin << ", "
            << tkn.end << "]";
}

// Lightweight helpers
namespace detail {

std::string_view trim(std::string_view sv) {
  std::size_t b = 0, e = sv.size();
  while (b < e &&
         (sv[b] == ' ' || sv[b] == '\t' || sv[b] == '\r' || sv[b] == '\n'))
    ++b;
  while (e > b && (sv[e - 1] == ' ' || sv[e - 1] == '\t' || sv[e - 1] == '\r' ||
                   sv[e - 1] == '\n'))
    --e;
  return sv.substr(b, e - b);
}

bool splitTab9(std::string_view line, std::array<std::string_view, 9> &out) {
  std::size_t pos = 0;
  for (int i = 0; i < 8; ++i) {
    std::size_t tab = line.find('\t', pos);
    if (tab == std::string_view::npos)
      return false;
    out[i] = line.substr(pos, tab - pos);
    pos = tab + 1;
  }
  out[8] = line.substr(pos);
  return true;
}

// parse a positive integer into int64 (no locale, no spaces)
bool parseInt64(std::string_view sv, int64_t &v) {
  sv = trim(sv);
  if (sv.empty())
    return false;
  int64_t x = 0;
  for (char c : sv) {
    if (c < '0' || c > '9')
      return false;
    int d = c - '0';
    if (x > (std::numeric_limits<std::int64_t>::max() - d) / 10)
      return false; // overflow
    x = x * 10 + d;
  }
  v = x;
  return true;
}

// Finds the value for a key in attributes field (GFF3 style key=value; pairs).
// Returns empty optional if not found. Does not unescape percent-encoding.
std::optional<std::string_view> findAttrValue(std::string_view attrs,
                                              std::string_view key) {
  // Attributes are semicolon-separated; whitespace can appear around
  // separators.
  std::size_t start = 0;
  while (start < attrs.size()) {
    // find end of this attribute
    std::size_t semi = attrs.find(';', start);
    std::string_view kv = (semi == std::string_view::npos)
                              ? attrs.substr(start)
                              : attrs.substr(start, semi - start);
    kv = trim(kv);
    if (!kv.empty()) {
      // Check if kv starts with key and '=' after optional spaces
      // Strategy: find '=' and compare trimmed key part
      std::size_t eq = kv.find('=');
      if (eq != std::string_view::npos) {
        // left and right of '='
        std::string_view k = trim(kv.substr(0, eq));
        if (k == key) {
          std::string_view v = kv.substr(eq + 1);
          // do not trim right side aggressively; only strip surrounding spaces
          v = trim(v);
          if (!v.empty())
            return v;
          else
            return std::nullopt;
        }
      }
    }
    if (semi == std::string_view::npos)
      break;
    start = semi + 1;
  }
  return std::nullopt;
}

} // namespace detail

/// Parses one GFF line.
/// Throws std::invalid_argument on hard structural errors (e.g., not 9 fields,
/// invalid coordinates). Returns GffParsed with `ignore=true` if the attributes
/// contain no ID. `id` will be "<accession>-<ID>" when ID is present, otherwise
/// empty.
std::pair<gffRecord, std::optional<std::string>>
parseGffLine(std::string_view line) {
  using namespace detail;

  // Skip comments and blank lines quickly
  std::string_view raw = trim(line);
  if (raw.empty() || (!raw.empty() && raw[0] == '#')) {
    // Treat as ignorable without throwing
    return std::make_pair(gffRecord{/*accn*/ "", /*featureId*/ 0, /*begin*/ 0,
                                    /*end*/ 0, /*strand*/ '.'},
                          std::nullopt);
  }

  std::array<std::string_view, 9> f{};
  if (!splitTab9(raw, f)) {
    throw std::invalid_argument(
        "GFF parse error: expected 9 tab-separated fields.");
  }

  std::string_view seqid = trim(f[0]); // accession
  // f[1] source (ignored)
  // f[2] type   (ignored)
  std::string_view start_s = trim(f[3]);
  std::string_view end_s = trim(f[4]);
  std::string_view strand_s = trim(f[6]);
  std::string_view attrs = trim(f[8]);

  if (seqid.empty()) {
    throw std::invalid_argument("GFF parse error: empty seqid (field 1).");
  }

  int64_t beg = 0, ed = 0;
  if (!parseInt64(start_s, beg) || !parseInt64(end_s, ed) || beg <= 0 ||
      ed <= 0) {
    throw std::invalid_argument(
        "GFF parse error: start/end must be positive integers.");
  }
  if (beg > ed) {
    // Some pipelines accept this and swap; we enforce correctness:
    // You can relax by swapping if desired.
    throw std::invalid_argument("GFF parse error: start > end.");
  } else {
    // convert from 1-based to 0-based
    beg -= 1;
  }

  char strand = '.';
  if (!strand_s.empty()) {
    char c = strand_s[0];
    if (c == '+' || c == '-' || c == '.' || c == '?') {
      strand = c;
    } else {
      // be strict; you may also map anything else to '.'
      throw std::invalid_argument(
          "GFF parse error: invalid strand (expected '+', '-', '.', or '?').");
    }
  }

  // Extract ID=... from attributes (GFF3); if absent, mark ignore=true
  std::optional<std::string_view> idv = findAttrValue(attrs, "ID");
  bool ignore = !idv.has_value();

  std::string id;
  if (idv) {
    // Build "<accession>-<ID>" with one allocation
    id.reserve(seqid.size() + 1 + idv->size());
    id.append(seqid.begin(), seqid.end());
    id.push_back('-');
    id.append(idv->begin(), idv->end());
  }

  // Materialize accession as std::string for storage beyond the lifetime of
  // input view
  gffRecord out;
  out.accn = seqid;
  out.begin = beg;
  out.end = ed;
  out.strand = strand;

  if (ignore)
    return std::make_pair(std::move(out), std::nullopt);
  else
    return std::make_pair(std::move(out), std::move(id));
}

bool operator==(const gffRecord &lhs, const gffRecord &rhs) {
  return (lhs.accn == rhs.accn) && (lhs.featureId == rhs.featureId) &&
         (lhs.begin == rhs.begin) && (lhs.end == rhs.end) &&
         (lhs.strand == rhs.strand);
}

bool operator!=(const gffRecord &lhs, const gffRecord &rhs) {
  return (lhs.accn != rhs.accn) || (lhs.featureId != rhs.featureId) ||
         (lhs.begin != rhs.begin) || (lhs.end != rhs.end) ||
         (lhs.strand != rhs.strand);
}

bool GenomeAnnotationList::AnnotLess::operator()(const gffRecord &lhs,
                                                 const gffRecord &rhs) const {
  if (lhs.accn != rhs.accn)
    return lhs.accn < rhs.accn;
  else if (lhs.strand != rhs.strand)
    return lhs.strand < rhs.strand;
  else
    return lhs.begin < rhs.begin;
}

bool GenomeAnnotationList::ContigStrandLess::operator()(
    const gffRecord &lhs, const gffRecord &rhs) const {
  if (lhs.accn != rhs.accn)
    return lhs.accn < rhs.accn;
  else
    return lhs.strand < rhs.strand;
}

// Sort by (in order) accession, strand and start position.
void GenomeAnnotationList::sort() {
  std::sort(records.begin(), records.end(), AnnotLess{});
}

// Sort by (in order) accession and strand.
void GenomeAnnotationList::sortContigStrand() {
  std::sort(records.begin(), records.end(), ContigStrandLess{});
}

// Get range of annotations correponding to a contig (all strands).
AnnotRange GenomeAnnotationList::findContig(const std::string &id) {
  return {std::lower_bound(records.begin(), records.end(), id, ContigIdLess{}),
          std::upper_bound(records.begin(), records.end(), id, ContigIdLess{})};
}

// Sort range by start position.
void AnnotRange::sortByStart() { std::sort(begin, end, BeginPosLess{}); }

// Get annotations from a strand (assumes sorted by strand).
AnnotRange AnnotRange::getStrand(char c) const {
  return {std::lower_bound(begin, end, c, StrandLess{}),
          std::upper_bound(begin, end, c, StrandLess{})};
}

// Get the first annotation that starts after `z`.
RecordVectorIter AnnotRange::encloseStart(int64_t z) const {
  return std::lower_bound(begin, end, z, BeginPosLess{});
}

AnnotatedSequenceEdgeIterator::AnnotatedSequenceEdgeIterator(
    Dna4SequenceConstIter seqIter, int64_t startPos, std::size_t windowSize,
    uint64_t nullFeatureId, AnnotationTokenVectorConstIter annotIter,
    AnnotationTokenVectorConstIter annotEnd)
    : inputSeqIter(seqIter + startPos), inputSeqPos(startPos),
      k(static_cast<int64_t>(windowSize)), nullFeatureId(nullFeatureId),
      inputAnnotIter(annotIter), inputAnnotEnd(annotEnd), currentAnnots(),
      numCurrentAnnots(0) {
  if (inputSeqPos >= k) {
    while ((inputAnnotIter != inputAnnotEnd) &&
           (inputAnnotIter->begin <= inputSeqPos - k)) {
      currentAnnots.push_back(*inputAnnotIter);
      ++numCurrentAnnots;
      ++inputAnnotIter;
    }

    std::size_t numRemoved = std::erase_if(
        currentAnnots, [&inputSeqPos = inputSeqPos](const gffToken &rec) {
          return rec.end <= inputSeqPos;
        });
    numCurrentAnnots -= numRemoved;
  }
}

void AnnotatedSequenceEdgeIterator::operator++() {
  // increment sequence position

  ++inputSeqIter;
  ++inputSeqPos;

  // remove any annotations that ended before this contig

  std::size_t numRemoved = std::erase_if(
      currentAnnots, [&inputSeqPos = inputSeqPos](const gffToken &rec) {
        return rec.end == inputSeqPos;
      });
  numCurrentAnnots -= numRemoved;

  // get new annotations that started within this region

  if (inputSeqPos >= k) {
    while ((inputAnnotIter != inputAnnotEnd) &&
           (inputAnnotIter->begin == inputSeqPos - k)) {
      currentAnnots.push_back(*inputAnnotIter);
      ++numCurrentAnnots;
      ++inputAnnotIter;
    }
  }
}

std::string AnnotatedSequenceEdgeIterator::currentAnnotsToString() const {
  std::stringstream data;
  for (const auto &annot : currentAnnots) {
    if (data.tellp() > 0) {
      data << ", ";
    }
    data << annot;
  }
  return data.str();
}

// FeatureSegment::FeatureSegment(Dna4SequenceConstIter seq_, size_t begin_,
//                                size_t end_, uint64_t id_, bool terminal_)
//     : _it(seq_), _begin(begin_), _end(end_), _id(id_),
//       _endIsTerminal(terminal_) {}
//
// std::size_t FeatureSegment::numKmers(std::size_t k) const noexcept {
//   size_t kmers = numInternalKmers(k);
//   if (endIsTerminal())
//     ++kmers;
//   return kmers;
// }
//
// std::size_t FeatureSegment::numTerminals(std::size_t k) const noexcept {
//   if (_begin == 0) {
//     assert(_end >= k);
//     return k;
//   } else
//     return 0;
// }
//
// std::size_t FeatureSegment::numEdges(std::size_t k) const noexcept {
//   return numKmers(k) + numTerminals(k);
// }

AnnotatedSequenceSegmentIterator::AnnotatedSequenceSegmentIterator(
    Dna4SequenceConstIter seqIter, std::size_t seqLength, std::size_t k_,
    uint64_t nullFeatureId, AnnotationTokenVectorConstIter annotIter,
    AnnotationTokenVectorConstIter annotEnd)
    : inputSeqIter(seqIter), inputEdgePos(static_cast<int64_t>(k_)),
      inputSeqLength(static_cast<int64_t>(seqLength)),
      k(static_cast<int64_t>(k_)), nullFeatureId(nullFeatureId),
      inputAnnotIter(annotIter), inputAnnotEnd(annotEnd) {
  operator++();
}

bool operator==(const AnnotatedSequenceSegmentIterator &obj,
                const std::default_sentinel_t &) {
  return (obj.currentSegmentBegin == obj.inputSeqLength) &&
         (obj.currentSegmentEnd == obj.inputSeqLength);
}
bool operator==(const std::default_sentinel_t &snt,
                const AnnotatedSequenceSegmentIterator &obj) {
  return obj == snt;
}
bool operator!=(const AnnotatedSequenceSegmentIterator &obj,
                const std::default_sentinel_t &snt) {
  return !(obj == snt);
}
bool operator!=(const std::default_sentinel_t &snt,
                const AnnotatedSequenceSegmentIterator &obj) {
  return !(obj == snt);
}

SequenceFragment AnnotatedSequenceSegmentIterator::operator*() const {
  return SequenceFragment(inputSeqIter + currentSegmentBegin,
                          inputSeqIter + currentSegmentEnd,
                          currentSegmentFeatureId,
                          (currentSegmentEnd == inputSeqLength) &&
                              (currentSegmentFeatureId == nullFeatureId));
}

AnnotatedSequenceSegmentIterator &
AnnotatedSequenceSegmentIterator::operator++() {
  if (inputAnnotIter != inputAnnotEnd) {
    if (inputAnnotIter->begin + k > inputEdgePos) {
      currentSegmentBegin = inputEdgePos - k;
      inputEdgePos = inputAnnotIter->begin + k;
      currentSegmentEnd = inputEdgePos;
      currentSegmentFeatureId = /*unannotated*/ nullFeatureId;
    } else {
      currentSegmentBegin = inputAnnotIter->begin;
      currentSegmentEnd = inputAnnotIter->end;
      currentSegmentFeatureId = inputAnnotIter->featureId;
      inputEdgePos = std::max(inputEdgePos, currentSegmentEnd);
      ++inputAnnotIter;
    }
  } else if (inputEdgePos <= inputSeqLength) {
    currentSegmentBegin = inputEdgePos - k;
    currentSegmentEnd = inputSeqLength;
    currentSegmentFeatureId = /*unannotated*/ nullFeatureId;
    inputEdgePos = inputSeqLength + 1;
  } else {
    currentSegmentBegin = currentSegmentEnd = inputSeqLength;
  }
  return *this;
}

std::size_t AnnotatedSequence::numKmers(std::size_t k) const {
  std::size_t edgeCount = 0;

  for (auto it = colouredSegmentsBegin(k); it != colouredSegmentsEnd(); ++it) {
    edgeCount += it.segmentNumInternalEdges();
  }

  return edgeCount + 1 /* add terminal */;
}

AnnotatedSequence::FragmentView
AnnotatedSequence::fragments(std::size_t k_) const {
  return FragmentView{*this, k_};
}

SequenceFragment AnnotatedSequence::view() const {
  return SequenceFragment(sequence.cbegin(), sequence.cend(), nullFeatureId,
                          true);
}

namespace {

std::vector<int64_t> findNs(const seqan3::dna5_vector &sequence) {
  std::vector<int64_t> positions;
  int64_t size = static_cast<int64_t>(sequence.size());
  positions.reserve(static_cast<std::size_t>(size) / 1000);

  for (int64_t i = 0; i < size; ++i) {
    if (sequence[i] == 'N'_dna5) {
      positions.push_back(i);
    }
  }

  positions.shrink_to_fit();
  return positions;
}

} // namespace

void Dna4Contig::_insertOrientation(
    SequenceVector &vec, std::ranges::random_access_range auto &&sequence,
    uint64_t seqFeatureId, const std::vector<int64_t> &nPositions,
    AnnotRange annots, std::size_t minAnnotSize) {
  // sort annotations by start position
  annots.sortByStart();
  // insert fragments
  AnnotationTokenList currentAnnots;
  int64_t seqPos = 0;
  for (int64_t nPos : nPositions) {
    std::size_t nextContigSize = nPos - seqPos;
    if (nextContigSize >= minFragmentSize)
      _insertFragment(vec, sequence, seqPos, nPos, seqFeatureId, annots,
                      currentAnnots, minAnnotSize);
    seqPos += nextContigSize + 1;
  }
  if (std::ranges::size(sequence) - seqPos >= minFragmentSize)
    _insertFragment(vec, sequence, seqPos, std::ranges::size(sequence),
                    seqFeatureId, annots, currentAnnots, minAnnotSize);
}

void Dna4Contig::_insertFragment(
    SequenceVector &vec, std::ranges::random_access_range auto &&sequence,
    int64_t startPos, int64_t endPos, uint64_t seqFeatureId, AnnotRange &annots,
    AnnotationTokenList &currentAnnots, std::size_t minAnnotSize) {
  std::size_t fragmentSize = endPos - startPos;
  totalLength += fragmentSize;
  AnnotatedSequence &seqObj = vec.emplace_back(seqFeatureId, fragmentSize);
  for (auto &&c_ : std::ranges::subrange(sequence.begin() + startPos,
                                         sequence.begin() + endPos)) {
    seqObj.sequence.push_back(static_cast<seqan3::dna4>(c_));
  }
  // remove any annotations that ended before this contig
  std::erase_if(currentAnnots, [startPos](const gffToken &rec) {
    return rec.end <= startPos;
  });
  // get new annotations that started within this region
  auto annotEnd = annots.encloseStart(endPos);
  while (annots.begin != annotEnd) {
    const gffRecord &nextAnnot = *annots.begin;
    if (nextAnnot.end > startPos)
      currentAnnots.push_back(nextAnnot.asToken());
    ++annots.begin;
  }
  // add annotations to contig
  for (const gffToken &annotToPush : currentAnnots) {
    gffToken nextAnnot = annotToPush;
    nextAnnot.begin =
        nextAnnot.begin < startPos ? 0 : nextAnnot.begin - startPos;
    nextAnnot.end =
        nextAnnot.end > endPos ? fragmentSize : nextAnnot.end - startPos;
    if (nextAnnot.size() >= minAnnotSize)
      seqObj.annotations.emplace_back(std::move(nextAnnot));
  }
}

void Dna4Contig::insert(const seqan3::dna5_vector &sequence,
                        uint64_t seqFeatureId, AnnotRange annots,
                        std::size_t minAnnotSize) {
  int64_t seqLength = static_cast<int64_t>(sequence.size());
  auto nPositions = ::findNs(sequence);

  // insert forward
  AnnotRange fwdAnnots = annots.getStrand('+');
  _insertOrientation(sequences[0], sequence, seqFeatureId, nPositions,
                     fwdAnnots, minAnnotSize);

  // get reverse complement view on sequence
  std::reverse(nPositions.begin(), nPositions.end());
  for (auto &c : nPositions)
    c = seqLength - c - 1;

  // translate reverse annotation positions to complement index
  AnnotRange revAnnots = annots.getStrand('-');
  for (auto it = revAnnots.begin; it != revAnnots.end; ++it) {
    gffRecord &annot = *it;
    std::swap(annot.begin, annot.end);
    annot.begin = annot.begin > seqLength ? 0 : seqLength - annot.begin;
    annot.end = annot.end > seqLength ? 0 : seqLength - annot.end;
  }

  // insert reverse
  _insertOrientation(sequences[1],
                     sequence | std::views::reverse | seqan3::views::complement,
                     seqFeatureId, nPositions, revAnnots, minAnnotSize);

  // revert annotation positions - important to maintain consistency with
  // original annotation file
  for (auto it = revAnnots.begin; it != revAnnots.end; ++it) {
    gffRecord &annot = *it;
    std::swap(annot.begin, annot.end);
    annot.begin = seqLength - annot.begin;
    annot.end = seqLength - annot.end;
  }
}

std::size_t Dna4Contig::numKmers(std::size_t k) const {
  auto accumulateKmers = [k](std::size_t total, const AnnotatedSequence &seq) {
    return total + seq.numKmers(k);
  };
  std::size_t total =
      std::accumulate(sequences[0].cbegin(), sequences[0].cend(),
                      (std::size_t)0, accumulateKmers);
  return std::accumulate(sequences[1].cbegin(), sequences[1].cend(),
                         total, accumulateKmers);
}

std::size_t Dna4Genome::numTerminals(std::size_t k) const {
  return k * numFragments();
}

poly_input_range<SequenceFragment> Dna4Genome::terminals() const {
  return poly_input_range<SequenceFragment>(
      contigs |
      std::views::transform([](const Dna4Contig &ctg) { return ctg.view(); }) |
      std::views::join |
      std::views::transform(
          [](const AnnotatedSequence &seq) { return seq.view(); }));
}

std::size_t Dna4Genome::numKmers(std::size_t k) const {
  auto accumulateKmers = [k](auto total, const Dna4Contig &contig) {
    return total + contig.numKmers(k);
  };
  return std::accumulate(contigs.cbegin(), contigs.cend(), (std::size_t)0,
                         accumulateKmers);
}

poly_input_range<SequenceFragment> Dna4Genome::fragments(std::size_t k) const {
  return poly_input_range<SequenceFragment>(
      contigs |
      std::views::transform([](const Dna4Contig &ctg) { return ctg.view(); }) |
      std::views::join |
      std::views::transform(
          [k](const AnnotatedSequence &seq) { return seq.fragments(k); }) |
      std::views::join);
}

std::size_t Dna4Genome::length() const {
  auto accumulateLength = [](auto total, const Dna4Contig &contig) {
    return total + contig.totalLength;
  };
  return std::accumulate(contigs.cbegin(), contigs.cend(), (std::size_t)0,
                         accumulateLength);
}

std::size_t Dna4Genome::numFragments() const {
  auto accumulateFragments = [](std::size_t total,
                                const Dna4Contig &contig) -> std::size_t {
    return total + contig.numFragments();
  };
  return std::accumulate(contigs.cbegin(), contigs.cend(), (std::size_t)0,
                         accumulateFragments);
}

std::size_t Dna4Genome::medianContigSize() const {
  std::vector<size_t> contigSizes;

  for (Dna4Contig const &ctg : contigs) {
    for (SequenceVector const *drn :
         {&ctg.sequences[0], &ctg.sequences[1]}) {
      for (AnnotatedSequence const &rec : *drn) {
        contigSizes.push_back(rec.sequence.size());
      }
    }
  }

  std::size_t medianIndex = contigSizes.size() / 2u;
  std::nth_element(contigSizes.begin(), contigSizes.begin() + medianIndex,
                   contigSizes.end());

  return contigSizes[medianIndex];
}

std::size_t RestrictedSequenceFragment::numTerminals(std::size_t k) const {
  if (_begin == 0)
    return k;
  else
    return 0;
}

poly_input_range<SequenceFragment>
RestrictedSequenceFragment::terminals() const {
  if (_begin == 0)
    return poly_input_range<SequenceFragment>(
        SequenceFragment(_it + _begin, _it + _end, _id, _terminal));
  else
    return poly_input_range<SequenceFragment>::empty_range();
}

std::size_t RestrictedSequenceFragment::numKmers(std::size_t k) const {
  return size() - k + 1;
}

poly_input_range<SequenceFragment>
RestrictedSequenceFragment::fragments([[maybe_unused]] std::size_t k) const {
  return poly_input_range<SequenceFragment>(
      SequenceFragment(_it + _begin, _it + _end, _id, _terminal));
}

void ChunkedDna4Genome::chunk(std::size_t granularity, std::size_t overlap) {
  if (!chunks.empty())
    chunks.clear();

  std::size_t ntChunkSize =
      std::max(genome.medianContigSize() / granularity, overlap);

  for (const Dna4Contig &ctg : genome.contigs) {
    for (const SequenceVector *drn :
         {&ctg.sequences[0], &ctg.sequences[1]}) {
      for (const AnnotatedSequence &rec : *drn) {
        assert(rec.annotations.size() == 0 &&
               "chunking scheme can only be applied to unannotated sequences");
        assert(rec.sequence.size() >= overlap &&
               "sequence size is less than requested overlap");
        std::size_t i = overlap;
        std::size_t j;
        while (i < rec.sequence.size()) {
          j = std::min(rec.sequence.size(), i + ntChunkSize);
          chunks.emplace_back(rec.sequence.cbegin(), i - overlap, j,
                              rec.nullFeatureId);
          i = j;
        }
        chunks.back().setEndToTerminal();
      }
    }
  }
}

std::size_t chunks(const std::vector<ChunkedDna4Genome> &genomes) {
  return std::accumulate(genomes.begin(), genomes.end(), (std::size_t)0,
                         [](std::size_t rtot, const ChunkedDna4Genome &g) {
                           return rtot + g.chunks.size();
                         });
}

void parseAnnotationStream(GenomeAnnotationList &annots, Colours &colours,
                           zstr::ifstream &gffStream) {
  std::string tmpData;
  while (std::getline(gffStream, tmpData)) {
    if (tmpData.starts_with('#')) {
      if (tmpData == "##FASTA") {
        break;
      }
    } else {
      // parse GFF line
      auto [gff, seed] = parseGffLine(tmpData);
      if (seed) {
        auto &rec = annots.records.emplace_back(std::move(gff));
        // assign colour to seed (or get same colour if identical seed)
        rec.featureId = colours.getOrAssign(std::move(*seed));
      }
    }
  }
  assert(tmpData == "##FASTA");
}

/**
 * Parse FASTA stream without assigning colours to contig sequences.
 */
void parseFastaStream(Dna4Genome &genome, std::istream &fastaStream,
                      std::size_t minContigSize) {
  seqan3::sequence_file_input seqInput(fastaStream, seqan3::format_fasta{});
  auto myLengthFilter = std::views::filter([&](auto &&seq) {
    return std::ranges::size(std::forward<decltype(seq)>(seq)) >= minContigSize;
  });

  auto Dna5ToDna4 = std::views::transform(
      [](auto &&
             c) { // https://docs.seqan.de/seqan3/main_user/cookbook.html#cookbook_convert_alphabet_range
        return static_cast<seqan3::dna4>(std::forward<decltype(c)>(c));
      });

  for (auto &&rec : seqInput) {
    Dna4Contig &out = genome.contigs.emplace_back(
        rec.id().substr(0, rec.id().find(' ')), minContigSize);

    // insert forward and reverse sequences

    SequenceVector &fwdSeqs = out.sequences[0],
                   &revSeqs = out.sequences[1];

    for (auto &&subSeq :
         rec.sequence() | std::views::split('N'_dna5) | myLengthFilter) {
      out.totalLength += 2 * subSeq.size();

      // insert forward

      auto &fwd = fwdSeqs.emplace_back(/*no colour*/ 0u, subSeq.size());

      for (auto &&c : subSeq | Dna5ToDna4) {
        fwd.sequence.push_back(c);
      }

      // insert reverse

      auto &rev = revSeqs.emplace_back(/*no colour*/ 0u, subSeq.size());

      for (auto &&c : subSeq | std::views::reverse | seqan3::views::complement |
                          Dna5ToDna4) {
        rev.sequence.push_back(c);
      }
    }
  }
}

/**
 * Parse FASTA stream, assigning colours by contig accession.
 */
void parseFastaStream(Dna4Genome &genome, std::istream &fastaStream,
                      Colours &colours, std::size_t minContigSize) {
  seqan3::sequence_file_input seqInput(fastaStream, seqan3::format_fasta{});
  auto myLengthFilter = std::views::filter([&](auto &&seq) {
    return std::ranges::size(std::forward<decltype(seq)>(seq)) >= minContigSize;
  });

  auto Dna5ToDna4 = std::views::transform(
      [](auto &&
             c) { // https://docs.seqan.de/seqan3/main_user/cookbook.html#cookbook_convert_alphabet_range
        return static_cast<seqan3::dna4>(std::forward<decltype(c)>(c));
      });

  for (auto &&rec : seqInput) {
    Dna4Contig &out = genome.contigs.emplace_back(
        rec.id().substr(0, rec.id().find(' ')), minContigSize);
    uint64_t seqId = colours.getOrAssign(std::string(out.accn));

    // insert forward and reverse sequences

    SequenceVector &fwdSeqs = out.sequences[0],
                   &revSeqs = out.sequences[1];

    for (auto &&subSeq :
         rec.sequence() | std::views::split('N'_dna5) | myLengthFilter) {
      out.totalLength += 2 * subSeq.size();

      // insert forward

      auto &fwd = fwdSeqs.emplace_back(seqId, subSeq.size());

      for (auto &&c : subSeq | Dna5ToDna4) {
        fwd.sequence.push_back(c);
      }

      // insert reverse

      auto &rev = revSeqs.emplace_back(seqId, subSeq.size());

      for (auto &&c : subSeq | std::views::reverse | seqan3::views::complement |
                          Dna5ToDna4) {
        rev.sequence.push_back(c);
      }
    }
  }
}

Dna4Genome parseGFF(const std::string &gff3File, Colours &colours,
                    std::size_t k) {
  // base input file stream that reads bytes
  zstr::ifstream zis(gff3File);

  // process stream data
  GenomeAnnotationList annots;
  parseAnnotationStream(annots, colours, zis);
  annots.sortContigStrand();

  // read remaining data into string (and de-compress if necessary)
  std::string fastaData(std::istreambuf_iterator<char>(zis), {});

  // parse fasta records
  Dna4Genome genome{};
  std::istringstream fs(fastaData);
  seqan3::sequence_file_input seqInput(fs, seqan3::format_fasta{});

  for (auto &&rec : seqInput) {
    Dna4Contig &outContig =
        genome.contigs.emplace_back(rec.id().substr(0, rec.id().find(' ')), k);
    AnnotRange contigAnnots = annots.findContig(outContig.accn);
    outContig.numAnnotations = contigAnnots.size();
    outContig.insert(rec.sequence(), /*unannotated regions*/ 0ull, contigAnnots,
                     k + 1);
  }

  return genome;
}

ChunkedDna4Genome parseFilterFNA(const std::string &fastaFile, Colours &colours,
                                 std::size_t granularity, std::size_t k) {
  // base input file stream that reads bytes
  zstr::ifstream zis(fastaFile);
  std::string fastaData(std::istreambuf_iterator<char>(zis), {});

  // process stream data
  ChunkedDna4Genome obj;
  std::istringstream fs(fastaData);
  parseFastaStream(obj.genome, fs, colours, k);

  // create genome slices
  obj.chunk(granularity, k);

  return obj;
}

std::vector<const SequenceContainer *>
toView(const std::vector<Dna4Genome> &gff) {
  std::vector<const SequenceContainer *> views;
  views.reserve(gff.size());
  for (const auto &seq : gff) {
    views.push_back(&seq);
  }
  return views;
}

std::vector<const SequenceContainer *>
toView(const std::vector<ChunkedDna4Genome> &fna) {
  std::vector<const SequenceContainer *> views;
  views.reserve(chunks(fna));
  for (const auto &seq : fna) {
    for (const auto &chunk : seq.chunks) {
      views.push_back(&chunk);
    }
  }
  return views;
}
