#include "maki/io/fasta.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <limits>
#include <numeric>
#include <string_view>
#include <glog/logging.h>
#include <oneapi/tbb.h>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <indicators/progress_bar.hpp>

std::string toString(Dna4Sequence const &seq)
{
  std::string result;
  result.reserve(seq.size());

  for (auto c : seq)
    result.push_back(seqan3::to_char(c));

  return result;
}

std::string toString(std::vector<Dna4Sequence> const &v, const char *delim)
{
  size_t resultSize = 0;
  for (auto const &seq : v)
    resultSize += seq.size();
  resultSize += strlen(delim) * (v.size() - 1);

  std::string result;
  result.reserve(resultSize);

  for (auto const &seq : v)
  {
    if (!result.empty())
      result.append(delim);
    for (auto c : seq)
      result.push_back(seqan3::to_char(c));
  }

  return result;
}

// Lightweight helpers
namespace detail
{

  std::string_view trim(std::string_view sv)
  {
    size_t b = 0, e = sv.size();
    while (b < e && (sv[b] == ' ' || sv[b] == '\t' || sv[b] == '\r' || sv[b] == '\n'))
      ++b;
    while (e > b && (sv[e - 1] == ' ' || sv[e - 1] == '\t' || sv[e - 1] == '\r' || sv[e - 1] == '\n'))
      --e;
    return sv.substr(b, e - b);
  }

  bool splitTab9(std::string_view line, std::array<std::string_view, 9> &out)
  {
    size_t pos = 0;
    for (int i = 0; i < 8; ++i)
    {
      size_t tab = line.find('\t', pos);
      if (tab == std::string_view::npos)
        return false;
      out[i] = line.substr(pos, tab - pos);
      pos = tab + 1;
    }
    out[8] = line.substr(pos);
    return true;
  }

  // parse a positive integer into int64 (no locale, no spaces)
  bool parseInt64(std::string_view sv, int64_t &v)
  {
    sv = trim(sv);
    if (sv.empty())
      return false;
    int64_t x = 0;
    for (char c : sv)
    {
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
  std::optional<std::string_view> findAttrValue(std::string_view attrs, std::string_view key)
  {
    // Attributes are semicolon-separated; whitespace can appear around separators.
    size_t start = 0;
    while (start < attrs.size())
    {
      // find end of this attribute
      size_t semi = attrs.find(';', start);
      std::string_view kv = (semi == std::string_view::npos) ? attrs.substr(start)
                                                             : attrs.substr(start, semi - start);
      kv = trim(kv);
      if (!kv.empty())
      {
        // Check if kv starts with key and '=' after optional spaces
        // Strategy: find '=' and compare trimmed key part
        size_t eq = kv.find('=');
        if (eq != std::string_view::npos)
        {
          // left and right of '='
          std::string_view k = trim(kv.substr(0, eq));
          if (k == key)
          {
            std::string_view v = kv.substr(eq + 1);
            // do not trim right side aggressively; only strip surrounding spaces
            v = trim(v);
            if (!v.empty())
              return v;
            // If empty value, still considered present but empty
            return std::string_view{};
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
/// Throws std::invalid_argument on hard structural errors (e.g., not 9 fields, invalid coordinates).
/// Returns GffParsed with `ignore=true` if the attributes contain no ID.
/// `id` will be "<accession>-<ID>" when ID is present, otherwise empty.
std::pair<gffRecord,std::optional<std::string>> parseGffLine(std::string_view line)
{
  using namespace detail;

  // Skip comments and blank lines quickly
  std::string_view raw = trim(line);
  if (raw.empty() || (!raw.empty() && raw[0] == '#'))
  {
    // Treat as ignorable without throwing
    return std::make_pair(gffRecord{/*accn*/ "", /*featureId*/ 0, /*begin*/ 0, /*end*/ 0, /*strand*/ '.'}, std::nullopt);
  }

  std::array<std::string_view, 9> f{};
  if (!splitTab9(raw, f))
  {
    throw std::invalid_argument("GFF parse error: expected 9 tab-separated fields.");
  }

  std::string_view seqid = trim(f[0]); // accession
  // f[1] source (ignored)
  // f[2] type   (ignored)
  std::string_view start_s = trim(f[3]);
  std::string_view end_s = trim(f[4]);
  std::string_view strand_s = trim(f[6]);
  std::string_view attrs = trim(f[8]);

  if (seqid.empty())
  {
    throw std::invalid_argument("GFF parse error: empty seqid (field 1).");
  }

  int64_t beg = 0, ed = 0;
  if (!parseInt64(start_s, beg) || !parseInt64(end_s, ed) || beg <= 0 || ed <= 0)
  {
    throw std::invalid_argument("GFF parse error: start/end must be positive integers.");
  }
  if (beg > ed)
  {
    // Some pipelines accept this and swap; we enforce correctness:
    // You can relax by swapping if desired.
    throw std::invalid_argument("GFF parse error: start > end.");
  }
  else
  {
    // convert from 1-based to 0-based
    beg -= 1;
  }

  char strand = '.';
  if (!strand_s.empty())
  {
    char c = strand_s[0];
    if (c == '+' || c == '-' || c == '.' || c == '?')
    {
      strand = c;
    }
    else
    {
      // be strict; you may also map anything else to '.'
      throw std::invalid_argument("GFF parse error: invalid strand (expected '+', '-', '.', or '?').");
    }
  }

  // Extract ID=... from attributes (GFF3); if absent, mark ignore=true
  std::optional<std::string_view> idv = findAttrValue(attrs, "ID");
  bool ignore = !idv.has_value();

  std::string id;
  if (idv)
  {
    // Build "<accession>-<ID>" with one allocation
    id.reserve(seqid.size() + 1 + idv->size());
    id.append(seqid.begin(), seqid.end());
    id.push_back('-');
    id.append(idv->begin(), idv->end());
  }

  // Materialize accession as std::string for storage beyond the lifetime of input view
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

std::ostream &operator<<(std::ostream &os, const gffToken &obj)
{
  return os
         << "gffRecord("
         << obj.accn << ", "
         << obj.begin << ", "
         << obj.end << ", "
         << obj.strand << ")";
}

bool operator==(const gffRecord &lhs, const gffRecord &rhs)
{
  return (lhs.accn == rhs.accn) &&
         (lhs.featureId == rhs.featureId) && 
         (lhs.begin == rhs.begin) &&
         (lhs.end == rhs.end) &&
         (lhs.strand == rhs.strand);
}

bool operator!=(const gffRecord &lhs, const gffRecord &rhs)
{
  return (lhs.id != rhs.id) ||
         (lhs.featureId != rhs.featureId) ||
         (lhs.begin != rhs.begin) ||
         (lhs.end != rhs.end) ||
         (lhs.strand != rhs.strand);
}

void GenomeAnnotationList::sortContigStrand()
{
  std::sort(records.begin(), records.end(), ContigStrandLess{});
}

void GenomeAnnotationList::sort()
{
  std::sort(records.begin(), records.end(), AnnotLess{});
}

void GenomeAnnotationList::ensureDisjoint()
{
  auto
      it = records.cbegin(),
      end = records.cend();

  while (it != end)
  {
    auto next = std::upper_bound(it, end, *it, ContigStrandLess{});
    ensureDisjointRegion(it, next);
    it = std::next(next);
  }
}

void GenomeAnnotationList::ensureDisjointRegion(RecordVectorConstIter it, RecordVectorConstIter end)
{
  gffRecord prev = *it++;

  while (it != end)
  {
    gffRecord curr = *it++;
    assert(prev.end <= curr.begin);
    std::swap(prev, curr);
  }
}

AnnotRange GenomeAnnotationList::findContig(const std::string &id)
{
  return {
      std::lower_bound(records.begin(), records.end(), id, ContigIdLess{}),
      std::upper_bound(records.begin(), records.end(), id, ContigIdLess{})};
}

bool ContigStrandLess::operator()(const gffRecord &lhs, const gffRecord &rhs) const
{
  if (lhs.id != rhs.id)
    return lhs.id < rhs.id;
  else
    return lhs.strand < rhs.strand;
}

bool AnnotLess::operator()(const gffRecord &lhs, const gffRecord &rhs) const
{
  if (lhs.id != rhs.id)
    return lhs.id < rhs.id;
  else if (lhs.strand != rhs.strand)
    return lhs.strand < rhs.strand;
  else
    return lhs.begin < rhs.begin;
}

void AnnotRange::sortByStart()
{
  std::sort(begin, end, BeginPosLess{});
}

AnnotRange AnnotRange::getStrand(char c) const
{
  return {
      std::lower_bound(begin, end, c, StrandLess{}),
      std::upper_bound(begin, end, c, StrandLess{})};
}

RecordVectorIter AnnotRange::encloseStart(size_t z) const
{
  return std::lower_bound(begin, end, z, BeginPosLess{});
}

void parseAnnotationStream(GenomeAnnotationList &annots, Colours &colours, zstr::ifstream &gffStream)
{
  std::string tmpData;
  while (std::getline(gffStream, tmpData))
  {
    if (tmpData.starts_with('#'))
    {
      if (tmpData == "##FASTA")
      {
        break;
      }
    }
    else
    {
      // parse GFF line
      auto [gff, seed] = parseGffLine(tmpData);
      if (seed)
      {
        auto &rec = annots.records.emplace_back(std::move(gff));
        // assign colour to seed (or get same colour if identical seed)
        rec.featureId = colours.getOrAssign(std::move(*seed));
      }
    }
  }
  assert(tmpData == "##FASTA");
}

AnnotatedSequenceEdgeIterator::AnnotatedSequenceEdgeIterator(
    Dna4SequenceConstIter seqIter,
    size_t startPos,
    uint8_t windowSize,
    AnnotationTokenVectorConstIter annotIter,
    AnnotationTokenVectorConstIter annotEnd)
    : inputSeqIter(seqIter + startPos),
      inputSeqPos(startPos),
      k(windowSize),
      inputAnnotIter(annotIter),
      inputAnnotEnd(annotEnd),
      currentAnnots(),
      numCurrentAnnots(0)
{
  if (inputSeqPos >= k)
  {
    while ((inputAnnotIter != inputAnnotEnd) && (inputAnnotIter->begin <= inputSeqPos - k))
    {
      currentAnnots.push_back(*inputAnnotIter);
      ++numCurrentAnnots;
      ++inputAnnotIter;
    }

    size_t numRemoved = std::erase_if(currentAnnots, [&inputSeqPos = inputSeqPos](const gffToken &rec)
                                      { return rec.end <= inputSeqPos; });
    numCurrentAnnots -= numRemoved;
  }
}

void AnnotatedSequenceEdgeIterator::operator++()
{
  // increment sequence position

  ++inputSeqIter;
  ++inputSeqPos;

  // remove any annotations that ended before this contig

  size_t numRemoved = std::erase_if(currentAnnots, [&inputSeqPos = inputSeqPos](const gffToken &rec)
                                    { return rec.end == inputSeqPos; });
  numCurrentAnnots -= numRemoved;

  // get new annotations that started within this region

  if (inputSeqPos >= k)
  {
    while ((inputAnnotIter != inputAnnotEnd) && (inputAnnotIter->begin == inputSeqPos - k))
    {
      currentAnnots.push_back(*inputAnnotIter);
      ++numCurrentAnnots;
      ++inputAnnotIter;
    }
  }
}

std::string AnnotatedSequenceEdgeIterator::currentAnnotsToString() const
{
  std::stringstream data;
  for (const auto &annot : currentAnnots)
  {
    if (data.tellp() > 0)
    {
      data << ", ";
    }
    data << annot;
  }
  return data.str();
}

AnnotatedSequenceSegment AnnotatedSequenceSegmentIterator::operator*() const
{
  return {
      inputSeqIter + currentSegmentBegin,
      inputSeqIter + currentSegmentEnd,
      currentSegmentFeatureId,
      (currentSegmentEnd == inputSeqLength) && (currentSegmentFeatureId == 0)};
}

void AnnotatedSequenceSegmentIterator::operator++()
{
  if (inputAnnotIter != inputAnnotEnd)
  {
    if (inputAnnotIter->begin + k > inputEdgePos)
    {
      currentSegmentBegin = inputEdgePos - k;
      inputEdgePos = inputAnnotIter->begin + k;
      currentSegmentEnd = inputEdgePos;
      currentSegmentFeatureId = /* UNANNOTATED */ 0;
    }
    else
    {
      currentSegmentBegin = inputAnnotIter->begin;
      currentSegmentEnd = inputAnnotIter->end;
      currentSegmentFeatureId = inputAnnotIter->featureId;

      inputEdgePos = std::max(inputEdgePos, currentSegmentEnd);
      ++inputAnnotIter;
    }
  }
  else if (inputEdgePos <= inputSeqLength)
  {
    currentSegmentBegin = inputEdgePos - k;
    currentSegmentEnd = inputSeqLength;
    currentSegmentFeatureId = /* UNANNOTATED */ 0;

    inputEdgePos = inputSeqLength + 1;
  }
  else
  {
    currentSegmentBegin = currentSegmentEnd = inputSeqLength;
  }
}

void parseFastaStream(Dna4Genome &genome, std::istream &fastaStream, size_t minContigSize)
{
  seqan3::sequence_file_input seqInput(fastaStream, seqan3::format_fasta{});
  auto myLengthFilter = std::views::filter(
      [&](auto &&seq)
      {
        return std::ranges::size(std::forward<decltype(seq)>(seq)) >= minContigSize;
      });

  auto Dna5ToDna4 = std::views::transform(
      [](auto &&c) { // https://docs.seqan.de/seqan3/main_user/cookbook.html#cookbook_convert_alphabet_range
        return static_cast<seqan3::dna4>(std::forward<decltype(c)>(c));
      });

  for (auto &&rec : seqInput)
  {
    Dna4Contig &out = genome.contigs.emplace_back(rec.id().substr(0, rec.id().find(' ')), minContigSize);

    // insert forward and reverse sequences

    SequenceVector
        &fwdSeqs = out.sequences.first,
        &revSeqs = out.sequences.second;

    for (auto &&subSeq : rec.sequence() | std::views::split('N'_dna5) | myLengthFilter)
    {
      out.totalLength += 2 * subSeq.size();

      // insert forward

      auto &fwd = fwdSeqs.emplace_back(subSeq.size());

      for (auto &&c : subSeq | Dna5ToDna4)
      {
        fwd.sequence.push_back(c);
      }

      // insert reverse

      auto &rev = revSeqs.emplace_back(subSeq.size());

      for (auto &&c : subSeq | std::views::reverse | seqan3::views::complement | Dna5ToDna4)
      {
        rev.sequence.push_back(c);
      }
    }
  }
}

Dna4Genome parseGenome(fs::path fastaFile, size_t minContigSize)
{
  // base input file stream that reads bytes

  zstr::ifstream zis(fastaFile);
  std::string fastaData(std::istreambuf_iterator<char>(zis), {});

  // process stream data

  Dna4Genome genome{};

  std::istringstream fs(fastaData);
  parseFastaStream(genome, fs, minContigSize);

  return genome;
}

Dna4Genome parseAnnotatedGenome(fs::path gff3File, Colours &colours, size_t minContigSize)
{
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

  for (auto &&rec : seqInput)
  {
    Dna4Contig &outContig = genome.contigs.emplace_back(rec.id().substr(0, rec.id().find(' ')), minContigSize);
    AnnotRange contigAnnots = annots.findContig(outContig.contigName);
    outContig.insert(rec.sequence(), contigAnnots, minContigSize + 1);
  }

  return genome;
}

size_t AnnotatedSequence::numKmers(uint8_t k) const
{
  size_t edgeCount = 0;

  for (auto it = colouredSegmentsBegin(k); it != colouredSegmentsEnd(); ++it)
  {
    edgeCount += it.segmentNumInternalEdges();
  }

  return edgeCount + 1 /* add terminal */;
}

void Dna4Contig::insert(const seqan3::dna5_vector &sequence, AnnotRange annots, size_t minAnnotSize)
{
  size_t seqLength = sequence.size();
  std::vector<size_t> nPositions = findNs(sequence);

  // insert forward
  AnnotRange fwdAnnots = annots.getStrand('+');
  _insertOrientation(sequences.first, sequence, nPositions, fwdAnnots, minAnnotSize);

  // get reverse complement view on sequence
  std::reverse(nPositions.begin(), nPositions.end());
  for (auto &c : nPositions)
    c = seqLength - c - 1;

  // translate reverse annotation positions to complement index
  AnnotRange revAnnots = annots.getStrand('-');
  for (auto it = revAnnots.begin; it != revAnnots.end; ++it)
  {
    gffRecord &annot = *it;
    std::swap(annot.begin, annot.end);
    annot.begin = annot.begin > seqLength ? 0 : seqLength - annot.begin;
    annot.end = annot.end > seqLength ? 0 : seqLength - annot.end;
  }

  // insert reverse
  auto revCompSequence = sequence | std::views::reverse | seqan3::views::complement;
  _insertOrientation(sequences.second, revCompSequence, nPositions, revAnnots, minAnnotSize);

  // revert annotation positions - important to maintain consistency with original annotation file
  for (auto it = revAnnots.begin; it != revAnnots.end; ++it)
  {
    gffRecord &annot = *it;
    std::swap(annot.begin, annot.end);
    annot.begin = seqLength - annot.begin;
    annot.end = seqLength - annot.end;
  }
}

std::vector<size_t> Dna4Contig::findNs(const seqan3::dna5_vector &sequence)
{
  std::vector<size_t> positions;
  positions.reserve(sequence.size() / 1000);

  for (size_t i = 0; i < sequence.size(); ++i)
  {
    if (sequence[i] == 'N'_dna5)
    {
      positions.push_back(i);
    }
  }

  positions.shrink_to_fit();
  return positions;
}

void Dna4Contig::_insertOrientation(SequenceVector &vec, R &&seqView, uint64_t seqFeatureId, const std::vector<size_t> &nPos, AnnotRange annots, size_t minAnnotSize)
{
  // sort annotations by start position
  annots.sortByStart();
  // insert fragments
  AnnotationTokenList currentAnnots;
  size_t seqPos = 0;
  std::ranges::iterator_t<R> seqIter = seqView.begin();
  for (size_t nPos : nPositions) {
      size_t nextContigSize = nPos - seqPos;
      if (nextContigSize >= minFragmentSize)
          _insertFragment(vec, seqIter, seqPos, nPos, seqFeatureId, annots, currentAnnots, minAnnotSize);
      seqPos += nextContigSize + 1;
  }
  if (seqView.size() - seqPos >= minFragmentSize)
      _insertFragment(vec, seqIter, seqPos, seqView.size(), seqFeatureId, annots, currentAnnots, minAnnotSize);
}

void Dna4Contig::_insertFragment(SequenceVector &vec, Iter seqIter, size_t startPos, size_t endPos, uint64_t seqFeatureId, AnnotRange &annots, AnnotationTokenList &currentAnnots, size_t minAnnotSize)
{
  size_t fragmentSize = endPos - startPos;
  totalLength += fragmentSize;
  AnnotatedSequence &seqObj = vec.emplace_back(seqFeatureId, fragmentSize);
  for (auto it = seqIter + startPos; it != seqIter + endPos; ++it) {
    auto c = static_cast<seqan3::dna4>(*it);
    seqObj.sequence.push_back(c);
  }
  // remove any annotations that ended before this contig
  std::erase_if(currentAnnots, [startPos](const Gff3Token &rec) { return rec.endPos <= startPos; });
  // get new annotations that started within this region
  auto annotEnd = annots.encloseStart(endPos);
  while (annots.begin != annotEnd) {
    const Gff3Record &nextAnnot = *annots.begin;
    if (nextAnnot.endPos > startPos)
      currentAnnots.push_back(nextAnnot.asToken());
    ++annots.begin;
  }
  // add annotations to contig
  for (const Gff3Token &annotToPush : currentAnnots) {
    Gff3Token nextAnnot = annotToPush;
    nextAnnot.beginPos = nextAnnot.beginPos < startPos ? 0 : nextAnnot.beginPos - startPos;
    nextAnnot.endPos = nextAnnot.endPos > endPos ? fragmentSize : nextAnnot.endPos - startPos;
    if (nextAnnot.size() >= minAnnotSize)
      seqObj.annotations.emplace_back(std::move(nextAnnot));
  }
}

size_t Dna4Contig::numKmers(uint8_t k) const
{
  auto accumulateKmers = [k](size_t total, const AnnotatedSequence &seq)
  { return total + seq.numKmers(k); };
  size_t total = std::accumulate(sequences.first.cbegin(), sequences.first.cend(), (size_t)0, accumulateKmers);
  return std::accumulate(sequences.second.cbegin(), sequences.second.cend(), total, accumulateKmers);
}

size_t Dna4Genome::length() const
{
  auto accumulateLength = [](auto total, const Dna4Contig &contig)
  { return total + contig.totalLength; };
  return std::accumulate(contigs.cbegin(), contigs.cend(), (size_t)0, accumulateLength);
}

size_t Dna4Genome::numKmers(uint8_t k) const
{
  auto accumulateKmers = [k](auto total, const Dna4Contig &contig)
  { return total + contig.numKmers(k); };
  return std::accumulate(contigs.cbegin(), contigs.cend(), (size_t)0, accumulateKmers);
}

size_t Dna4Genome::numFragments() const
{
  auto accumulateFragments = [](size_t total, const Dna4Contig &contig) -> size_t
  { return total + contig.numFragments(); };
  return std::accumulate(contigs.cbegin(), contigs.cend(), (size_t)0, accumulateFragments);
}

size_t Dna4Genome::numTerminals(uint8_t k) const
{
  return size_t{k} * numFragments();
}

size_t Dna4Genome::medianContigSize() const
{
  std::vector<size_t> contigSizes;

  for (Dna4Contig const &ctg : contigs)
  {
    for (SequenceVector const *drn : {&ctg.sequences.first, &ctg.sequences.second})
    {
      for (AnnotatedSequence const &rec : *drn)
      {
        contigSizes.push_back(rec.sequence.size());
      }
    }
  }

  size_t medianIndex = contigSizes.size() / 2u;
  std::nth_element(
      contigSizes.begin(),
      contigSizes.begin() + medianIndex,
      contigSizes.end());

  return contigSizes[medianIndex];
}

size_t Dna4Genome::rss() const
{
  size_t currentSize = 0;
  for (const Dna4Contig &ctg : contigs)
  {
    for (const SequenceVector *ornt : {&ctg.sequences.first, &ctg.sequences.second})
    {
      for (const AnnotatedSequence &rec : *ornt)
      {
        currentSize += ((2ul /*bits-per-nt Dna4Sequence*/ * rec.sequence.capacity()) + 7ul) / 8ul;
        currentSize += sizeof(gffToken) * rec.annotations.capacity();
      }
    }
  }
  return currentSize;
}

void Dna4Genome::reset_memory()
{
  std::vector<Dna4Contig>{}.swap(contigs);
}

std::vector<fs::path> readFilePaths(const char *manifest_file)
{
  std::ifstream manifest(manifest_file);
  // read lines from file
  std::vector<fs::path> files;
  std::copy(std::istream_iterator<fs::path>(manifest),
            std::istream_iterator<fs::path>(),
            std::back_inserter(files));
  // validate all these files exist
  size_t missing = 0;
  for (const fs::path &file : files)
    if (!fs::exists(file))
    {
      ++missing;
      LOG(ERROR) << file << " not found";
    }
  if (missing)
  {
    LOG(ERROR) << missing << " files not found";
    throw "input files not found";
  }
  else
  {
    LOG(INFO) << "parsed " << files.size() << " files";
  }
  return files;
}

size_t Dna4GenomeVector::totalLength() const
{
  auto accumulateLength = [](auto total, const Dna4Genome &genome)
  { return total + genome.length(); };
  return std::accumulate(genomes.cbegin(), genomes.cend(), (uint64_t)0, accumulateLength);
}

size_t Dna4GenomeVector::numKmers(uint8_t k) const
{
  auto accumulateKmers = [k](auto total, const Dna4Genome &genome)
  { return total + genome.numKmers(k); };
  return std::accumulate(genomes.cbegin(), genomes.cend(), (size_t)0, accumulateKmers);
}

Dna4GenomeVector loadGenomes(const std::vector<fs::path> &fastaFiles,
                             Colours &colours,
                             size_t minContigSize)
{
  size_t numGenomes = fastaFiles.size();
  Dna4GenomeVector genomes(numGenomes);

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
      indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
      indicators::option::MaxProgress{numGenomes}};

  tbb::parallel_for(tbb::blocked_range<size_t>(0, genomes.numGenomes(), 5),
                    [&](tbb::blocked_range<size_t> &r)
                    {
                      for (size_t i = r.begin(); i != r.end(); ++i)
                      {
                        switch (fileutils::detectFileType(fastaFiles[i]))
                        {
                        case fileutils::FastaFileType:
                          genomes[i] = parseGenome(fastaFiles[i], minContigSize);
                          break;

                        case fileutils::Gff3FileType:
                          genomes[i] = parseAnnotatedGenome(fastaFiles[i], colours, minContigSize);
                          break;

                        case fileutils::UnknownFileType:
                          LOG(ERROR) << "could not detect sequence file type: " << fastaFiles[i];
                          break;

                        default:
                          LOG(ERROR) << "unsuported file type: " << fastaFiles[i];
                          break;
                        }

                        pbar.tick();
                      }
                    });
  return genomes;
}

Dna4GenomeVector loadGenomes(const std::vector<fs::path> &fastaFiles,
                             size_t minContigSize)
{
  Colours colours{};
  return loadGenomes(fastaFiles, colours, minContigSize);
}

RandomGenomeIds::RandomGenomeIds(uint64_t maxId, size_t seed)
    : _maxId(maxId), _seed(seed), gen(_seed)
{
}

void RandomGenomeIds::generate(std::unordered_set<uint64_t> &out, size_t n)
{
  out.reserve(n);
  double sparsity = static_cast<double>(n) / static_cast<double>(_maxId);
  if (sparsity < _sparseLimit)
    generateSparse(out, n);
  else
    generateDense(out, n);
}

// fill values by random selection
void RandomGenomeIds::generateSparse(std::unordered_set<uint64_t> &out, size_t n)
{
  std::uniform_int_distribution<uint64_t> ids(1, _maxId);
  while (out.size() < n)
    out.insert(ids(gen));
}

// choose values to ignore by random selection, and then fill sequentially with complement
void RandomGenomeIds::generateDense(std::unordered_set<uint64_t> &out, size_t n)
{
  size_t np = _maxId /* =max num. IDs */ - n;
  std::unordered_set<uint64_t> ignore;
  ignore.reserve(np);

  std::uniform_int_distribution<uint64_t> ids(1, _maxId);
  while (ignore.size() < np)
    ignore.insert(ids(gen));

  for (size_t i = 1; i <= _maxId; ++i)
    if (!ignore.contains(i))
      out.insert(i);
}

std::vector<std::unordered_set<uint64_t>>
generateRandomGenomeIdSets(std::vector<size_t> const &setSizes, size_t reps, uint64_t maxId, size_t seed)
{
  RandomGenomeIds Gen(maxId, seed);
  std::vector<std::unordered_set<uint64_t>> buffers(setSizes.size() * reps);

  auto out = buffers.begin();
  for (size_t const &n : setSizes)
  {
    for (size_t r = 0; r < reps; ++r)
    {
      Gen.generate(*out, n);
      ++out;
    }
  }
  return buffers;
}

size_t SequenceRange::numKmers(uint8_t k) const noexcept
{
  size_t kmers = numInternalKmers(k);
  if (endIsTerminal())
    ++kmers;
  return kmers;
}

size_t SequenceRange::numTerminals(size_t k) const noexcept
{
  if (_begin == 0)
  {
    assert(_end >= k);
    return k;
  }
  else
    return 0;
}

size_t SequenceRange::numEdges(size_t k) const noexcept
{
  return numKmers(k) + numTerminals(k);
}

ChunkedDna4Genome::ChunkedDna4Genome(Dna4Genome &&genome, size_t granularity, size_t overlap)
    : _genome(std::move(genome)), _chunks()
{
  size_t ntChunkSize = std::max(_genome.medianContigSize() / granularity,
                                overlap);

  // chunk contigs

  for (Dna4Contig const &ctg : _genome.contigs)
  {
    for (SequenceVector const *drn : {&ctg.sequences.first, &ctg.sequences.second})
    {
      for (AnnotatedSequence const &rec : *drn)
      {
        assert(rec.annotations.size() == 0 &&
               "chunking scheme can only be applied to unannotated sequences");
        assert(rec.sequence.size() >= overlap &&
               "sequence size is less than requested overlap");
        size_t
            i = overlap,
            j;
        while (i < rec.sequence.size())
        {
          j = std::min(rec.sequence.size(), i + ntChunkSize);
          _chunks.emplace_back(std::ref(rec.sequence), i - overlap, j);
          i = j;
        }
        _chunks.back().setEndToTerminal();
      }
    }
  }
}

void ChunkedDna4Genome::reset_memory()
{
  _genome.reset_memory();
  std::vector<SequenceRange>{}.swap(_chunks);
}
