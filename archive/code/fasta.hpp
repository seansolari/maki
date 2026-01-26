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
#include <maki/utils.hpp>
#include <maki/colourEncoder.hpp>

namespace fs = std::filesystem;

using colour::encoding::FeatureIdAssigner;
using colour::encoding::FeatureIdGenerator;

using namespace seqan3::literals;

struct Gff3Token {
    size_t beginPos, endPos;
    uint64_t featureId;
    bool operator==(const Gff3Token &) const =default;
    bool operator!=(const Gff3Token &) const =default;
    friend std::ostream &operator<<(std::ostream &, const Gff3Token &obj);
    friend void PrintTo(const Gff3Token &obj, std::ostream *os) { (*os) << obj; }
    size_t size() const noexcept { return endPos - beginPos; }
};

// based on https://seqan.readthedocs.io/en/main/Tutorial/InputOutput/GffAndGtfIO.html
struct Gff3Record {
    std::string ref;                // sequence name
    size_t beginPos;                // begin position of the interval
    size_t endPos;                  // end position of the interval
    char strand;                    // the strand
    uint64_t featureId = 0;         // (internal) feature ID
    Gff3Record(std::string &tmpData);
    Gff3Record(std::string_view ref_, size_t beginPos_, size_t endPos_, char strand_);
    Gff3Record(std::string &&ref_, size_t beginPos_, size_t endPos_, char strand_);
    static void writeHeader(std::ostream&);
    friend bool operator==(const Gff3Record &, const Gff3Record &);
    friend bool operator!=(const Gff3Record &, const Gff3Record &);
    friend std::ostream &operator<<(std::ostream &, const Gff3Record &obj);
    friend void PrintTo(const Gff3Record &obj, std::ostream *os) { (*os) << obj; return; }
    inline Gff3Token asToken() const noexcept { return {beginPos, endPos, featureId}; }
    operator Gff3Token() const noexcept { return asToken(); }
};

namespace std {
    template <>
    struct hash<Gff3Record> {
        size_t operator()(const Gff3Record& obj) const {
            size_t h_ = 0;
            h_ ^= hash<std::string>{}(obj.ref) + 0x9e3779b9;
            h_ ^= hash<size_t>{}(obj.beginPos) + 0x9e3779b9 + (h_ << 6) + (h_ >> 2);
            h_ ^= hash<size_t>{}(obj.endPos) + 0x9e3779b9 + (h_ << 6) + (h_ >> 2);
            h_ ^= hash<char>{}(obj.strand) + 0x9e3779b9 + (h_ << 6) + (h_ >> 2);
            return h_;
        }
    };

    template <>
    struct equal_to<Gff3Record> {
        size_t operator()(const Gff3Record& lhs, const Gff3Record& rhs) const {
            return (lhs.ref == rhs.ref) && (lhs.beginPos == rhs.beginPos)
                && (lhs.endPos == rhs.endPos) && (lhs.strand == rhs.strand);
        }
    };
} // namespace std

struct GenomeStats;
struct GenomeStatsSummary;

using RecordVector = std::vector<Gff3Record>;
using RecordVectorIter = RecordVector::iterator;
using RecordVectorConstIter = RecordVector::const_iterator;

using AnnotationTokenVector = std::vector<Gff3Token>;
using AnnotationTokenVectorIter = AnnotationTokenVector::iterator;
using AnnotationTokenVectorConstIter = AnnotationTokenVector::const_iterator;

struct AnnotRange {
    RecordVectorIter begin, end;
    void sortByStart();
    AnnotRange getStrand(char c) const;
    RecordVectorIter encloseStart(size_t z) const;
    size_t size() const noexcept { return end - begin; }
};

struct GenomeAnnotationList {
    RecordVector records;
    size_t size() const noexcept { return records.size(); }
    // sort by contig, then strand of features
    void sortContigStrand();
    // sort by contig, strand and start position of features
    void sort();
    // assign feature Ids to each annotation
    void assignIds(FeatureIdAssigner<GenomeStatsSummary, Gff3Record> encoder);
    // ensure that consecutive features are not overlapping
    void ensureDisjoint();
    AnnotRange findContig(const std::string &id);
protected:
    static void ensureDisjointRegion(RecordVectorConstIter it, RecordVectorConstIter end);
};

struct ContigIdLess {
    bool operator()(const Gff3Record &lhs, const Gff3Record &rhs) const { return lhs.ref < rhs.ref; }
    bool operator()(const Gff3Record &lhs, const std::string &rhs) const { return lhs.ref < rhs; }
    bool operator()(const std::string &lhs, const Gff3Record &rhs) const { return lhs < rhs.ref; }
};

struct StrandLess {
    bool operator()(const Gff3Record &lhs, const Gff3Record &rhs) const { return lhs.strand < rhs.strand; }
    bool operator()(const Gff3Record &lhs, char rhs) const { return lhs.strand < rhs; }
    bool operator()(char lhs, const Gff3Record &rhs) const { return lhs < rhs.strand; }
};

struct BeginPosLess {
    bool operator()(const Gff3Record &lhs, const Gff3Record &rhs) const { return lhs.beginPos < rhs.beginPos; }
    bool operator()(const Gff3Record &lhs, size_t rhs) const { return lhs.beginPos < rhs; }
    bool operator()(size_t lhs, const Gff3Record &rhs) const { return lhs < rhs.beginPos; }
};

struct ContigStrandLess { bool operator()(const Gff3Record &lhs, const Gff3Record &rhs) const; };

struct AnnotLess { bool operator()(const Gff3Record &lhs, const Gff3Record &rhs) const; };

using Dna4Sequence = seqan3::bitpacked_sequence<seqan3::dna4>;

std::string toString(Dna4Sequence const &);
std::string toString(std::vector<Dna4Sequence> const&, const char*);

struct RandomDna4 {
protected:
    std::array<seqan3::dna4, 4> _ALPHA = { 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4 };
    std::mt19937 gen;
    std::uniform_int_distribution<> distr;
public:
    template <class SeedSeq> RandomDna4(SeedSeq &&seed) : gen(std::forward<SeedSeq>(seed)), distr(0, 3) {}
    seqan3::dna4 next() { return _ALPHA[distr(gen)]; }
};

using Dna4SequenceConstIter = Dna4Sequence::const_iterator;
using AnnotationTokenList = std::vector<Gff3Token>;

class AnnotatedSequenceEdgeIterator {
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

struct AnnotatedSequenceSegment {
    Dna4SequenceConstIter seqBegin, seqEnd;
    uint64_t featureId;
    bool endIsTerminal;
    size_t size() const noexcept { return seqEnd - seqBegin; }
};

class AnnotatedSequenceSegmentSentinel {};

class AnnotatedSequenceSegmentIterator {
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
    inline size_t segmentNumInternalEdges() const noexcept { return currentSegmentEnd - currentSegmentBegin - k; }
    void operator++();
};

struct AnnotatedSequence {
    Dna4Sequence sequence;
    AnnotationTokenVector annotations;
    uint64_t nullFeatureId;
    AnnotatedSequence(Dna4Sequence const &seq) : sequence(seq), annotations(), nullFeatureId(1u) {}
    AnnotatedSequence(Dna4Sequence const &seq, uint64_t id) : sequence(seq), annotations(), nullFeatureId(id) {}
    AnnotatedSequence(Dna4Sequence &&seq) : sequence(std::move(seq)), annotations(), nullFeatureId(1u) {}
    AnnotatedSequence(Dna4Sequence &&seq, uint64_t id) : sequence(std::move(seq)), annotations(), nullFeatureId(id) {}
    AnnotatedSequence(uint64_t id, size_t length) : sequence(), annotations(), nullFeatureId(id){ sequence.reserve(length); }
    inline AnnotatedSequenceEdgeIterator colouredSegmentsEdgeBegin(uint8_t windowSize) const { return AnnotatedSequenceEdgeIterator(sequence.cbegin(), windowSize, nullFeatureId, annotations.cbegin(), annotations.cend()); }
    inline AnnotatedSequenceEdgeIterator colouredSegmentsEdgeBegin(size_t startPos, uint8_t windowSize) const { return AnnotatedSequenceEdgeIterator(sequence.cbegin(), static_cast<int32_t>(startPos), windowSize, nullFeatureId, annotations.cbegin(), annotations.cend()); }
    inline AnnotatedSequenceEdgeIterator colouredSegmentsEdgeEnd() const { return AnnotatedSequenceEdgeIterator(sequence.cend()); }
    inline AnnotatedSequenceSegmentIterator colouredSegmentsBegin(uint8_t windowSize) const { return AnnotatedSequenceSegmentIterator(sequence.cbegin(), sequence.size(), windowSize, nullFeatureId, annotations.cbegin(), annotations.cend()); }
    inline AnnotatedSequenceSegmentSentinel colouredSegmentsEnd() const { return AnnotatedSequenceSegmentSentinel(); }
    size_t numKmers(uint8_t k) const;
};

using SequenceVector = std::vector<AnnotatedSequence>;

struct Dna4Contig {
    std::pair</* FORWARD */ SequenceVector, /* REVERSE */ SequenceVector> sequences;
    std::string contigName;
    size_t minFragmentSize = 0, totalLength = 0, numAnnotations = 0, gcCount = 0; 
    Dna4Contig() = default;
    Dna4Contig(std::string contigName_, size_t minContigSize_) : contigName(contigName_), minFragmentSize(minContigSize_) {}
    void insert(const seqan3::dna5_vector &, uint64_t seqFeatureId, AnnotRange annots, size_t minAnnotSize);
    std::vector<size_t> findNs(const seqan3::dna5_vector &);
    
    template <typename R>
        requires std::ranges::input_range<R> && std::same_as<std::ranges::range_value_t<R>, seqan3::dna5>
    void insertOrientation(SequenceVector &, R &&seqView, uint64_t seqFeatureId, const std::vector<size_t> &nPos, AnnotRange annots, size_t minAnnotSize);

    template <typename Iter>
    void insertFragment(SequenceVector &vec, Iter seqIter, size_t startPos, size_t endPos, uint64_t seqFeatureId, AnnotRange &annots, AnnotationTokenList &currentAnnots, size_t minAnnotSize);

    // includes concluding terminal edge, but not starting terminals
    size_t numKmers(uint8_t k) const;
    inline size_t numFragments() const noexcept { return sequences.first.size() + sequences.second.size(); }
};

struct GenomeStats {
    std::string genomeName, sourceFile;
    size_t minObsContigLength = 0, maxObsContigLength = 0, numAnnotations = 0;
    GenomeStats() = default;
    GenomeStats(std::string_view genomeName_, std::string_view sourceFile_) : genomeName(genomeName_), sourceFile(sourceFile_) {}
    friend std::ostream &operator<<(std::ostream &, const GenomeStats &obj);
};

struct Dna4Genome {
    GenomeStats stats;
    std::vector<Dna4Contig> contigs;
    Dna4Genome() = default;
    Dna4Genome(std::string file) : stats(fileutils::extractSequenceName(file), file) {}
    size_t length() const;
    // includes concluding terminal edge, but not starting terminals
    size_t numKmers(uint8_t k) const;
    inline size_t numContigs() const noexcept { return contigs.size(); }
    size_t numFragments() const;
    // number of starting terminals, does not include trailing terminal edge
    size_t numTerminals(uint8_t k) const;
    size_t gcCount() const;
    size_t medianContigSize() const;
    // number of bytes used
    size_t rss() const;
    void reset_memory();
};

template <typename Fn>
concept CallsRange = requires(Fn f_, Dna4SequenceConstIter begin, Dna4SequenceConstIter end, bool b)
{ { f_(begin, end, b) } -> std::same_as<void>; };

template <typename Fn>
concept CallsColouredRange = requires(Fn f_, Dna4SequenceConstIter begin, Dna4SequenceConstIter end, uint64_t colour, bool b)
{ { f_(begin, end, colour, b) } -> std::same_as<void>; };

template <class ...Args>
struct UnwindGenome { UnwindGenome(Args ...args); };

template<>
struct UnwindGenome<uint8_t> {
    UnwindGenome(uint8_t windowSize) : windowSize(windowSize) {}

    template <typename Fn_>
        requires CallsRange<Fn_> || CallsColouredRange<Fn_>
    void forEach(Dna4Genome const &g, Fn_ f) const {
        for (const Dna4Contig &contig : g.contigs){
            for (const SequenceVector *orientation : { &contig.sequences.first, &contig.sequences.second }){
                for (const AnnotatedSequence &seqObj : *orientation){
                    for (AnnotatedSequenceSegmentIterator it = seqObj.colouredSegmentsBegin(windowSize); it != seqObj.colouredSegmentsEnd(); ++it) {
                        AnnotatedSequenceSegment segment = *it;
                        if constexpr ( CallsRange<Fn_> ) f(segment.seqBegin, segment.seqEnd, segment.endIsTerminal);
                        else f(segment.seqBegin, segment.seqEnd, segment.featureId, segment.endIsTerminal);
                    }
                }
            }
        }
    }

    template <typename Fn_>
    void forEachBegin(Dna4Genome const &g, Fn_ f) const {
        for (const Dna4Contig &contig : g.contigs) {
            for (const SequenceVector *orientation : { &contig.sequences.first, &contig.sequences.second }) {
                for (const AnnotatedSequence &seqObj : *orientation) {
                    for (AnnotatedSequenceSegmentIterator it = seqObj.colouredSegmentsBegin(windowSize); it != seqObj.colouredSegmentsEnd(); ++it) {
                        AnnotatedSequenceSegment segment = *it;
                        f(segment.seqBegin);
                    }
                }
            }
        }
    }

protected:
    uint8_t  windowSize;
};

template<>
struct UnwindGenome<> {
    UnwindGenome() {}

    template <typename Fn_>
        requires CallsRange<Fn_> || CallsColouredRange<Fn_>
    void forEach(Dna4Genome const &g, Fn_ f) const {
        for (const Dna4Contig &contig : g.contigs) {
            for (const SequenceVector *orientation : { &contig.sequences.first, &contig.sequences.second }) {
                for (const AnnotatedSequence &seqObj : *orientation) {
                    if constexpr ( CallsRange<Fn_> ) f(seqObj.sequence.cbegin(), seqObj.sequence.cend(), true);
                    else f(seqObj.sequence.cbegin(), seqObj.sequence.cend(), seqObj.nullFeatureId, true);
                }
            }
        }
    }

    template <typename Fn_>
    void forEachBegin(Dna4Genome const &g, Fn_ f) const {
        for (const Dna4Contig &contig : g.contigs) {
            for (const SequenceVector *orientation : { &contig.sequences.first, &contig.sequences.second }) {
                for (const AnnotatedSequence &seqObj : *orientation)
                    f(seqObj.sequence.cbegin());
            }
        }
    }
};

template <typename T>
struct UnwindTraits {};

template <typename T>
using range_t = typename UnwindTraits<T>::range_type;

template <class ...Args>
struct UnwindTraits<UnwindGenome<Args...>> {
    typedef Dna4Genome range_type;
};

struct GenomeStatsSummary : public GenomeStats {
    size_t sequenceLength, numContigs, numFragments, gcCount;
    GenomeStatsSummary() = default;
    GenomeStatsSummary(const Dna4Genome &genome);
    static void writeHeader(std::ostream&);
    friend std::ostream &operator<<(std::ostream &, const GenomeStatsSummary &obj);
};

namespace parsing
{
    template <typename O>
    inline O dna4ToRank(seqan3::dna4 &&nt)
    { return static_cast<O>(seqan3::to_rank(std::forward<seqan3::dna4>(nt))); }

    inline uint8_t dna4ToShort(seqan3::dna4 &&nt)
    { return dna4ToRank<uint8_t>(std::forward<seqan3::dna4>(nt)); }

    inline uint64_t dna4ToLong(seqan3::dna4 &&nt)
    { return dna4ToRank<uint64_t>(std::forward<seqan3::dna4>(nt)); }

    inline uint64_t dna4ToDna5(uint64_t &&nt)
    { return (nt + 1) & 0b00000111; }

    inline seqan3::dna4 dna5ToDna4(uint64_t &&rnk)
    { return seqan3::assign_rank_to((rnk - 1) & 0b11, seqan3::dna4{}); }

} // namespace parsing

void parseAnnotationStream(GenomeAnnotationList &annots, zstr::ifstream &gffStream);
void parseFastaStream(Dna4Genome &genome, std::istream &fastaStream, uint64_t seqId, size_t minContigSize);
Dna4Genome parseGenome(fs::path fastaFile, FeatureIdAssigner<GenomeStatsSummary, Gff3Record> encoder, size_t minContigSize);
Dna4Genome parseAnnotatedGenome(fs::path gffFile, FeatureIdAssigner<GenomeStatsSummary, Gff3Record> encoder, size_t minContigSize);
std::vector<fs::path> readFilePaths(const char *manifest_file);

class Dna4GenomeVector {
public:
    using BaseGenomeVector = std::vector<Dna4Genome>;
    using BaseGenomeVectorIter = BaseGenomeVector::iterator;
    using BaseGenomeVectorConstIter = BaseGenomeVector::const_iterator;
protected:
    BaseGenomeVector genomes;
public:
    uint64_t maxFeatureId;
    uint8_t baseFeatureWidth; // number of LSB that correspond to genome ID
    Dna4GenomeVector(size_t n) : genomes(n) {}
    operator BaseGenomeVector const&() const noexcept { return genomes; }
    BaseGenomeVector const &buffer() const noexcept { return genomes; }
    // sends current data to temporary object, whose deconstructor will then be called
    void reset_memory() { BaseGenomeVector().swap(genomes); }
    size_t numGenomes() const { return genomes.size(); }
    void setMaxFeatureId(uint64_t maxId_) { maxFeatureId = maxId_; }
    void setBaseFeatureWidth(uint8_t w_) { baseFeatureWidth = w_; }
    uint64_t getMaxFeatureId() const { return maxFeatureId; }
    uint8_t getFeatureWidth() const { return ceil_log2(maxFeatureId); }
    uint8_t getBaseFeatureWidth() const noexcept { return baseFeatureWidth; }
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

uint64_t getMaxFeatureId(const Dna4GenomeVector &);
Dna4GenomeVector loadGenomes(const std::vector<fs::path> &fastaFiles, FeatureIdGenerator<GenomeStatsSummary, Gff3Record> &encoder, size_t minContigSize);
Dna4GenomeVector loadGenomes(const std::vector<fs::path> &fastaFiles, size_t minContigSize);

class RandomGenomeIds {
    constexpr static double _sparseLimit = 0.5;
public:
    RandomGenomeIds(uint64_t maxId, size_t seed);
    void generate(std::unordered_set<uint64_t> &out, size_t n);
protected:
    uint64_t     _maxId;
    size_t       _seed;
    std::mt19937 gen;
    void generateSparse(std::unordered_set<uint64_t> &out, size_t n);
    void generateDense(std::unordered_set<uint64_t> &out, size_t n);
};

std::vector<std::unordered_set<uint64_t>>
generateRandomGenomeIdSets(std::vector<size_t> const &setSizes, size_t reps, uint64_t maxId, size_t seed);

struct SequenceRange {
    SequenceRange(Dna4Sequence const &seq_, size_t begin_, size_t end_) : _seq(seq_) , _begin(begin_) , _end(end_) , _endIsTerminal(false) {}
    Dna4SequenceConstIter begin() const { return _seq.cbegin() + _begin; }
    Dna4SequenceConstIter end() const { return _seq.cbegin() + _end; }
    size_t size() const noexcept { return _end - _begin; }
    size_t numInternalKmers(uint8_t k) const noexcept { assert(size() >= k); return size() - k; }
    size_t numKmers(uint8_t k) const noexcept;
    size_t numTerminals(/* size_t for cmp */size_t k) const noexcept;
    size_t numEdges(size_t k) const noexcept;
    bool endIsTerminal() const noexcept { return _endIsTerminal; }
    void setEndToTerminal() { _endIsTerminal = true; }
    friend struct UnwindChunkedGenome;
protected:
    Dna4Sequence const &_seq;
    size_t _begin, _end;
    bool _endIsTerminal;
};

class ChunkedDna4Genome {
public:
    typedef std::vector<SequenceRange>::const_iterator iterator;
    ChunkedDna4Genome(Dna4Genome &&genome, size_t granularity, size_t overlap);
    auto begin() const { return _chunks.begin(); }
    auto end() const { return _chunks.end(); }
    operator Dna4Genome const&() const noexcept { return _genome; }
    operator std::vector<SequenceRange>&() noexcept { return _chunks; }
    operator std::vector<SequenceRange> const&() const noexcept { return _chunks; }
    SequenceRange const& operator[](size_t i) const { return _chunks[i]; }
    size_t chunks() const noexcept { return _chunks.size(); }
    size_t numEdges() const noexcept { return _genome.length() + _genome.numFragments(); }
    size_t sequenceRSS() const { return _genome.rss(); }
    size_t chunksRSS() const { return sizeof(SequenceRange) * _chunks.capacity(); }
    void reset_memory();
protected:
    Dna4Genome _genome;
    std::vector<SequenceRange> _chunks;
};

struct UnwindChunkedGenome {
    template <typename Fn_> requires CallsRange<Fn_> void forEach(SequenceRange const &rng, Fn_ f) const { f(rng.begin(), rng.end(), rng.endIsTerminal()); }
    template <typename Fn_> void forEachBegin(SequenceRange const &rng, Fn_ f) const { if (rng._begin == 0) f(rng.begin()); }
};

template <> struct UnwindTraits<UnwindChunkedGenome> { typedef SequenceRange range_type; };

namespace fasta_utils {

    Dna4Sequence randomSequence(size_t n, RandomDna4 &g);
    Dna4Sequence randomSequence(size_t n, size_t seed);
    Dna4Sequence randomSequence(size_t low, size_t high, size_t seed);
    Dna4GenomeVector randomGenomes(size_t n, size_t seed, size_t low, size_t high);
    Dna4GenomeVector randomMutations(Dna4Sequence const &seq, size_t n, float rateLimit, size_t seed);

} // namespace fasta_utils

// template definitions

template <typename R>
    requires std::ranges::input_range<R> && std::same_as<std::ranges::range_value_t<R>, seqan3::dna5>
void Dna4Contig::insertOrientation(SequenceVector &vec, R &&seqView, uint64_t seqFeatureId, const std::vector<size_t> &nPositions, AnnotRange annots, size_t minAnnotSize) {
    // sort annotations by start position
    annots.sortByStart();
    // insert fragments
    AnnotationTokenList currentAnnots;
    size_t seqPos = 0;
    std::ranges::iterator_t<R> seqIter = seqView.begin();
    for (size_t nPos : nPositions) {
        size_t nextContigSize = nPos - seqPos;
        if (nextContigSize >= minFragmentSize)
            insertFragment(vec, seqIter, seqPos, nPos, seqFeatureId, annots, currentAnnots, minAnnotSize);
        seqPos += nextContigSize + 1;
    }
    if (seqView.size() - seqPos >= minFragmentSize)
        insertFragment(vec, seqIter, seqPos, seqView.size(), seqFeatureId, annots, currentAnnots, minAnnotSize);
}

template <typename Iter>
void Dna4Contig::insertFragment(SequenceVector &vec, Iter seqIter, size_t startPos, size_t endPos, uint64_t seqFeatureId, AnnotRange &annots, AnnotationTokenList &currentAnnots, size_t minAnnotSize) {
    size_t fragmentSize = endPos - startPos;
    totalLength += fragmentSize;
    std::array<size_t, 4> ACGT = {0, 0, 0, 0};
    AnnotatedSequence &seqObj = vec.emplace_back(seqFeatureId, fragmentSize);
    for (auto it = seqIter + startPos; it != seqIter + endPos; ++it) {
        auto c = static_cast<seqan3::dna4>(*it);
        ++ACGT[c.to_rank()];
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
    gcCount += ACGT[1] + ACGT[2];
}