#include <maki/fasta.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <numeric>
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
    for (auto const &seq : v) resultSize += seq.size();
    resultSize += strlen(delim) * (v.size() - 1);

    std::string result;
    result.reserve(resultSize);

    for (auto const & seq : v)
    {
        if (!result.empty()) result.append(delim);
        for (auto c : seq)
            result.push_back(seqan3::to_char(c));
    }

    return result;
}

std::ostream &operator<<(std::ostream &os, const Gff3Token &obj)
{
    return os
           << "Gff3Record("
           << obj.beginPos << ", "
           << obj.endPos << ", "
           << obj.featureId << ")";
}

Gff3Record::Gff3Record(std::string &tmpData)
{
    // ensure integrity
    size_t numEntries = detail::countTsvWords(tmpData);
    if (numEntries != 9) { LOG(ERROR) << "invalid GFF3 line (only found " << numEntries << " entries)" << "\n\t" << tmpData; return; }
    size_t l = 0, r;
    // reference name
    r = tmpData.find('\t', l);
    if (r > l) { ref = tmpData.substr(l, r - l); }
    l = r + 1;
    // source free text descriptor and type of the feature (not recorded)
    l = tmpData.find('\t', tmpData.find('\t', l) + 1) + 1;
    // begin position of the interval
    r = tmpData.find('\t', l);
    if (r > l) { beginPos = std::stoul(tmpData.substr(l, r - l)) - 1 /* convert from 1-based to 0-based */; }
    l = r + 1;
    // end position of the interval
    // Note: this value is interpreted as the RHS of an inclusive interval. The rest of the
    // code works in half-open intervals, so this value is left unmodified.
    r = tmpData.find('\t', l);
    if (r > l) { endPos = std::stoul(tmpData.substr(l, r - l)) /* - 1 + 1 */; }
    l = r + 1;
    // score of the annotation
    l = tmpData.find('\t', l) + 1;
    // the strand
    r = tmpData.find('\t', l);
    if (r > l) { assert(r - l == 1); strand = tmpData[l]; }
    l = r + 1;
}

Gff3Record::Gff3Record(std::string_view ref_, size_t beginPos_, size_t endPos_, char strand_)
    : ref(ref_)
    , beginPos(beginPos_)
    , endPos(endPos_)
    , strand(strand_)
{}

Gff3Record::Gff3Record(std::string &&ref_, size_t beginPos_, size_t endPos_, char strand_)
    : ref(std::move(ref_))
    , beginPos(beginPos_)
    , endPos(endPos_)
    , strand(strand_)
{}

void Gff3Record::writeHeader(std::ostream &os)
{ os << "annotID\taccession\tbegin\tend\tstrand\n"; }

bool operator==(const Gff3Record &lhs, const Gff3Record &rhs) {
    return (lhs.ref == rhs.ref) &&
           (lhs.beginPos == rhs.beginPos) &&
           (lhs.endPos == rhs.endPos) &&
           (lhs.strand == rhs.strand) &&
           (lhs.featureId == rhs.featureId);
}

bool operator!=(const Gff3Record &lhs, const Gff3Record &rhs) {
    return (lhs.ref != rhs.ref) ||
           (lhs.beginPos != rhs.beginPos) ||
           (lhs.endPos != rhs.endPos) ||
           (lhs.strand != rhs.strand) ||
           (lhs.featureId != rhs.featureId);
}

std::ostream &operator<<(std::ostream &os, const Gff3Record &obj)
{
    return os
           << obj.featureId << '\t'
           << obj.ref << '\t'
           << obj.beginPos << '\t'
           << obj.endPos << '\t'
           << obj.strand;
}

void GenomeAnnotationList::sortContigStrand()
{
    std::sort(records.begin(), records.end(), ContigStrandLess{});
}

void GenomeAnnotationList::sort()
{
    std::sort(records.begin(), records.end(), AnnotLess{});
}

void GenomeAnnotationList::assignIds(FeatureIdAssigner<GenomeStatsSummary, Gff3Record> encoder)
{
    for (auto &record : records)
    {
        record.featureId = encoder.nextId();
    }
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
    Gff3Record prev = *it++;

    while (it != end)
    {
        Gff3Record curr = *it++;
        assert(prev.endPos <= curr.beginPos);
        std::swap(prev, curr);
    }
}

AnnotRange GenomeAnnotationList::findContig(const std::string &id)
{
    return {
        std::lower_bound(records.begin(), records.end(), id, ContigIdLess{}),
        std::upper_bound(records.begin(), records.end(), id, ContigIdLess{})};
}

bool ContigStrandLess::operator()(const Gff3Record &lhs, const Gff3Record &rhs) const
{
    if (lhs.ref != rhs.ref) return lhs.ref < rhs.ref;
    else return lhs.strand < rhs.strand;
}

bool AnnotLess::operator()(const Gff3Record &lhs, const Gff3Record &rhs) const
{
    if (lhs.ref != rhs.ref) return lhs.ref < rhs.ref;
    else if (lhs.strand != rhs.strand) return lhs.strand < rhs.strand;
    else return lhs.beginPos < rhs.beginPos;
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

void parseAnnotationStream(GenomeAnnotationList &annots, zstr::ifstream &gffStream)
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
            annots.records.emplace_back(tmpData);
        }
    }
    assert(tmpData == "##FASTA");
}

AnnotatedSequenceEdgeIterator::AnnotatedSequenceEdgeIterator(
    Dna4SequenceConstIter seqIter,
    size_t startPos,
    uint8_t windowSize,
    uint64_t nullFeatureId,
    AnnotationTokenVectorConstIter annotIter,
    AnnotationTokenVectorConstIter annotEnd)
    : inputSeqIter(seqIter + startPos),
      inputSeqPos(startPos),
      k(windowSize),
      nullFeatureId(nullFeatureId),
      inputAnnotIter(annotIter),
      inputAnnotEnd(annotEnd),
      currentAnnots(),
      numCurrentAnnots(0)
{
    if (inputSeqPos >= k)
    {
        while ((inputAnnotIter != inputAnnotEnd) && (inputAnnotIter->beginPos <= inputSeqPos - k))
        {
            currentAnnots.push_back(*inputAnnotIter);
            ++numCurrentAnnots;
            ++inputAnnotIter;
        }

        size_t numRemoved = std::erase_if(currentAnnots, [&inputSeqPos = inputSeqPos](const Gff3Token &rec)
                                          { return rec.endPos <= inputSeqPos; });
        numCurrentAnnots -= numRemoved;
    }
}

void AnnotatedSequenceEdgeIterator::operator++()
{
    // increment sequence position

    ++inputSeqIter;
    ++inputSeqPos;

    // remove any annotations that ended before this contig

    size_t numRemoved = std::erase_if(currentAnnots, [&inputSeqPos = inputSeqPos](const Gff3Token &rec)
                                      { return rec.endPos == inputSeqPos; });
    numCurrentAnnots -= numRemoved;

    // get new annotations that started within this region

    if (inputSeqPos >= k)
    {
        while ((inputAnnotIter != inputAnnotEnd) && (inputAnnotIter->beginPos == inputSeqPos - k))
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
        (currentSegmentEnd == inputSeqLength) && (currentSegmentFeatureId == nullFeatureId)};
}

void AnnotatedSequenceSegmentIterator::operator++()
{
    if (inputAnnotIter != inputAnnotEnd)
    {
        if (inputAnnotIter->beginPos + k > inputEdgePos)
        {
            currentSegmentBegin = inputEdgePos - k;
            inputEdgePos = inputAnnotIter->beginPos + k;
            currentSegmentEnd = inputEdgePos;
            currentSegmentFeatureId = /* UNANNOTATED */ nullFeatureId;
        }
        else
        {
            currentSegmentBegin = inputAnnotIter->beginPos;
            currentSegmentEnd = inputAnnotIter->endPos;
            currentSegmentFeatureId = inputAnnotIter->featureId;

            inputEdgePos = std::max(inputEdgePos, currentSegmentEnd);
            ++inputAnnotIter;
        }
    }
    else if (inputEdgePos <= inputSeqLength)
    {
        currentSegmentBegin = inputEdgePos - k;
        currentSegmentEnd = inputSeqLength;
        currentSegmentFeatureId = /* UNANNOTATED */ nullFeatureId;

        inputEdgePos = inputSeqLength + 1;
    }
    else
    {
        currentSegmentBegin = currentSegmentEnd = inputSeqLength;
    }
}

void parseFastaStream(Dna4Genome &genome, std::istream &fastaStream, uint64_t seqId, size_t minContigSize)
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

        // update genome stats

        genome.stats.minObsContigLength = genome.stats.minObsContigLength == 0 ? rec.sequence().size() : std::min(genome.stats.minObsContigLength, rec.sequence().size());
        genome.stats.maxObsContigLength = std::max(genome.stats.maxObsContigLength, rec.sequence().size());

        std::array<size_t, 4> ACGT = {0, 0, 0, 0};

        // insert forward and reverse sequences

        SequenceVector
            &fwdSeqs = out.sequences.first,
            &revSeqs = out.sequences.second;

        for (auto &&subSeq : rec.sequence() | std::views::split('N'_dna5) | myLengthFilter)
        {
            out.totalLength += 2 * subSeq.size();

            // insert forward

            auto &fwd = fwdSeqs.emplace_back(seqId, subSeq.size());

            for (auto &&c : subSeq | Dna5ToDna4)
            {
                ++ACGT[c.to_rank()];
                fwd.sequence.push_back(c);
            }

            // insert reverse

            auto &rev = revSeqs.emplace_back(seqId, subSeq.size());

            for (auto &&c : subSeq | std::views::reverse | seqan3::views::complement | Dna5ToDna4)
            {
                ++ACGT[c.to_rank()];
                rev.sequence.push_back(c);
            }
        }

        out.gcCount = ACGT[1] + ACGT[2];
    }
}

Dna4Genome parseGenome(fs::path fastaFile, FeatureIdAssigner<GenomeStatsSummary, Gff3Record> encoder, size_t minContigSize)
{
    // base input file stream that reads bytes

    zstr::ifstream zis(fastaFile);
    std::string fastaData(std::istreambuf_iterator<char>(zis), {});

    // process stream data

    Dna4Genome genome(fastaFile);

    std::istringstream fs(fastaData);
    parseFastaStream(genome, fs, encoder.baseId(), minContigSize);

    encoder.assignName(GenomeStatsSummary(genome));

    return genome;
}

Dna4Genome parseAnnotatedGenome(fs::path gff3File, FeatureIdAssigner<GenomeStatsSummary, Gff3Record> encoder, size_t minContigSize)
{
    // base input file stream that reads bytes

    zstr::ifstream zis(gff3File);

    // process stream data

    GenomeAnnotationList annots;

    parseAnnotationStream(annots, zis);
    annots.sortContigStrand();
    annots.assignIds(encoder);

    // read remaining data into string (and de-compress if necessary)

    std::string fastaData(std::istreambuf_iterator<char>(zis), {});

    // parse fasta records

    Dna4Genome genome(gff3File);
    genome.stats.numAnnotations = annots.size();

    std::istringstream fs(fastaData);
    seqan3::sequence_file_input seqInput(fs, seqan3::format_fasta{});

    for (auto &&rec : seqInput)
    {
        Dna4Contig &outContig = genome.contigs.emplace_back(rec.id().substr(0, rec.id().find(' ')), minContigSize);

        // update genome stats

        genome.stats.minObsContigLength = genome.stats.minObsContigLength == 0 ? rec.sequence().size() : std::min(genome.stats.minObsContigLength, rec.sequence().size());
        genome.stats.maxObsContigLength = std::max(genome.stats.maxObsContigLength, rec.sequence().size());

        AnnotRange contigAnnots = annots.findContig(outContig.contigName);
        outContig.numAnnotations = contigAnnots.size();

        outContig.insert(rec.sequence(), encoder.baseId(), contigAnnots, minContigSize + 1);
    }

    // store annotation metadata

    encoder.assignName(GenomeStatsSummary(genome));
    encoder.assignFeatures(std::move(annots.records));

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

void Dna4Contig::insert(const seqan3::dna5_vector &sequence, uint64_t seqFeatureId, AnnotRange annots, size_t minAnnotSize)
{
    size_t seqLength = sequence.size();
    std::vector<size_t> nPositions = findNs(sequence);

    // insert forward
    AnnotRange fwdAnnots = annots.getStrand('+');
    insertOrientation(sequences.first, sequence, seqFeatureId, nPositions, fwdAnnots, minAnnotSize);

    // get reverse complement view on sequence
    std::reverse(nPositions.begin(), nPositions.end());
    for (auto &c : nPositions) c = seqLength - c - 1;

    // translate reverse annotation positions to complement index
    AnnotRange revAnnots = annots.getStrand('-');
    for (auto it = revAnnots.begin; it != revAnnots.end; ++it) {
        Gff3Record &annot = *it;
        std::swap(annot.beginPos, annot.endPos);
        annot.beginPos = annot.beginPos > seqLength ? 0 : seqLength - annot.beginPos;
        annot.endPos = annot.endPos > seqLength ? 0 : seqLength - annot.endPos;
    }

    // insert reverse
    auto revCompSequence = sequence | std::views::reverse | seqan3::views::complement;
    insertOrientation(sequences.second, revCompSequence, seqFeatureId, nPositions, revAnnots, minAnnotSize);

    // revert annotation positions - important to maintain consistency with original annotation file
    for (auto it = revAnnots.begin; it != revAnnots.end; ++it) {
        Gff3Record &annot = *it;
        std::swap(annot.beginPos, annot.endPos);
        annot.beginPos = seqLength - annot.beginPos;
        annot.endPos = seqLength - annot.endPos;
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

size_t Dna4Contig::numKmers(uint8_t k) const
{
    auto accumulateKmers = [k](size_t total, const AnnotatedSequence &seq)
        { return total + seq.numKmers(k); };
    size_t total = std::accumulate(sequences.first.cbegin(), sequences.first.cend(), (size_t)0, accumulateKmers);
    return std::accumulate(sequences.second.cbegin(), sequences.second.cend(), total, accumulateKmers);
}

std::ostream &operator<<(std::ostream &os, const GenomeStats &obj)
{
    return os
           << obj.genomeName << '\t'
           << obj.sourceFile << '\t'
           << obj.numAnnotations << '\t'
           << obj.minObsContigLength << '\t'
           << obj.maxObsContigLength;
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
    return size_t{ k } * numFragments();
}

size_t Dna4Genome::gcCount() const
{
    auto accumulateGcContent = [](size_t total, const Dna4Contig &contig) -> size_t
    { return total + contig.gcCount; };
    return std::accumulate(contigs.cbegin(), contigs.cend(), (size_t)0, accumulateGcContent);
}

size_t Dna4Genome::medianContigSize() const
{
    std::vector<size_t> contigSizes;

    for (Dna4Contig const &ctg : contigs)
    {
        for (SequenceVector const *drn : { &ctg.sequences.first, &ctg.sequences.second })
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
        for (const SequenceVector *ornt : { &ctg.sequences.first, &ctg.sequences.second})
        {
            for (const AnnotatedSequence &rec : *ornt)
            {
                currentSize += ((2ul /*bits-per-nt Dna4Sequence*/ * rec.sequence.capacity()) + 7ul) / 8ul;
                currentSize += sizeof(Gff3Token) * rec.annotations.capacity();
            }
        }
    }
    return currentSize;
}

void Dna4Genome::reset_memory() {
    std::vector<Dna4Contig>{}.swap(contigs);
}

GenomeStatsSummary::GenomeStatsSummary(const Dna4Genome &genome)
    : GenomeStats(genome.stats)
    , sequenceLength(genome.length())
    , numContigs(genome.numContigs())
    , numFragments(genome.numFragments())
    , gcCount(genome.gcCount())
{
}

void GenomeStatsSummary::writeHeader(std::ostream &os) {
    os << "sequenceID\tgenomeName\tsourceFile\tnumAnnotations\tminContigLength\tmaxContigLength\tsequenceLength\tnumContigs\tnumFragments\tgcCount\n";
}

std::ostream &operator<<(std::ostream &os, const GenomeStatsSummary &obj) {
    return os
           << obj.genomeName << '\t'
           << obj.sourceFile << '\t'
           << obj.numAnnotations << '\t'
           << obj.minObsContigLength << '\t'
           << obj.maxObsContigLength << '\t'
           << obj.sequenceLength << '\t'
           << obj.numContigs << '\t'
           << obj.numFragments << '\t'
           << obj.gcCount;
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
        if (!fs::exists(file)) { ++missing; LOG(ERROR) << file << " not found"; }
    if (missing) { LOG(ERROR) << missing << " files not found"; throw "input files not found"; }
    else { LOG(INFO) << "parsed " << files.size() << " files"; }
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

uint64_t getMaxFeatureId(const Dna4GenomeVector &genomes)
{
    return tbb::parallel_reduce(
        tbb::blocked_range(genomes.cbegin(), genomes.cend()),
        (uint64_t)0,
        [](const tbb::blocked_range<Dna4GenomeVector::BaseGenomeVectorConstIter> &r, uint64_t maxColour) -> uint64_t
        {
            for (const Dna4Genome &genome : r)
            {
                for (const Dna4Contig &contig : genome.contigs)
                {
                    for (const SequenceVector *sequences : { &contig.sequences.first, &contig.sequences.second })
                    {
                        for (const AnnotatedSequence &sequence : *sequences)
                        {
                            maxColour = std::max(maxColour, sequence.nullFeatureId);

                            for (const Gff3Token &annot : sequence.annotations)
                            {
                                maxColour = std::max(maxColour, annot.featureId);
                            }
                        }
                    }
                }
            }
            return maxColour;
        },
        [](uint64_t lhs, uint64_t rhs) -> uint64_t
        { return std::max(lhs, rhs); });
}

Dna4GenomeVector loadGenomes(const std::vector<fs::path> &fastaFiles,
                             FeatureIdGenerator<GenomeStatsSummary, Gff3Record> &encoder,
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
                                  genomes[i] = parseGenome(fastaFiles[i], encoder.specify(i + 1, i), minContigSize);
                                  break;

                              case fileutils::Gff3FileType:
                                  genomes[i] = parseAnnotatedGenome(fastaFiles[i], encoder.specify(i + 1, i), minContigSize);
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

    // set colour width

    genomes.setMaxFeatureId(getMaxFeatureId(genomes));
    genomes.setBaseFeatureWidth(encoder.getBaseFeatureWidth());

    return genomes;
}

Dna4GenomeVector loadGenomes(const std::vector<fs::path> &fastaFiles,
                             size_t minContigSize)
{
    FeatureIdGenerator<GenomeStatsSummary, Gff3Record> encoder(fastaFiles.size());
    return loadGenomes(fastaFiles, encoder, minContigSize);
}

RandomGenomeIds::RandomGenomeIds(uint64_t maxId, size_t seed)
    : _maxId(maxId)
    , _seed(seed)
    , gen(_seed)
{}

void RandomGenomeIds::generate(std::unordered_set<uint64_t> &out, size_t n) {
    out.reserve(n);
    double sparsity = static_cast<double>(n) / static_cast<double>(_maxId);
    if (sparsity < _sparseLimit)
        generateSparse(out, n);
    else
        generateDense(out, n);
}

// fill values by random selection
void RandomGenomeIds::generateSparse(std::unordered_set<uint64_t> &out, size_t n) {
    std::uniform_int_distribution<uint64_t> ids(1, _maxId);
    while (out.size() < n)
        out.insert(ids(gen));
}

// choose values to ignore by random selection, and then fill sequentially with complement
void RandomGenomeIds::generateDense(std::unordered_set<uint64_t> &out, size_t n) {
    size_t np = _maxId/* =max num. IDs */ - n;
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
generateRandomGenomeIdSets(std::vector<size_t> const &setSizes, size_t reps, uint64_t maxId, size_t seed) {
    RandomGenomeIds Gen(maxId, seed);
    std::vector<std::unordered_set<uint64_t>> buffers(setSizes.size() * reps);
    
    auto out = buffers.begin();
    for (size_t const &n : setSizes) {
        for (size_t r = 0; r < reps; ++r) {
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
    : _genome(std::move(genome))
    , _chunks()
{
    size_t ntChunkSize = std::max(_genome.medianContigSize() / granularity,
                                  overlap);

    // chunk contigs

    for (Dna4Contig const &ctg : _genome.contigs)
    {
        for (SequenceVector const *drn : { &ctg.sequences.first, &ctg.sequences.second })
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

void ChunkedDna4Genome::reset_memory() {
    _genome.reset_memory();
    std::vector<SequenceRange>{}.swap(_chunks);
}

namespace fasta_utils
{

    Dna4Sequence randomSequence(size_t n, RandomDna4 &g)
    {
        Dna4Sequence seq;
        seq.reserve(n);

        for (size_t i = 0; i < n; ++i)
        {
            seq.push_back(g.next());
        }

        return seq;
    }

    Dna4Sequence randomSequence(size_t n, size_t seed)
    {
        RandomDna4 g(seed);
        return randomSequence(n, g);
    }

    Dna4Sequence randomSequence(size_t low, size_t high, size_t seed)
    {
        // sequence length generator

        std::mt19937 gen(seed);
        std::uniform_int_distribution<size_t> len(low, high);

        // random sequence generator

        RandomDna4 dna(seed);

        // generate sequence

        return randomSequence(len(gen), dna);
    }

    Dna4GenomeVector randomGenomes(size_t n, size_t seed, size_t low, size_t high)
    {
        Dna4GenomeVector colln(n);

        // sequence length generator

        std::mt19937 gen(seed);
        std::uniform_int_distribution<size_t> len(low, high);

        // random sequence generator

        RandomDna4 dna(seed);

        // generate sequences

        for (size_t i = 0; i < n; ++i)
        {
            Dna4Contig &ctg = colln[i].contigs.emplace_back();
            AnnotatedSequence const &rec = ctg.sequences.first.emplace_back(randomSequence(len(gen), dna), i + 1);
            ctg.totalLength += rec.sequence.size();
        }

        colln.setMaxFeatureId(getMaxFeatureId(colln));
        colln.setBaseFeatureWidth(ceil_log2(n + 1 /* null item */));

        return colln;
    }

    Dna4GenomeVector randomMutations(Dna4Sequence const &seq, size_t n, float rateLimit, size_t seed)
    {
        Dna4GenomeVector colln(n);

        // random generators

        std::mt19937 gen(seed);

        RandomDna4 dna(seed);

        size_t mutLimit = std::max((size_t)(rateLimit * (float)seq.size()), (size_t)1u);
        std::uniform_int_distribution<size_t> numMuts(1, mutLimit);
        std::uniform_int_distribution<size_t> mutPos(0, seq.size() - 1);

        // generate sequences

        for (size_t i = 0; i < n; ++i)
        {
            Dna4Contig &ctg = colln[i].contigs.emplace_back();
            AnnotatedSequence &rec = ctg.sequences.first.emplace_back(seq, i + 1);
            ctg.totalLength += rec.sequence.size();

            size_t
                z = 0,
                M = numMuts(gen);

            for (; z < M; ++z)
                rec.sequence[mutPos(gen)] = dna.next();
        }

        colln.setMaxFeatureId(getMaxFeatureId(colln));
        colln.setBaseFeatureWidth(ceil_log2(n + 1 /* null item */));

        return colln;
    }

} // namespace fasta_utils