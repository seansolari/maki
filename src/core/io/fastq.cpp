#include <maki/reads.hpp>

#include <algorithm>
#include <numeric>
#include <ranges>
#include <oneapi/tbb.h>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/utility/views/zip.hpp>
#include <zstr.hpp>

using namespace seqan3::literals;

namespace reads
{

    std::ostream &operator<<(std::ostream &os, const DataFilePair &obj)
    {
        return os
               << "PairedReads("
               << obj.forwardFile << ", "
               << obj.reverseFile
               << ")";
    }

    bool DataFilePair::gzCompressed() const {
        bool
            lCx = forwardFile.ends_with(".gz"),
            rCx = reverseFile.ends_with(".gz");
        
        if (lCx != rCx)
            throw std::runtime_error("paired files " + forwardFile + " and " + reverseFile + " have different compression");
        else
            return lCx;
    }

    std::vector<DataFilePair> pairInputFiles(const std::vector<fs::path> &files)
    {
        // check even number of files

        assert((files.size() & 1u) == 0u);
        size_t numPairs = files.size() >> 1;

        // sort files lexicographically

        std::vector<size_t> fileIndex(files.size());
        std::iota(fileIndex.begin(), fileIndex.end(), 0);
        std::sort(fileIndex.begin(), fileIndex.end(), [&](size_t i, size_t j)
                  { return files[i] < files[j]; });

        // pair and check files

        std::vector<DataFilePair> pairs;
        pairs.reserve(numPairs);

        for (size_t i = 0; i < numPairs; ++i)
        {
            std::string
                fwd = files[fileIndex[2 * i]].string(),
                rev = files[fileIndex[(2 * i) + 1]].string();

            std::string_view
                fwdName = fileutils::extractSequenceName(fwd),
                revName = fileutils::extractSequenceName(rev);

            if (fwdName.back() != '1')
            {
                throw std::runtime_error("Paired files " + fwd + " and " + rev + ", but " + fwd + " does not end in 1.");
            }
            else if (revName.back() != '2')
            {
                throw std::runtime_error("Paired files " + fwd + " and " + rev + ", but " + rev + " does not end in 2.");
            }
            else if (fwdName.substr(0, fwdName.size() - 1) != revName.substr(0, revName.size() - 1))
            {
                throw std::runtime_error("Paired files " + fwd + " and " + rev + ", but prefixes do not match.");
            }
            else
            {
                pairs.emplace_back(std::move(fwd), std::move(rev));
            }
        }

        return pairs;
    }

    void
    detail::CharBuffer::flushTo(detail::CharBuffer &dest)
    {
        assert(_size < dest.bytesRemaining());
        std::copy(_data, _data + _size, dest.back());
        dest._size += _size;
        _size = 0;
    }

    void
    detail::CharBuffer::swap(detail::CharBuffer &rhs)
    {
        std::swap(_data, rhs._data);
        std::swap(_capacity, rhs._capacity);
        std::swap(_size, rhs._size);
    }

    const char *
    detail::CharBuffer::prev(char c) const
    {
        auto
            _rend = std::make_reverse_iterator(begin()),
            _it = std::find(std::make_reverse_iterator(back()), _rend, c);
        if (_it == _rend)
            return back();
        else
            return begin() + (_rend - _it) - 1;
    }

    const char *
    detail::CharBuffer::nPrevFrom(const char *end_, size_t n, char c) const
    {
        auto
            _rend = std::make_reverse_iterator(begin()),
            _it = std::find_if(std::make_reverse_iterator(end_),
                               _rend,
                               [&](char const &z_) { return (z_ == c) && !(--n); });
        if (_it == _rend)
            return back();
        else
            return begin() + (_rend - _it) - 1;
    }

    void
    detail::CharBuffer::write(const char *src, size_t n)
    {
        assert(_size + n < _capacity);
        std::copy_n(src, n, back());
        _size += n;
    }

    void
    detail::CharBuffer::read(std::istream &is)
    {
        is.read(back(), bytesRemaining());
        _size += is.gcount();
    }

    void
    detail::CharBuffer::flushOverflow(const char *end_, detail::CharBuffer &out)
    {
        size_t newsize_ = end_ - _data;
        assert(newsize_ < _size);
        out.write(end_, _size - newsize_);
        _size = newsize_;
    }

    std::ostream &detail::operator<<(std::ostream &os, const detail::CharBuffer &obj)
    {
        os << "CharBuffer[\n";
        for (size_t i = 0; i < obj.size(); ++i)
            os << obj[i];
        os << "]";
        return os;
    }

    detail::SharedBufferPair *
    detail::Operator_ReadChunk::operator()(tbb::detail::d1::flow_control &fc) const
    {
        if (ifwd->eof() && irev->eof())
        {
            if (!ovfl->empty())
                throw std::runtime_error("paired-end structure is compromised");

            fc.stop();
            return nullptr;
        }
        else return readChunk();
    }

    detail::SharedBufferPair *
    detail::Operator_ReadChunk::readChunk() const
    {
        // find next available buffer

        for (SharedBufferPair *bp = buffers; bp < buffers + numBuffers; ++bp)
        {
            if (!bp->avail())
                continue;
            bp->setUnavail();

            // flush overflow data

            if (!ovfl->empty())
                bp->swap(*ovfl);

            // fill chars left

            bp->fwd.read(*ifwd);
            bp->rev.read(*irev);

            // truncate sequences and splash to overflow

            const char
                *fend = completedRec(ifwd) ? bp->fwd.back() : bp->fwd.prev(DELIM),
                *rend = completedRec(irev) ? bp->rev.back() : bp->rev.prev(DELIM);

            if ((fend == bp->fwd.begin()) || (rend == bp->rev.begin()))
                throw std::runtime_error("buffer size not begin enough to parse read record");

            size_t
                fc = bp->fwd.countTo(fend, DELIM),
                rc = bp->rev.countTo(rend, DELIM),
                recs = std::min(fc, rc);

            if (fc < recs)
                fend = bp->fwd.nPrevFrom(fend, recs - fc, DELIM);
            else if (rc < recs)
                rend = bp->rev.nPrevFrom(rend, recs - rc, DELIM);

            bp->fwd.flushOverflow(fend, ovfl->fwd);
            bp->rev.flushOverflow(rend, ovfl->rev);

            return bp;
        }

        throw std::runtime_error("no buffers available");
    }

    void
    Read::push(seqan3::dna5_vector &&data, size_t minFragLength)
    {
        auto lengthFilter = std::views::filter(
            [&](auto &&seq)
            {
                return std::ranges::size(std::forward<decltype(seq)>(seq)) >= minFragLength;
            });

        auto Dna5ToDna4 = std::views::transform(
            [](auto &&c) { // https://docs.seqan.de/seqan3/main_user/cookbook.html#cookbook_convert_alphabet_range
                return static_cast<seqan3::dna4>(std::forward<decltype(c)>(c));
            });

        for (auto &&subSeq :
             data | std::views::split('N'_dna5) | lengthFilter)
        {
            Dna4Sequence
                &f = raw.emplace_back(),
                &r = rcomp.emplace_back();

            f.reserve(subSeq.size());
            for (auto &&c : subSeq | Dna5ToDna4)
                f.push_back(c);
            
            r.reserve(subSeq.size());
            for (auto &&c : subSeq | std::views::reverse | seqan3::views::complement | Dna5ToDna4)
                r.push_back(c);
        }
    }

    size_t Read::numFragments() const noexcept
    {
        return raw.size() + rcomp.size();
    }

    size_t Read::rss() const
    {
        size_t tot = (id.capacity() * sizeof(std::string::value_type))
            + (raw.capacity() * sizeof(Dna4Sequence))
            + (rcomp.capacity() * sizeof(Dna4Sequence));
        for (Dna4Sequence const &seq : raw)
            tot += ((seq.capacity() + 63) / 64) * 8;
        for (Dna4Sequence const &seq : rcomp)
            tot += ((seq.capacity() + 63) / 64) * 8;
        return tot;
    }

    size_t ReadChunks::numFragments() const
    {
        auto accumulateFragments = [](size_t total, const ReadVector &reads) -> size_t
        {
            for (const Read &r : reads)
                total += r.numFragments();
            return total;
        };
        return std::accumulate(_data.cbegin(), _data.cend(), (size_t)0, accumulateFragments);
    }

    size_t ReadChunks::numTerminals(uint8_t k) const
    {
        return size_t{ k } * numFragments();
    }

    size_t ReadChunks::rss() const
    {
        size_t tot = 0;
        for (auto const &ch : _data)
            for (Read const &r : ch)
                tot += r.rss();
        return tot;
    }

    void
    parse::Operator_ParseChunk::operator()(detail::SharedBufferPair *in) const
    {
        ReadVector &out = buffers->local();

        // parse paired data

        detail::basic_spanbuf
            fwdView = in->viewForward(),
            revView = in->viewReverse();

        std::istream
            fwdStream(&fwdView),
            revStream(&revView);

        seqan3::sequence_file_input
            ifwd{fwdStream, seqan3::format_fastq{}},
            irev{revStream, seqan3::format_fastq{}};

        for (auto &&[rec1, rec2] : seqan3::views::zip(ifwd, irev))
        {
            if ((rec1.sequence().size() >= minReadLength) &&
                (rec2.sequence().size() >= minReadLength))
            {
                Read &v = out.emplace_back(std::move(rec1.id()));
                v.push(std::move(rec1).sequence(), minFragLength);
                v.push(std::move(rec2).sequence(), minFragLength);
            }
        }

        // free data

        in->clear();
        in->setAvail();
    }

    ReadChunks
    parse::parsePairedFastq(const DataFilePair &fp, size_t minReadLength, size_t minFragLength, size_t threads, size_t buffersize)
    {
        if (fp.gzCompressed())
        {
            zstr::ifstream
                ifwd(fp.forwardFile),
                irev(fp.reverseFile);
            return parsePairedFastq(&ifwd, &irev, minReadLength, minFragLength, threads, buffersize);
        }
        else
        {
            std::basic_ifstream<char>
                ifwd(fp.forwardFile),
                irev(fp.reverseFile);
            return parsePairedFastq(&ifwd, &irev, minReadLength, minFragLength, threads, buffersize);
        }
    }

    ReadChunks
    parse::parsePairedFastq(std::basic_istream<char> *ifwd, std::basic_istream<char> *irev, size_t minReadLength, size_t minFragLength, size_t threads, size_t buffersize)
    {
        detail::SharedBufferPairVector buffers(threads, buffersize);
        detail::BufferPair ovflw(buffersize);

        tbb::enumerable_thread_specific<ReadVector> threadOutput;

        detail::Operator_ReadChunk ReadChunk(ifwd, irev,
                                             '@',
                                             buffers.data(), buffers.numBuffers(),
                                             &ovflw);
        Operator_ParseChunk ParseChunk(&threadOutput, minReadLength, minFragLength);

        oneapi::tbb::parallel_pipeline(threads,
                                       oneapi::tbb::make_filter<void, detail::SharedBufferPair *>(serial_in_order, ReadChunk),
                                       oneapi::tbb::make_filter<detail::SharedBufferPair *, void>(parallel, ParseChunk));

        return ReadChunks(std::make_move_iterator(threadOutput.begin()),
                          std::make_move_iterator(threadOutput.end()));
    }

    size_t detail::countBps(const ReadVector &data)
    {
        return std::accumulate(
            data.cbegin(),
            data.cend(),
            (size_t)0,
            [](size_t total, Read const &rec) -> size_t
            {
                for (Dna4Sequence const &seq : rec.raw)
                    total += seq.size();
                for (Dna4Sequence const &seq : rec.rcomp)
                    total += seq.size();
                return total;
            });
    }

    size_t detail::countBps(const ReadChunks &data)
    {
        return std::accumulate(
            data.cbegin(),
            data.cend(),
            (size_t)0,
            [](size_t total, ReadVector const &v)->size_t { return total + detail::countBps(v); });
    }

    size_t detail::numReads(const ReadVector &data)
    {
        return std::accumulate(
            data.cbegin(),
            data.cend(),
            (size_t)0,
            [](size_t total, Read const &rec)->size_t { return total + rec.raw.size() + rec.rcomp.size(); });
    }

    size_t detail::numReads(const ReadChunks &data)
    {
        return std::accumulate(
            data.cbegin(),
            data.cend(),
            (size_t)0,
            [](size_t total, ReadVector const &v)->size_t { return total + detail::numReads(v); });
    }

    size_t detail::numEdges(const ReadChunks &data)
    {
        return detail::countBps(data) + /* end terminals */detail::numReads(data);
    }

    ReadVector
    flatten(ReadChunks &&data)
    {
        size_t totalSize = std::accumulate(data.cbegin(),
                                           data.cend(),
                                           (size_t)0,
                                           [](size_t l, ReadVector const& rv)->size_t { return l + rv.size(); });

        ReadVector out;
        out.reserve(totalSize);

        for (auto it = std::make_move_iterator(data.begin()),
                  end = std::make_move_iterator(data.end());
             it != end;
             ++it)
        {
            out.insert(out.end(),
                       std::make_move_iterator(it->begin()),
                       std::make_move_iterator(it->end()));
        }

        return out;
    }

    std::vector<ReadChunks>
    chunk(ReadChunks &&data, size_t granularity)
    {
        // index bps per chunk

        std::vector<size_t> sizes(data.size(), (size_t)0);
        std::transform(data.cbegin(), data.cend(),
                       sizes.begin(),
                       [](ReadVector const &chunk) -> size_t
                       {
                           return detail::countBps(chunk);
                       });

        size_t ntChunkSize = std::accumulate(sizes.cbegin(), sizes.cend(), (size_t)0) / granularity;

        // create chunks

        std::vector<ReadChunks> xdata;
        xdata.reserve(granularity);

        ReadChunks &currentChunk = xdata.emplace_back();
        size_t currentChunkSize = 0;

        auto nextChunkSizeIt = sizes.cbegin();
        for (auto it = std::make_move_iterator(data.begin()),
                  end = std::make_move_iterator(data.end());
             it != end;
             (void)++it, (void)++nextChunkSizeIt)
        {
            if (currentChunkSize >= ntChunkSize)
            {
                currentChunk = xdata.emplace_back();
                currentChunkSize = 0;
            }

            currentChunk.emplace_back(*it);
            currentChunkSize += *nextChunkSizeIt;
        }

        return xdata;
    }

} // namespace reads
