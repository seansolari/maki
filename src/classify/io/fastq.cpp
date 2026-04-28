
#include "maki/classify/io/fastq.hpp"
#include "maki/core/seq/io.hpp"
#include <algorithm>
#include <cassert>
#include <numeric>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/utility/views/zip.hpp>
#include <zstr.hpp>

using namespace seqan3::literals;

std::ostream &operator<<(std::ostream &os, const DataFilePair &obj) {
  return os << "PairedReads(" << obj.forwardFile << ", " << obj.reverseFile
            << ")";
}

bool DataFilePair::gzCompressed() const {
  bool lCx = forwardFile.ends_with(".gz"), rCx = reverseFile.ends_with(".gz");

  if (lCx != rCx)
    throw std::runtime_error("paired files " + forwardFile + " and " +
                             reverseFile + " have different compression");
  else
    return lCx;
}

std::vector<DataFilePair> pairInputFiles(const std::vector<fs::path> &files) {
  // check even number of files

  assert((files.size() & 1u) == 0u);
  std::size_t numPairs = files.size() >> 1;

  // sort files lexicographically

  std::vector<std::size_t> fileIndex(files.size());
  std::iota(fileIndex.begin(), fileIndex.end(), 0);
  std::sort(fileIndex.begin(), fileIndex.end(),
            [&](std::size_t i, std::size_t j) { return files[i] < files[j]; });

  // pair and check files

  std::vector<DataFilePair> pairs;
  pairs.reserve(numPairs);

  for (std::size_t i = 0; i < numPairs; ++i) {
    std::string fwd = files[fileIndex[2 * i]].string(),
                rev = files[fileIndex[(2 * i) + 1]].string();

    std::string_view fwdName = extractSequenceName(fwd),
                     revName = extractSequenceName(rev);

    if (fwdName.back() != '1') {
      throw std::runtime_error("Paired files " + fwd + " and " + rev +
                               ", but " + fwd + " does not end in 1.");
    } else if (revName.back() != '2') {
      throw std::runtime_error("Paired files " + fwd + " and " + rev +
                               ", but " + rev + " does not end in 2.");
    } else if (fwdName.substr(0, fwdName.size() - 1) !=
               revName.substr(0, revName.size() - 1)) {
      throw std::runtime_error("Paired files " + fwd + " and " + rev +
                               ", but prefixes do not match.");
    } else {
      pairs.emplace_back(std::move(fwd), std::move(rev));
    }
  }

  return pairs;
}

detail::CharBuffer::CharBuffer(std::size_t size)
    : _data(reinterpret_cast<char *>(std::malloc(size * sizeof(char)))),
      _capacity(size), _size(0) {}

void detail::CharBuffer::flushTo(detail::CharBuffer &dest) {
  assert(_size < dest.bytesRemaining());
  std::copy(_data, _data + _size, dest.back());
  dest._size += _size;
  _size = 0;
}

void detail::CharBuffer::swap(detail::CharBuffer &rhs) {
  std::swap(_data, rhs._data);
  std::swap(_capacity, rhs._capacity);
  std::swap(_size, rhs._size);
}

const char *detail::CharBuffer::prev(char c) const {
  auto _rend = std::make_reverse_iterator(begin()),
       _it = std::find(std::make_reverse_iterator(back()), _rend, c);
  if (_it == _rend)
    return back();
  else
    return begin() + (_rend - _it) - 1;
}

const char *detail::CharBuffer::nPrevFrom(const char *end_, std::size_t n,
                                          char c) const {
  auto _rend = std::make_reverse_iterator(begin()),
       _it = std::find_if(std::make_reverse_iterator(end_), _rend,
                          [&](char const &z_) { return (z_ == c) && !(--n); });
  if (_it == _rend)
    return back();
  else
    return begin() + (_rend - _it) - 1;
}

void detail::CharBuffer::write(const char *src, std::size_t n) {
  assert(_size + n < _capacity);
  std::copy_n(src, n, back());
  _size += n;
}

void detail::CharBuffer::read(std::istream &is) {
  is.read(back(), bytesRemaining());
  _size += is.gcount();
}

void detail::CharBuffer::flushOverflow(const char *end_,
                                       detail::CharBuffer &out) {
  std::size_t newsize_ = end_ - _data;
  assert(newsize_ < _size);
  out.write(end_, _size - newsize_);
  _size = newsize_;
}

std::ostream &detail::operator<<(std::ostream &os,
                                 const detail::CharBuffer &obj) {
  os << "CharBuffer[\n";
  for (std::size_t i = 0; i < obj.size(); ++i)
    os << obj[i];
  os << "]";
  return os;
}

constexpr std::size_t detail::CharBuffer::count(const char c) const {
  return std::count(_data, _data + _size, c);
}

constexpr std::size_t detail::CharBuffer::countTo(const char *end_,
                                                  const char c) const {
  return std::count(const_cast<const char *>(_data), end_, c);
}

void detail::BufferPair::swap(BufferPair &rhs) {
  fwd.swap(rhs.fwd);
  rev.swap(rhs.rev);
}

bool detail::BufferPair::empty() const { return fwd.empty() && rev.empty(); }

void detail::BufferPair::clear() {
  fwd.clear();
  rev.clear();
}

bool detail::SharedBufferPair::avail() {
  std::lock_guard<std::mutex> lock(mtx);
  return !inUse;
}

void detail::SharedBufferPair::setUnavail() {
  std::lock_guard<std::mutex> lock(mtx);
  inUse = true;
}

void detail::SharedBufferPair::setAvail() {
  std::lock_guard<std::mutex> lock(mtx);
  inUse = false;
}

detail::SharedBufferPairVector::SharedBufferPairVector(std::size_t numBuffers,
                                                       std::size_t bufferSize)
    : _p(reinterpret_cast<SharedBufferPair *>(
          std::malloc(numBuffers * sizeof(SharedBufferPair)))),
      _width(numBuffers), _length(bufferSize) {
  for (std::size_t i = 0; i < _width; ++i)
    std::construct_at(_p + i, bufferSize);
}

detail::Operator_ReadChunk::Operator_ReadChunk(std::basic_istream<char> *fp_,
                                               std::basic_istream<char> *rp_,
                                               const char delim_,
                                               SharedBufferPair *buffers_,
                                               std::size_t numBuffers_,
                                               BufferPair *ovfl_)
    : ifwd(fp_), irev(rp_), DELIM(delim_), buffers(buffers_),
      numBuffers(numBuffers_), ovfl(ovfl_) {}

detail::SharedBufferPair *detail::Operator_ReadChunk::operator()(
    tbb::detail::d1::flow_control &fc) const {
  if (ifwd->eof() && irev->eof()) {
    if (!ovfl->empty())
      throw std::runtime_error("paired-end structure is compromised");

    fc.stop();
    return nullptr;
  } else
    return readChunk();
}

detail::SharedBufferPair *detail::Operator_ReadChunk::readChunk() const {
  // find next available buffer

  for (SharedBufferPair *bp = buffers; bp < buffers + numBuffers; ++bp) {
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

    const char *fend =
                   completedRec(ifwd) ? bp->fwd.back() : bp->fwd.prev(DELIM),
               *rend =
                   completedRec(irev) ? bp->rev.back() : bp->rev.prev(DELIM);

    if ((fend == bp->fwd.begin()) || (rend == bp->rev.begin()))
      throw std::runtime_error(
          "buffer size not begin enough to parse read record");

    std::size_t fc = bp->fwd.countTo(fend, DELIM),
                rc = bp->rev.countTo(rend, DELIM), recs = std::min(fc, rc);

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

bool detail::Operator_ReadChunk::completedRec(
    std::basic_istream<char> *const is) const {
  auto c_ = is->peek();
  return (c_ == DELIM) || (c_ == std::basic_istream<char>::traits_type::eof());
}

void Read::push(seqan3::dna5_vector &&data, std::size_t minFragLength) {
  auto lengthFilter = std::views::filter([&](auto &&seq) {
    return std::ranges::size(std::forward<decltype(seq)>(seq)) >= minFragLength;
  });

  auto Dna5ToDna4 = std::views::transform(
      [](auto &&
             c) { // https://docs.seqan.de/seqan3/main_user/cookbook.html#cookbook_convert_alphabet_range
        return static_cast<seqan3::dna4>(std::forward<decltype(c)>(c));
      });

  for (auto &&subSeq : data | std::views::split('N'_dna5) | lengthFilter) {
    Dna4Sequence &f = sequences[0].emplace_back(),
                 &r = sequences[1].emplace_back();

    f.reserve(subSeq.size());
    for (auto &&c : subSeq | Dna5ToDna4)
      f.push_back(c);

    r.reserve(subSeq.size());
    for (auto &&c :
         subSeq | std::views::reverse | seqan3::views::complement | Dna5ToDna4)
      r.push_back(c);
  }
}

std::size_t Read::numKmers(std::size_t k) const {
  std::size_t total = 0;
  for (const auto &seq : sequences[0]) {
    assert(seq.size() >= k);
    total += seq.size() - k + 1;
  }
  for (const auto &seq : sequences[1]) {
    assert(seq.size() >= k);
    total += seq.size() - k + 1;
  }
  return total;
}

std::size_t Read::numFragments() const noexcept {
  return sequences[0].size() + sequences[1].size();
}

std::size_t Read::rss() const {
  std::size_t tot = (id.capacity() * sizeof(std::string::value_type)) +
                    (sequences[0].capacity() * sizeof(Dna4Sequence)) +
                    (sequences[1].capacity() * sizeof(Dna4Sequence));
  for (Dna4Sequence const &seq : sequences[0])
    tot += ((seq.capacity() + 63) / 64) * 8;
  for (Dna4Sequence const &seq : sequences[1])
    tot += ((seq.capacity() + 63) / 64) * 8;
  return tot;
}

std::size_t ReadChunks::numTerminals(std::size_t k) const {
  return (std::size_t)k * numFragments();
}

poly_input_range<SequenceFragment> ReadChunks::terminals() const {
  return poly_input_range<SequenceFragment>(
      data | std::views::join |
      std::views::transform([](const Read &r) { return r.view(); }) |
      std::views::join | std::views::transform([](const Dna4Sequence &seq) {
        return SequenceFragment(seq.cbegin(), seq.cend(), 0, true);
      }));
}

std::size_t ReadChunks::numKmers(std::size_t k) const {
  return std::accumulate(
      data.cbegin(), data.cend(), (std::size_t)0,
      [k](std::size_t total, const ReadVector &v) -> std::size_t {
        return total + detail::numKmers(v, k);
      });
}

poly_input_range<SequenceFragment>
ReadChunks::fragments([[maybe_unused]] std::size_t k) const {
  return poly_input_range<SequenceFragment>(
      data | std::views::join |
      std::views::transform([](const Read &r) { return r.view(); }) |
      std::views::join | std::views::transform([](const Dna4Sequence &seq) {
        return SequenceFragment(seq.cbegin(), seq.cend(), 0, true);
      }));
}

std::size_t ReadChunks::countBps() const {
  return std::accumulate(
      data.cbegin(), data.cend(), (std::size_t)0,
      [](std::size_t total, const ReadVector &v) -> std::size_t {
        return total + detail::countBps(v);
      });
}

std::size_t ReadChunks::numFragments() const {
  return std::accumulate(
      data.cbegin(), data.cend(), (std::size_t)0,
      [](std::size_t total, const ReadVector &v) -> std::size_t {
        return total + detail::numFragments(v);
      });
}

size_t ReadChunks::rss() const {
  size_t tot = 0;
  for (const auto &ch : data)
    for (const Read &r : ch)
      tot += r.rss();
  return tot;
}

detail::Operator_ParseChunk::Operator_ParseChunk(
    tbb::enumerable_thread_specific<ReadVector> *b_, std::size_t minReadLength_,
    std::size_t minFragLength_)
    : buffers(b_), minReadLength(minReadLength_),
      minFragLength(minFragLength_) {}

void detail::Operator_ParseChunk::operator()(
    detail::SharedBufferPair *in) const {
  ReadVector &out = buffers->local();

  // parse paired data

  detail::basic_spanbuf fwdView = in->viewForward(),
                        revView = in->viewReverse();

  std::istream fwdStream(&fwdView), revStream(&revView);

  seqan3::sequence_file_input ifwd{fwdStream, seqan3::format_fastq{}},
      irev{revStream, seqan3::format_fastq{}};

  for (auto &&[rec1, rec2] : seqan3::views::zip(ifwd, irev)) {
    if ((rec1.sequence().size() >= minReadLength) &&
        (rec2.sequence().size() >= minReadLength)) {
      Read &v = out.emplace_back(std::move(rec1.id()));
      v.push(std::move(rec1).sequence(), minFragLength);
      v.push(std::move(rec2).sequence(), minFragLength);
    }
  }

  // free data

  in->clear();
  in->setAvail();
}

ReadChunks detail::parsePairedFastq(const DataFilePair &fp,
                                    std::size_t minReadLength,
                                    std::size_t minFragLength,
                                    std::size_t threads,
                                    std::size_t buffersize) {
  if (fp.gzCompressed()) {
    zstr::ifstream ifwd(fp.forwardFile), irev(fp.reverseFile);
    return parsePairedFastq(&ifwd, &irev, minReadLength, minFragLength, threads,
                            buffersize);
  } else {
    std::basic_ifstream<char> ifwd(fp.forwardFile), irev(fp.reverseFile);
    return parsePairedFastq(&ifwd, &irev, minReadLength, minFragLength, threads,
                            buffersize);
  }
}

ReadChunks detail::parsePairedFastq(std::basic_istream<char> *ifwd,
                                    std::basic_istream<char> *irev,
                                    std::size_t minReadLength,
                                    std::size_t minFragLength,
                                    std::size_t threads,
                                    std::size_t buffersize) {
  detail::SharedBufferPairVector buffers(threads, buffersize);
  detail::BufferPair ovflw(buffersize);

  tbb::enumerable_thread_specific<ReadVector> threadOutput;

  detail::Operator_ReadChunk ReadChunk(ifwd, irev, '@', buffers.data(),
                                       buffers.numBuffers(), &ovflw);
  Operator_ParseChunk ParseChunk(&threadOutput, minReadLength, minFragLength);

  oneapi::tbb::parallel_pipeline(
      threads,
      oneapi::tbb::make_filter<void, detail::SharedBufferPair *>(
          serial_in_order, ReadChunk),
      oneapi::tbb::make_filter<detail::SharedBufferPair *, void>(parallel,
                                                                 ParseChunk));

  return ReadChunks(std::make_move_iterator(threadOutput.begin()),
                    std::make_move_iterator(threadOutput.end()));
}

std::size_t detail::countBps(const ReadVector &data) {
  return std::accumulate(data.cbegin(), data.cend(), (std::size_t)0,
                         [](std::size_t total, Read const &rec) -> std::size_t {
                           for (Dna4Sequence const &seq : rec.sequences[0])
                             total += seq.size();
                           for (Dna4Sequence const &seq : rec.sequences[1])
                             total += seq.size();
                           return total;
                         });
}

std::size_t detail::numKmers(const ReadVector &v, std::size_t k) {
  return std::accumulate(v.cbegin(), v.cend(), (std::size_t)0,
                         [k](std::size_t total, const Read &r) -> std::size_t {
                           return total + r.numKmers(k);
                         });
}

std::size_t detail::numFragments(const ReadVector &data) {
  return std::accumulate(data.cbegin(), data.cend(), (std::size_t)0,
                         [](std::size_t total, const Read &rec) -> std::size_t {
                           return total + rec.sequences[0].size() +
                                  rec.sequences[1].size();
                         });
}

ReadVector flattenReads(ReadChunks &&chunks) {
  std::size_t totalSize =
      std::accumulate(chunks.data.cbegin(), chunks.data.cend(), (std::size_t)0,
                      [](std::size_t l, ReadVector const &rv) -> std::size_t {
                        return l + rv.size();
                      });

  ReadVector out;
  out.reserve(totalSize);

  for (auto it = std::make_move_iterator(chunks.data.begin()),
            end = std::make_move_iterator(chunks.data.end());
       it != end; ++it) {
    out.insert(out.end(), std::make_move_iterator(it->begin()),
               std::make_move_iterator(it->end()));
  }

  return out;
}

std::vector<ReadChunks> chunkReads(ReadChunks &&chunks,
                                   std::size_t granularity) {
  // index bps per chunk

  std::vector<std::size_t> sizes(chunks.data.size(), (std::size_t)0);
  std::transform(chunks.data.cbegin(), chunks.data.cend(), sizes.begin(),
                 [](ReadVector const &chunk) -> std::size_t {
                   return detail::countBps(chunk);
                 });

  std::size_t ntChunkSize =
      std::accumulate(sizes.cbegin(), sizes.cend(), (std::size_t)0) /
      granularity;

  // create chunks

  std::vector<ReadChunks> xdata;
  xdata.reserve(granularity);

  ReadChunks &currentChunk = xdata.emplace_back();
  std::size_t currentChunkSize = 0;

  auto nextChunkSizeIt = sizes.cbegin();
  for (auto it = std::make_move_iterator(chunks.data.begin()),
            end = std::make_move_iterator(chunks.data.end());
       it != end; (void)++it, (void)++nextChunkSizeIt) {
    if (currentChunkSize >= ntChunkSize) {
      currentChunk = xdata.emplace_back();
      currentChunkSize = 0;
    }

    currentChunk.data.emplace_back(*it);
    currentChunkSize += *nextChunkSizeIt;
  }

  return xdata;
}
