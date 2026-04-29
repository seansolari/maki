#pragma once
#include "maki/core/seq/concepts.hpp"
#include "maki/core/seq/io.hpp"
#include <filesystem>
#include <mutex>
#include <string>
#include <vector>

#include <oneapi/tbb/enumerable_thread_specific.h>
#include <oneapi/tbb/parallel_pipeline.h>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>

namespace fs = std::filesystem;
using oneapi::tbb::filter_mode::parallel;
using oneapi::tbb::filter_mode::serial_in_order;

struct DataFilePair {
  std::string forwardFile, reverseFile;

  DataFilePair(std::string &&fwd, std::string &&rev)
      : forwardFile(std::move(fwd)), reverseFile(std::move(rev)) {}

  friend std::ostream &operator<<(std::ostream &, const DataFilePair &obj);
  inline bool operator==(DataFilePair const &) const = default;
  bool gzCompressed() const;
};

std::vector<DataFilePair> pairInputFiles(const std::vector<fs::path> &files);

namespace detail {

class CharBuffer {
public:
  CharBuffer(std::size_t size);
  ~CharBuffer() { std::free(_data); }
  inline char *begin() noexcept { return _data; }
  inline const char *begin() const noexcept { return _data; }
  inline char *end() noexcept { return _data + _capacity; }
  inline const char *end() const noexcept { return _data + _capacity; }
  inline char *back() noexcept { return _data + _size; }
  inline const char *back() const noexcept { return _data + _size; }
  inline std::size_t size() const noexcept { return _size; }
  inline std::size_t capacity() const noexcept { return _capacity; }
  std::size_t bytesRemaining() const noexcept;
  inline bool empty() const { return _size == 0; }
  inline void clear() { _size = 0; }
  void flushTo(CharBuffer &dest);
  void swap(CharBuffer &rhs);
  const char *prev(const char c) const;
  const char *nPrevFrom(const char *end_, std::size_t n, const char c) const;
  void write(const char *src, std::size_t n);
  void read(std::istream &);
  void flushOverflow(const char *end_, CharBuffer &out);
  inline char operator[](std::size_t i_) const { return _data[i_]; }
  std::size_t count(const char c) const;
  std::size_t countTo(const char *end_, const char c) const;

private:
  char *_data;
  std::size_t _capacity, _size;
};

std::ostream &operator<<(std::ostream &, const CharBuffer &obj);

template <typename CharT, typename TraitsT = std::char_traits<CharT>>
class basic_spanbuf : public std::basic_streambuf<CharT, TraitsT> {
public:
  basic_spanbuf(CharT *p, std::size_t n_) { this->setg(p, p, p + n_); }

  basic_spanbuf(CharBuffer &buffer)
    requires std::is_same_v<CharT, char>
      : basic_spanbuf(buffer.begin(), buffer.size()) {}
};

class BufferPair {
  friend class Operator_ReadChunk;

public:
  BufferPair(std::size_t n) : fwd(n), rev(n) {}
  void swap(BufferPair &rhs);
  bool empty() const;
  void clear();
  inline CharBuffer const &forward() const noexcept { return fwd; }
  inline CharBuffer const &reverse() const noexcept { return rev; }
  inline basic_spanbuf<char> viewForward() { return basic_spanbuf<char>(fwd); }
  inline basic_spanbuf<char> viewReverse() { return basic_spanbuf<char>(rev); }

protected:
  CharBuffer fwd, rev;
};

class SharedBufferPair : public BufferPair {
public:
  SharedBufferPair(std::size_t n) : BufferPair(n), inUse(false), mtx() {}
  bool avail();
  void setUnavail();
  void setAvail();

protected:
  bool inUse;

private:
  std::mutex mtx;
};

class SharedBufferPairVector {
  SharedBufferPair *_p;
  std::size_t _width, _length;

public:
  SharedBufferPairVector(std::size_t numBuffers, std::size_t bufferSize);
  ~SharedBufferPairVector() { std::free(_p); }
  inline std::size_t numBuffers() const noexcept { return _width; }
  inline std::size_t bufferSize() const noexcept { return _length; }
  inline SharedBufferPair *data() const noexcept { return _p; }
};

class Operator_ReadChunk {
public:
  Operator_ReadChunk(std::basic_istream<char> *fp_,
                     std::basic_istream<char> *rp_, const char delim_,
                     SharedBufferPair *buffers_, std::size_t numBuffers_,
                     BufferPair *ovfl_);
  SharedBufferPair *operator()(oneapi::tbb::flow_control &fc) const;
  SharedBufferPair *readChunk() const;
  bool completedRec(std::basic_istream<char> *const is) const;

private:
  std::basic_istream<char> *const ifwd, *const irev;
  const char DELIM; /* '@' - fastq, '>' - fasta */
  SharedBufferPair *const buffers;
  std::size_t numBuffers;
  BufferPair *const ovfl;
};

} // namespace detail

struct Read {
  std::string id;
  std::array<std::vector<Dna4Sequence>, 2> sequences;
  Read(std::string &&id_) : id(std::move(id_)), sequences() {}
  void push(seqan3::dna5_vector &&vec, std::size_t minFragLength);
  auto view() const { return sequences | std::views::join; }
  inline bool operator==(const Read &) const = default;
  std::size_t numKmers(std::size_t k) const;
  std::size_t numFragments() const noexcept;
  std::size_t rss() const;
};

using ReadVector = std::vector<Read>;

struct ReadChunks final : public SequenceContainer {
  std::vector<ReadVector> data;

  template <class... Args>
  ReadChunks(Args &&...args) : data(std::forward<Args>(args)...) {}

  virtual std::size_t numTerminals(std::size_t k) const override final;
  virtual poly_input_range<SequenceFragment> terminals() const override final;
  virtual std::size_t numKmers(std::size_t k) const override final;
  virtual poly_input_range<SequenceFragment>
  fragments(std::size_t k) const override final;

  std::size_t numFragments() const;
  std::size_t countBps() const;
  std::size_t rss() const;
};

namespace detail {

struct Operator_ParseChunk {
  Operator_ParseChunk(tbb::enumerable_thread_specific<ReadVector> *b_,
                      std::size_t minReadLength_, std::size_t minFragLength_);

public:
  void operator()(detail::SharedBufferPair *) const;

protected:
  tbb::enumerable_thread_specific<ReadVector> *buffers;
  std::size_t minReadLength, minFragLength;
};

ReadChunks parsePairedFastq(const DataFilePair &fp, std::size_t minReadLength,
                            std::size_t minFragLength, std::size_t threads,
                            std::size_t buffersize = 10 * 1024 *
                                                     1024 /*10 Mb buffer*/);
ReadChunks
parsePairedFastq(std::basic_istream<char> *ifwd, std::basic_istream<char> *irev,
                 std::size_t minReadLength, std::size_t minFragLength,
                 std::size_t threads,
                 std::size_t buffersize = 10 * 1024 * 1024 /*10 Mb buffer*/);

std::size_t numKmers(const ReadVector &v, std::size_t k);
std::size_t countBps(const ReadVector &data);

// counts number of contained fragments (a 'read' can be split into
// multiple fragments by an 'N').
std::size_t numFragments(const ReadVector &data);

} // namespace detail

ReadVector flattenReads(ReadChunks &&data);

// `granularity` is roughly the length of the resulting vector
std::vector<ReadChunks> chunkReads(ReadChunks &&data, std::size_t granularity);
