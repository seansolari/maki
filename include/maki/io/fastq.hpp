#pragma once
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <mutex>
#include <streambuf>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <oneapi/tbb/parallel_pipeline.h>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <maki/utils.hpp>
#include <maki/fasta.hpp>

namespace fs = std::filesystem;

using oneapi::tbb::filter_mode::parallel;
using oneapi::tbb::filter_mode::serial_in_order;

namespace reads
{

    struct DataFilePair
    {
        std::string
            forwardFile,
            reverseFile;

        DataFilePair(std::string &&fwd, std::string &&rev)
            : forwardFile(std::move(fwd)),
              reverseFile(std::move(rev))
        {
        }

        friend std::ostream &operator<<(std::ostream &, const DataFilePair &obj);

        inline bool operator==(DataFilePair const&) const =default;

        bool gzCompressed() const;
    };

    std::vector<DataFilePair> pairInputFiles(const std::vector<fs::path> &files);

    namespace detail
    {

        class CharBuffer
        {
        public:
            CharBuffer(size_t size)
                : _data(reinterpret_cast<char *>(std::malloc(size * sizeof(char)))),
                  _capacity(size),
                  _size(0)
            {
            }

            ~CharBuffer() { std::free(_data); }

            inline char *begin() noexcept { return _data; }

            inline const char *begin() const noexcept { return _data; }

            inline char *end() noexcept { return _data + _capacity; }

            inline const char *end() const noexcept { return _data + _capacity; }

            inline char *back() noexcept { return _data + _size; }

            inline const char *back() const noexcept { return _data + _size; }

            inline size_t size() const noexcept { return _size; }

            inline size_t capacity() const noexcept { return _capacity; }

            inline size_t bytesRemaining() const noexcept { return _capacity - _size; }

            inline bool empty() const { return _size == 0; }

            inline void clear() { _size = 0; }

            void flushTo(CharBuffer &dest);

            void swap(CharBuffer &rhs);

            inline constexpr size_t count(const char c) const { return std::count(_data, _data + _size, c); }

            inline constexpr size_t countTo(const char *end_, const char c) const { return std::count(const_cast<const char *>(_data), end_, c); }

            const char *prev(const char c) const;

            const char *nPrevFrom(const char *end_, size_t n, const char c) const;

            void write(const char *src, size_t n);

            void read(std::istream &);

            void flushOverflow(const char *end_, CharBuffer &out);

            inline char operator[](size_t i_) const { return _data[i_]; }

        private:
            char *_data;

            size_t
                _capacity,
                _size;
        };

        std::ostream &operator<<(std::ostream &, const CharBuffer &obj);

        template <typename CharT, typename TraitsT = std::char_traits<CharT>>
        class basic_spanbuf : public std::basic_streambuf<CharT, TraitsT>
        {
        public:
            basic_spanbuf(CharT *p, size_t n_) { this->setg(p, p, p + n_); }

            basic_spanbuf(CharBuffer &buffer)
                requires std::is_same_v<CharT, char>
                : basic_spanbuf(buffer.begin(), buffer.size())
            {
            }
        };

        class BufferPair
        {
            friend class Operator_ReadChunk;

        public:
            BufferPair(size_t n) : fwd(n), rev(n) {}

            inline void swap(BufferPair &rhs)
            {
                fwd.swap(rhs.fwd);
                rev.swap(rhs.rev);
            }

            inline bool empty() const { return fwd.empty() && rev.empty(); }

            inline void clear()
            {
                fwd.clear();
                rev.clear();
            }

            inline CharBuffer const& forward() const noexcept { return fwd; }

            inline CharBuffer const& reverse() const noexcept { return rev; }

            inline basic_spanbuf<char> viewForward() { return basic_spanbuf<char>(fwd); }

            inline basic_spanbuf<char> viewReverse() { return basic_spanbuf<char>(rev); }

        protected:
            CharBuffer
                fwd,
                rev;
        };

        class SharedBufferPair : public BufferPair
        {
        public:
            SharedBufferPair(size_t n) : BufferPair(n), inUse(false), mtx() {}

            inline bool avail()
            {
                std::lock_guard<std::mutex> lock(mtx);
                return !inUse;
            }

            inline void setUnavail()
            {
                std::lock_guard<std::mutex> lock(mtx);
                inUse = true;
            }

            inline void setAvail()
            {
                std::lock_guard<std::mutex> lock(mtx);
                inUse = false;
            }

        protected:
            bool inUse;

        private:
            std::mutex mtx;
        };

        class SharedBufferPairVector
        {
            SharedBufferPair *_p;
            size_t
                _width,
                _length;

        public:
            SharedBufferPairVector(size_t numBuffers, size_t bufferSize)
                : _p(reinterpret_cast<SharedBufferPair *>(std::malloc(numBuffers * sizeof(SharedBufferPair)))),
                  _width(numBuffers),
                  _length(bufferSize)
            {
                for (size_t i = 0; i < _width; ++i)
                    std::construct_at(_p + i, bufferSize);
            }

            ~SharedBufferPairVector() { std::free(_p); }

            inline size_t numBuffers() const noexcept { return _width; }

            inline size_t bufferSize() const noexcept { return _length; }

            inline SharedBufferPair *data() const noexcept { return _p; }
        };

        class Operator_ReadChunk
        {
        public:
            Operator_ReadChunk(std::basic_istream<char> *fp_,
                               std::basic_istream<char> *rp_,
                               const char delim_,
                               SharedBufferPair *buffers_,
                               size_t numBuffers_,
                               BufferPair *ovfl_)
                : ifwd(fp_),
                  irev(rp_),
                  DELIM(delim_),
                  buffers(buffers_),
                  numBuffers(numBuffers_),
                  ovfl(ovfl_)
            {
            }

            SharedBufferPair * operator()(oneapi::tbb::flow_control &fc) const;

            SharedBufferPair * readChunk() const;

            inline bool completedRec(std::basic_istream<char> *const is) const {
                auto c_ = is->peek();
                return (c_ == DELIM) || (c_ == std::basic_istream<char>::traits_type::eof());
            }

        private:
            std::basic_istream<char>
                *const ifwd,
                *const irev;

            const char DELIM; /* '@' - fastq, '>' - fasta */

            SharedBufferPair *const buffers;
            size_t numBuffers;

            BufferPair *const ovfl;
        };

    } // namespace detail

    struct Read
    {
        std::string id;
        std::vector<Dna4Sequence>
            raw,
            rcomp;

        Read(std::string &&id_) : id(std::move(id_)), raw(), rcomp() {}

        void push(seqan3::dna5_vector &&vec, size_t minFragLength);

        inline bool operator==(Read const&) const =default;

        size_t numFragments() const noexcept;

        size_t rss() const;
    };

    using ReadVector = std::vector<Read>;

    class ReadChunks {
        std::vector<ReadVector> _data;

    public:
        template <class ...Args>
        ReadChunks(Args&& ...args)
            : _data(std::forward<Args>(args)...)
        {
        }

        inline auto cbegin() const noexcept { return _data.cbegin(); }

        inline auto cend() const noexcept { return _data.cend(); }

        inline auto begin() const noexcept { return _data.begin(); }

        inline auto end() const noexcept { return _data.end(); }

        inline auto begin() noexcept { return _data.begin(); }

        inline auto end() noexcept { return _data.end(); }

        inline size_t size() const noexcept { return _data.size(); }

        inline ReadVector & operator[](size_t i) { return _data[i]; }

        inline const ReadVector & operator[](size_t i) const { return _data[i]; }

        template <class ...Args>
        inline ReadVector& emplace_back(Args&& ...args) {
            return _data.emplace_back(std::forward<Args>(args)...);
        }

        size_t numFragments() const noexcept;

        size_t numTerminals(uint8_t k) const;

        size_t rss() const;
    };

    namespace parse
    {

        struct Operator_ParseChunk
        {
            Operator_ParseChunk(tbb::enumerable_thread_specific<ReadVector> *b_,
                                size_t minReadLength_,
                                size_t minFragLength_)
                : buffers(b_),
                  minReadLength(minReadLength_),
                  minFragLength(minFragLength_)
            {
            }

        public:
            void operator()(detail::SharedBufferPair *) const;

        protected:
            tbb::enumerable_thread_specific<ReadVector> *buffers;
            size_t
                minReadLength,
                minFragLength;
        };

        ReadChunks parsePairedFastq(DataFilePair const &fp, size_t minReadLength, size_t minFragLength, size_t threads, size_t buffersize = 10 * 1024 * 1024 /*10 Mb buffer*/);
        ReadChunks parsePairedFastq(std::basic_istream<char> *ifwd, std::basic_istream<char> *irev, size_t minReadLength, size_t minFragLength, size_t threads, size_t buffersize = 10 * 1024 * 1024 /*10 Mb buffer*/);

    } // namespace parse

    namespace detail {

        size_t countBps(ReadVector const &data);

        size_t countBps(ReadChunks const &data);

        // counts number of contained fragments (a 'read' can be split into
        // multiple fragments by an 'N').
        size_t numReads(ReadVector const &data);

        // counts number of contained fragments
        size_t numReads(ReadChunks const &data);

        size_t numEdges(ReadChunks const &data);

    } // namespace detail

    struct UnwindReads
    {
        template <typename Fn_>
            requires CallsRange<Fn_>
        void forEach(ReadChunks const &data, Fn_ f) const
        {
            for (ReadVector const &rv : data)
            {
                for (Read const &read : rv)
                {
                    for (Dna4Sequence const &seq : read.raw)
                        f(seq.cbegin(), seq.cend(), true);
                    for (Dna4Sequence const &seq : read.rcomp)
                        f(seq.cbegin(), seq.cend(), true);
                }
            }
        }

        template <typename Fn_>
        void forEachBegin(ReadChunks const &data, Fn_ f) const
        {
            for (ReadVector const &rv : data)
            {
                for (Read const &read : rv)
                {
                    for (Dna4Sequence const &seq : read.raw)
                        f(seq.cbegin());
                    for (Dna4Sequence const &seq : read.rcomp)
                        f(seq.cbegin());
                }
            }
        }
    };

} // namespace reads

template <>
struct UnwindTraits<reads::UnwindReads>
{
    typedef reads::ReadChunks range_type;
};

namespace reads
{

    ReadVector flatten(ReadChunks &&data);

    std::vector<ReadChunks> chunk(ReadChunks &&data, size_t granularity);

} // namespace reads
