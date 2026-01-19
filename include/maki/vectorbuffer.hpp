#pragma once

// Modified from:
// SPDX-FileCopyrightText: 2023 Stephan Lachnit
// SPDX-License-Identifier: MIT

#include <cassert>
#include <memory>
#include <sstream>
#include <streambuf>
#include <string>
#include <ostream>
#include <vector>
#include <zstd.h> // presumes zstd library is installed
#include <sdsl/vectors.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>
#include <glog/logging.h>
#include <maki/maki.h>

template <class CharT = char, class Traits = std::char_traits<CharT>>
class vectorbuf : public std::basic_streambuf<CharT>
{
public:
    using streambuf = std::basic_streambuf<CharT, Traits>;
    using char_type = typename streambuf::char_type;
    using int_type = typename streambuf::int_type;
    using traits_type = typename streambuf::traits_type;
    using vector = std::vector<char_type>;
    using value_type = typename vector::value_type;
    using size_type = typename vector::size_type;

    // Constructor for vecbuf with optional initial capacity
    vectorbuf(size_type capacity = 0) : vector_() { reserve(capacity); }

    // Clear vector contents and set output pointers
    constexpr void clear()
    {
        vector_.clear();
        setp_from_vector();
    }

    // Forwarder for std::vector::reserve
    constexpr void reserve(size_type capacity)
    {
        vector_.reserve(capacity);
        setp_from_vector();
    }

    // Increase the capacity of the buffer by reserving the current_size + additional_capacity
    constexpr void reserve_additional(size_type additional_capacity) { reserve(size() + additional_capacity); }

    // Forwarder for std::vector::data
    constexpr const value_type *data() const { return vector_.data(); }

    // Forwarder for std::vector::size
    constexpr size_type size() const { return vector_.size(); }

    // Forwarder for std::vector::capacity
    constexpr size_type capacity() const { return vector_.capacity(); }

    // Implements std::basic_streambuf::xsputn
    std::streamsize xsputn(const char_type *s, std::streamsize count) override
    {
        try
        {
            reserve_additional(count);
        }
        catch (const std::bad_alloc &error)
        {
            // reserve did not work, use slow algorithm
            return xsputn_slow(s, count);
        }
        // reserve worked, use fast algorithm
        return xsputn_fast(s, count);
    }

protected:
    // Calculates value to std::basic_streambuf::pbase from vector
    constexpr value_type *pbase_from_vector() const { return const_cast<value_type *>(vector_.data()); }

    // Calculates value to std::basic_streambuf::pptr from vector
    constexpr value_type *pptr_from_vector() const { return const_cast<value_type *>(vector_.data() + vector_.size()); }

    // Calculates value to std::basic_streambuf::epptr from vector
    constexpr value_type *epptr_from_vector() const { return const_cast<value_type *>(vector_.data()) + vector_.capacity(); }

    // Sets the values for std::basic_streambuf::pbase, std::basic_streambuf::pptr and std::basic_streambuf::epptr from vector
    constexpr void setp_from_vector()
    {
        streambuf::setp(pbase_from_vector(), epptr_from_vector());
        streambuf::pbump(size());
    }

private:
    // std::vector containing the data
    vector vector_;

    // Fast implementation of std::basic_streambuf::xsputn if reserve_additional(count) succeeded
    std::streamsize xsputn_fast(const char_type *s, std::streamsize count)
    {
        // store current pptr (end of vector location)
        auto *old_pptr = pptr_from_vector();
        // resize the vector, does not move since space already reserved
        vector_.resize(vector_.size() + count);
        // directly memcpy new content to old pptr (end of vector before it was resized)
        traits_type::copy(old_pptr, s, count);
        // reserve() already calls setp_from_vector(), only adjust pptr to new epptr
        streambuf::pbump(count);

        return count;
    }

    // Slow implementation of std::basic_streambuf::xsputn if reserve_additional(count) did not succeed, might calls std::basic_streambuf::overflow()
    std::streamsize xsputn_slow(const char_type *s, std::streamsize count)
    {
        // reserving entire vector failed, emplace char for char
        std::streamsize written = 0;
        while (written < count)
        {
            try
            {
                // copy one char, should throw eventually std::bad_alloc
                vector_.emplace_back(s[written]);
            }
            catch (const std::bad_alloc &error)
            {
                // try overflow(), if eof return, else continue
                int_type c = this->overflow(traits_type::to_int_type(s[written]));
                if (traits_type::eq_int_type(c, traits_type::eof()))
                {
                    return written;
                }
            }
            // update pbase, pptr and epptr
            setp_from_vector();
            written++;
        }
        return written;
    }
};

struct IntVectorStats
{
    size_t
        numElements,
        elementWidth,
        numWords;
};

struct CompressedIntVector
{
    // original size
    IntVectorStats info;

    // compressed size and data
    size_t cxSize;
    std::vector<unsigned char> cxData;

    template <typename Vec>
    inline void setInfo(const Vec &v)
    {
        info.numElements = v.size();
        info.elementWidth = v.width();
        info.numWords = (v.bit_size() + 63) >> 6;
    }

    // compression ratio
    inline float cxBitsPerWord() const
    {
        return static_cast<double>(cxSize) / (8.0l * static_cast<double>(info.numWords));
    }

    constexpr inline size_t size() const noexcept { return cxSize; }

    constexpr inline void resize(size_t newSize)
    {
        cxSize = newSize;
        cxData.resize(newSize);
    }

    constexpr inline const unsigned char *data() const noexcept { return cxData.data(); }

    constexpr inline unsigned char *data() noexcept { return cxData.data(); }

    inline void shrink(size_t newSize)
    {
        assert(newSize <= cxSize);
        resize(newSize);
    }
};

struct ZstdDeleter
{
    inline void operator()(ZSTD_CCtx *cctx) const { ZSTD_freeCCtx(cctx); }
    inline void operator()(ZSTD_DCtx *dctx) const { ZSTD_freeDCtx(dctx); }
};

void checkZstd(const size_t err, const char *fail_msg);

struct ZstdCompressor
{
    vectorbuf<> bytes;
    std::unique_ptr<ZSTD_CCtx, ZstdDeleter> cctx;
    CompressedIntVector rsl;

    ZstdCompressor();

    ZstdCompressor(const ZstdCompressor &) = delete;
    ZstdCompressor &operator=(const ZstdCompressor &) = delete;

    template <uint8_t t_width>
    void compress(sdsl::int_vector<t_width> const &vec)
    {
        bytes.clear();

        // serialise to bytes
        {
            // initialise stream state

            std::ostream os(&bytes);

            // serialise archive

            cereal::BinaryOutputArchive oarchive(os);
            oarchive(vec);
        } // archive goes out of scope, ensuring all contents are flushed

        // reset output buffer

        rsl.setInfo(vec);
        rsl.resize(ZSTD_compressBound(bytes.size()));

        // compress data

        const size_t cSize = ZSTD_compressCCtx(cctx.get(), rsl.data(), rsl.size(), bytes.data(), bytes.size(), ZSTD_COMPRESSION_LEVEL);
        checkZstd(cSize, "ZSTD_compressCCtx failed!");
        rsl.shrink(cSize);
    }
};

struct ZstdDecompressor
{
    CompressedIntVector src;
    std::unique_ptr<ZSTD_DCtx, ZstdDeleter> dctx;
    std::string serialData;

    ZstdDecompressor()
        : src(),
          dctx(ZSTD_createDCtx(), ZstdDeleter()),
          serialData()
    {
        if (dctx.get() == NULL)
        {
            LOG(ERROR) << "ZSTD_createDCtx() failed!";
            throw std::bad_alloc();
        }
    }

    template <uint8_t t_width>
    void decompress(sdsl::int_vector<t_width> &arr)
    {
        // allocate decompression buffer

        size_t serialLen = ZSTD_getDecompressedSize(src.data(), src.size());
        serialData.resize(serialLen);

        // decompress

        const size_t dSize = ZSTD_decompressDCtx(dctx.get(), serialData.data(), serialLen, src.data(), src.size());
        checkZstd(dSize, "ZSTD_decompressDCtx() failed!");

        // deserialise - move data into stream

        std::istringstream is(std::move(serialData));
        {
            cereal::BinaryInputArchive iarchive(is);
            iarchive(arr);
        } // archive goes out of scope, ensuring all contents are flushed

        // retake ownership of data stream

        serialData = std::move(is).str();
    }
};
