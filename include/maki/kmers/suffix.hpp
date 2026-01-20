#pragma once

#include <cassert>
#include <cstddef>
#include <numeric>
#include <vector>
#include <oneapi/tbb/parallel_for.h>
#include <maki/fasta.hpp>
#include <maki/utils.hpp>

template <typename ChunkType_, typename Proj_>
std::vector<size_t>
accumulatePartials(std::vector<ChunkType_> const &chunks,
                   Proj_ proj)
{
    std::vector<size_t> partials(chunks.size(), 0);
    auto out = partials.begin();
    for (ChunkType_ const &obj : chunks)
    {
        *out = proj(obj);
        ++out;
    }

    std::partial_sum(partials.cbegin(), partials.cend(),
                     partials.begin());
    return partials;
}

namespace suffix
{

    /**
     * Represent nt using the 2-bit encoding: `[A, C, G, T] <=> [0, 1, 2, 3]`.
     * These are used for efficiently representing suffixes in a single cell
     * so that they can be quickly counted or compared. Most significant
     * base-pair of each window is the one most recently inserted.
     */
    struct SmallRollingNuclSeq
    {
        friend class SuffixTable;

    public:
        using window_t = uint64_t;
        static constexpr uint8_t bitsPerNucl = 2;

        // constructors

        SmallRollingNuclSeq(size_t _size, window_t init) : _s(_size),
                                                           _data(init)
        {
            assert((sizeof(window_t) * 8) / bitsPerNucl >= _s);
        }

        SmallRollingNuclSeq(size_t _size) : SmallRollingNuclSeq(_size, 0)
        {
        }

        SmallRollingNuclSeq(size_t _size, Dna4SequenceConstIter _it) : SmallRollingNuclSeq(_size, 0)
        {
            for (size_t i = 0; i < _s; ++i)
                _data |= parsing::dna4ToLong(*_it++) << (2 * i);
        }

        SmallRollingNuclSeq(Dna4SequenceConstIter begin, Dna4SequenceConstIter end) : SmallRollingNuclSeq(std::distance(begin, end), 0)
        {
            for (size_t i = 0; i < _s; ++i)
                _data |= parsing::dna4ToLong(*begin++) << (2 * i);
        }

        SmallRollingNuclSeq(const SmallRollingNuclSeq &other) = default;
        SmallRollingNuclSeq(SmallRollingNuclSeq &&other) = default;

        SmallRollingNuclSeq &operator=(const SmallRollingNuclSeq &other) =default;
        SmallRollingNuclSeq &operator=(SmallRollingNuclSeq &&other) =default;

        bool operator==(const SmallRollingNuclSeq &other) const =default;
        bool operator!=(const SmallRollingNuclSeq &other) const =default;

        // rolls back base-pair to most significant position
        inline void rollBack(window_t nt)
        {
            _data = (_data >> bitsPerNucl) | (nt << (bitsPerNucl * (_s - 1)));
        }

        // appends least significant base-pair
        inline SmallRollingNuclSeq operator+(window_t nt) const
        {
            return SmallRollingNuclSeq(_s + 1, (_data << bitsPerNucl) | (nt & 0b11));
        }

        // shift in units of base-pairs
        inline SmallRollingNuclSeq operator<<=(size_t w_)
        {
            _s += w_;
            _data <<= bitsPerNucl * w_;
        }

        // shift in units of base-pairs
        inline SmallRollingNuclSeq operator<<(size_t w_) const
        {
            return SmallRollingNuclSeq(_s + w_, _data << (bitsPerNucl * w_));
        }

        // shift in units of base-pairs
        inline SmallRollingNuclSeq operator>>=(size_t w_)
        {
            assert(w_ <= _s);
            _s -= w_;
            _data >>= bitsPerNucl * w_;
        }

        // shift in units of base-pairs
        inline SmallRollingNuclSeq operator>>(size_t w_) const
        {
            return SmallRollingNuclSeq(_s - w_, _data >> (bitsPerNucl * w_));
        }

        /**
         * Iterates from most-to-least significant byte of small sequence.
         */
        class Iterator
        {
        public:
            // initialise from value and width in bits
            Iterator(uint64_t &v_, uint8_t w_) : _value(v_),
                                                 i(static_cast<int>((w_ - 1) / 8))
            {
                assert(w_ > 0);
            }

            Iterator(const Iterator &rhs) : _value(rhs._value), i(rhs.i) {}

            Iterator(Iterator &&rhs) : _value(rhs._value), i(rhs.i) {}

        protected:
            uint64_t &_value;
            int i;

        public:
            using difference_type = int;
            using value_type = uint8_t;

            inline Iterator &operator=(const Iterator &rhs_)
            {
                _value = rhs_._value;
                i = rhs_.i;
                return *this;
            }

            inline Iterator &operator=(Iterator &&rhs_)
            {
                _value = rhs_._value;
                i = rhs_.i;
                return *this;
            }

            inline uint8_t operator*() const
            {
                assert(i >= 0);
                return static_cast<uint8_t>((_value >> (8 * i)) & 0xFF);
            }

            inline Iterator &operator++()
            {
                --i;
                return *this;
            }

            void operator++(int) { ++*this; }
        };

    protected:
        size_t _s;
        window_t _data;

    public:
        // return iterator to most significant byte
        inline Iterator bytes() { return Iterator(_data, _s * bitsPerNucl); }

        window_t data() const { return _data; }

        size_t size() const { return _s; }

        std::string toString() const;

        inline uint8_t msb() const noexcept { return _data >> (bitsPerNucl * (_s - 1)); }
    };

    class SuffixTable
    {
    public:
        using count_t = uint64_t;
        using CountVector = std::vector<count_t>;
        using CountVectorIter = typename CountVector::iterator;
        using CountVectorConstIter = typename CountVector::const_iterator;
        using CountVectorPointer = typename CountVector::pointer;

        SuffixTable() = default;

        SuffixTable(const SuffixTable &other) = default;

        SuffixTable(SuffixTable &&other) = default;

        SuffixTable(size_t suffixSize_)
            : s(suffixSize_),
              _data((size_t)1u << (2 * s), 0)
        {
        }

    private:
        size_t s;
        CountVector _data;

    public:
        inline void resize(size_t s_)
        {
            s = s_;
            _data.resize((size_t)1u << (2 * s), 0);
        }

        inline SuffixTable &operator=(const SuffixTable &other) = default;

        inline SuffixTable &operator=(SuffixTable &&other) = default;

        inline size_t suffixSize() const noexcept { return s; }

        const CountVector &cdata() const { return _data; }

        inline CountVectorIter begin() { return _data.begin(); }

        inline CountVectorIter end() { return _data.end(); }

        inline CountVectorConstIter begin() const { return _data.begin(); }

        inline CountVectorConstIter end() const { return _data.end(); }

        inline CountVectorConstIter cbegin() const { return _data.cbegin(); }

        inline CountVectorConstIter cend() const { return _data.cend(); }

        void count(Dna4SequenceConstIter it,
                   Dna4SequenceConstIter end);

        inline count_t &operator[](const SmallRollingNuclSeq &sfx) { return _data.operator[](sfx._data); }

        inline const count_t &operator[](uint64_t sfx) const { return _data.operator[](sfx); }

        inline const count_t &operator[](const SmallRollingNuclSeq &sfx) const { return _data.operator[](sfx._data); }

        size_t maxValue() const;

        inline size_t size() const noexcept { return _data.size(); }

        friend SuffixTable &inplaceAdd(SuffixTable &dest, SuffixTable const &rhs);
    };

    SuffixTable &inplaceAdd(SuffixTable &dest, SuffixTable const &rhs);

    template <class Unwind_>
    std::vector<SuffixTable>
    countSuffixes(std::vector<range_t<Unwind_>> const &ranges, Unwind_ Apply, size_t suffixSize, size_t offset = 0u)
    {
        std::vector<SuffixTable> tables(ranges.size());
        oneapi::tbb::parallel_for((size_t)0, ranges.size(), (size_t)1u, [&](size_t i)
            {
                tables[i].resize(suffixSize);
                Apply.forEach(ranges[i], [&](Dna4SequenceConstIter begin, Dna4SequenceConstIter end, bool endIsTerminal)
                    {
                        if (!endIsTerminal)
                            --end;
                        size_t rangeSize = end - begin;
                        if (rangeSize >= offset + suffixSize)
                            tables[i].count(begin + offset, end);
                    });
            });
        return tables;
    }

    SuffixTable accumulate(std::vector<SuffixTable> const &tables);

    class SuffixBV
    {
    protected:
        static constexpr size_t _w = 64;

    public:
        SuffixBV(size_t s);

    protected:
        static size_t offset(size_t s);

        static size_t bitIndex(SmallRollingNuclSeq const &sfx);

    public:
        void insert(SmallRollingNuclSeq const &sfx);

        bool contains(SmallRollingNuclSeq const &sfx) const;

    protected:
        std::vector<copyable_atomic_uint64_t> _data;
    };

    template <class Unwind_>
    SuffixBV detectSuffixes(std::vector<range_t<Unwind_>> const &ranges, Unwind_ Apply, size_t suffixSize)
    {
        SuffixBV bv(suffixSize);
        oneapi::tbb::parallel_for((size_t)0, ranges.size(), (size_t)1u, [&](size_t i)
            {
                Apply.forEachBegin(ranges[i], [&](Dna4SequenceConstIter begin)
                    {
                        SmallRollingNuclSeq sfx(suffixSize, begin);
                        while (sfx.size() > 0)
                        {
                            bv.insert(sfx);
                            sfx >>= 1;
                        }
                        bv.insert(sfx);
                    } );
            } );
        return bv;
    }

} // namespace suffix
