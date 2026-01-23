#pragma once
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <execution>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>
#include <oneapi/tbb.h>
#include <maki/colourEncoder.hpp>
#include <maki/fasta.hpp>
#include <maki/matrix.hpp>
#include <maki/suffix.hpp>
#include <maki/mers.hpp>
#include <maki/utils.hpp>

namespace terminals
{

    /**
     * Subclass on contiguous rolling sequence that tracks how many
     * base-pairs have been inserted so far. This is used to represent
     * terminals at the beginning of contigs.
     */
    struct RollingTerminal : public mers::ContiguousRollingNuclSeq
    {
        RollingTerminal(size_t _length)
            : mers::ContiguousRollingNuclSeq(_length)
            , _size(0)
        {
        }

        template <class... Args>
        void rollBack(Args &&...args)
        {
            mers::ContiguousRollingNuclSeq::rollBack(std::forward<Args>(args)...);
            if (_size < _num_nts)
                ++_size;
        }

        // Returns number of base-pairs currently in terminal (starts with 0).
        inline size_t terminalLength() const noexcept { return _size; }

    public:
        size_t _size;
    };

    class TerminalRange;

    namespace _detail
    {

        void ser(size_t src,
                 uint8_t *dst, uint32_t bytes);

        size_t deser(uint8_t const *src, uint32_t bytes);

    } // namespace _detail

    template <typename T>
    class TerminalBufferRandomAccessIterator : public ndim::StridedIteratorBase<T, TerminalBufferRandomAccessIterator<T>>
    {
        using BaseIterType = ndim::StridedIteratorBase<T, TerminalBufferRandomAccessIterator<T>>;
        using difference_type = std::ptrdiff_t;

    protected:
        uint32_t
            _lengthBytes,
            _seqBytes,
            _valueOffset;
        uint8_t _k_eff;

    public:
        TerminalBufferRandomAccessIterator(uint32_t lengthBytes, uint32_t seqBytes, uint32_t valueOffset, uint8_t kEff) 
            : BaseIterType()
            , _lengthBytes(lengthBytes)
            , _seqBytes(seqBytes)
            , _valueOffset(valueOffset)
            , _k_eff(kEff)
        {
        }

        TerminalBufferRandomAccessIterator(T *p, size_t w, uint32_t lengthBytes, uint32_t seqBytes, uint32_t valueOffset, uint8_t kEff)
            : BaseIterType(p, w)
            , _lengthBytes(lengthBytes)
            , _seqBytes(seqBytes)
            , _valueOffset(valueOffset)
            , _k_eff(kEff)
        {
        }

    public:
        constexpr inline void writeTerminal(RollingTerminal const &terminal) const
            requires(!std::is_const_v<T>)
        {
            std::memcpy(this->_data + _lengthBytes, terminal.data(), _seqBytes);
        }

        constexpr inline void writeSize(size_t v_) const
            requires(!std::is_const_v<T>)
        {
            _detail::ser(v_, this->_data, _lengthBytes);
        }

        constexpr inline void writeKey(RollingTerminal const &tl) const
            requires(!std::is_const_v<T>)
        {
            writeSize(tl.terminalLength());
            writeTerminal(tl);
        }

        constexpr inline void writeKey(Dna4SequenceConstIter in, std::size_t size_) const
            requires(!std::is_const_v<T>)
        {
            writeSize(size_);
            mers::reverseWriteMerTo(in, this->_data + _lengthBytes, _k_eff, size_);
        }

        constexpr inline void writeEdge(T edge_) const
            requires(!std::is_const_v<T>)
        {
            *(this->_data + _valueOffset) = edge_;
        }

        inline uint32_t lengthBytes() const noexcept { return _lengthBytes; }

        inline uint32_t keyBytes() const noexcept { return _seqBytes + _lengthBytes; }

        inline size_t readSize() const { return _detail::deser(this->_data, _lengthBytes); }

        inline T readEdge() const { return *(this->_data + _valueOffset); }

        inline T readValue() const { return readEdge(); }

        inline T *data() { return this->_data; }

        inline const T *cdata() const { return this->_data; }

        // most significant byte of key
        inline T keyMSB() const { return *(this->_data + keyBytes() - 1); }
    };

    /**
     * Values are serialised in orientation (LSB --> MSB):
     *
     *      | --- length --- | --- sequence --- | 0 ... 0 | edge |
     */
    class TerminalBuffer : public mers::BaseKmerBuffer<TerminalBufferRandomAccessIterator>
    {
        friend struct TerminalsLessThan;

    protected:
        using BaseBufferType = BaseKmerBuffer<TerminalBufferRandomAccessIterator>;

    public:
        static constexpr struct autofit_tag_t{} autofit_tag{};

        explicit TerminalBuffer(uint64_t numRecords, uint8_t k_, uint8_t k_eff_, uint32_t lengthBytes)
            : BaseBufferType(numRecords, _valueBytes, k_, k_eff_, lengthBytes)
            , _lengthBytes(lengthBytes)
            , _seqBytes(_key_bytes - _lengthBytes)
        {
        }

        explicit TerminalBuffer(uint64_t numRecords, uint8_t k_, uint8_t k_eff_, autofit_tag_t)
            : TerminalBuffer(numRecords, k_, k_eff_, required_bytes(k_eff_))
        {
        }

        explicit TerminalBuffer(uint64_t numRecords, uint8_t k_, uint32_t lengthBytes)
            : BaseBufferType(numRecords, _valueBytes, k_, k_, lengthBytes)
            , _lengthBytes(lengthBytes)
            , _seqBytes(_key_bytes - _lengthBytes)
        {
        }

        explicit TerminalBuffer(uint64_t numRecords, uint8_t k_, autofit_tag_t)
            : TerminalBuffer(numRecords, k_, required_bytes(k_))
        {
        }

        explicit TerminalBuffer(uint8_t k_,
                                uint8_t k_eff_,
                                uint32_t lengthBytes,
                                std::initializer_list<std::initializer_list<uint8_t>> data)
            : BaseBufferType(_valueBytes, k_, k_eff_, lengthBytes, data)
            , _lengthBytes(lengthBytes)
            , _seqBytes(_key_bytes - _lengthBytes)
        {
        }

        explicit TerminalBuffer(uint8_t k_,
                                uint8_t k_eff_,
                                std::initializer_list<std::initializer_list<uint8_t>> data,
                                autofit_tag_t)
            : TerminalBuffer(k_, k_eff_, required_bytes(k_eff_), data)
        {
        }

        explicit TerminalBuffer(uint8_t k_,
                                uint32_t lengthBytes,
                                std::initializer_list<std::initializer_list<uint8_t>> data)
            : BaseBufferType(_valueBytes, k_, k_, lengthBytes, data)
            , _lengthBytes(lengthBytes)
            , _seqBytes(_key_bytes - _lengthBytes)
        {
        }

        explicit TerminalBuffer(uint8_t k_,
                                std::initializer_list<std::initializer_list<uint8_t>> data,
                                autofit_tag_t)
            : TerminalBuffer(k_, required_bytes(k_), data)
        {
        }

        template <std::input_iterator InputIt>
        TerminalBuffer(const TerminalBuffer& src_, InputIt it_, InputIt end_)
            : TerminalBuffer(std::distance(it_, end_), src_.getK(), src_.getEffK(), src_.lengthBytes())
        {
            tbb::parallel_for((size_t)0, (size_t)_num_records, [&](size_t i) {
                auto __srcspan = *src_.at(*(it_ + i));
                std::copy(__srcspan.begin(), __srcspan.end(), data(i));
            } );
        }

    private:
        uint32_t
            _lengthBytes,
            _seqBytes;
        static constexpr uint32_t _valueBytes = 1u;

    public:
        inline uint32_t lengthBytes() const noexcept { return _lengthBytes; }

        inline static constexpr uint32_t valueBytes() noexcept { return _valueBytes; }

        // Will not reallocate - can only reduce size. Does not modify values (even those now out of reach).
        void shrink(size_t new_size_);

    public:
        inline iterator begin() {
            return _data.begin(_lengthBytes, _seqBytes, valueOffset(), k_eff);
        }

        inline iterator end() {
            return _num_records == _record_capacity
                ? _data.end(_lengthBytes, _seqBytes, valueOffset(), k_eff)
                : _data.begin(_lengthBytes, _seqBytes, valueOffset(), k_eff) + _num_records;
        }

        inline const_iterator constBegin() const noexcept {
            return _data.cbegin(_lengthBytes, _seqBytes, valueOffset(), k_eff);
        }

        inline const_iterator constEnd() const noexcept {
            return _num_records == _record_capacity
                ? _data.cend(_lengthBytes, _seqBytes, valueOffset(), k_eff)
                : _data.cbegin(_lengthBytes, _seqBytes, valueOffset(), k_eff) + _num_records;
        }

        inline const_iterator begin() const noexcept { return constBegin(); }

        inline const_iterator end() const noexcept { return constEnd(); }

        inline iterator at(size_type idx) { return _data.at(idx, _lengthBytes, _seqBytes, valueOffset(), k_eff); }

        inline const_iterator constAt(size_type idx) const { return _data.at(idx, _lengthBytes, _seqBytes, valueOffset(), k_eff); }

        inline const_iterator at(size_type idx) const { return _data.at(idx, _lengthBytes, _seqBytes, valueOffset(), k_eff); }

    public:
        iterator insert(iterator it,
                        Dna4SequenceConstIter begin,
                        Dna4SequenceConstIter end,
                        bool endIsTerminal);

        iterator insert(iterator it,
                        Dna4SequenceConstIter begin,
                        Dna4SequenceConstIter end,
                        bool endIsTerminal,
                        suffix::SmallRollingNuclSeq key);

        
        void sort(TerminalBuffer *temp, uint32_t threads_);
        void sort(uint32_t threads_);
        void unique();
        TerminalBuffer OOPsort() const;
        TerminalRange asRange() const;
    };

    using TerminalConstIter = typename TerminalBuffer::const_iterator;

    class TerminalRange
    {
    public:
        using iterator = typename TerminalBuffer::const_iterator;
        using const_iterator = typename TerminalBuffer::const_iterator;

    public:
        TerminalRange(uint8_t k_eff_, long lengthBytes_, TerminalConstIter begin_, TerminalConstIter end_)
            : _k_eff(k_eff_)
            , _lengthBytes(lengthBytes_)
            , _begin(begin_)
            , _end(end_)
        {
        }

        /** Split constructor.
         */
        TerminalRange(TerminalRange &rhs, TerminalConstIter pivot)
            : _k_eff(rhs._k_eff)
            , _lengthBytes(rhs._lengthBytes)
            , _begin(rhs._begin)
            , _end(pivot)
        {
            rhs._begin = pivot;
        }

    protected:
        uint8_t _k_eff;

        long _lengthBytes;

        TerminalConstIter
            _begin,
            _end;

    public:
        inline TerminalConstIter constBegin() const noexcept { return _begin; }

        inline TerminalConstIter constEnd() const noexcept { return _end; }

        inline TerminalConstIter begin() const noexcept { return constBegin(); }

        inline TerminalConstIter end() const noexcept { return constEnd(); }

        inline uint8_t getEffK() const noexcept { return _k_eff; }

        inline uint32_t lengthBytes() const noexcept { return _lengthBytes; }

        inline size_t size() const noexcept { return _end - _begin; }

        inline bool empty() const noexcept { return _begin == _end; }

    public:
        inline void stepForward(size_t di_) { _begin += di_; }

    public:
        /**
         * Find lower bound `target` in the range `[begin, end)`.
         */
        TerminalConstIter lowerBound(ndim::Span<uint8_t *>, uint8_t) const;
        TerminalConstIter lowerBound(ndim::Span<uint8_t *>) const;

        /**
         * Find upper bound `target` in the range `[begin, end)`.
         */
        TerminalConstIter upperBound(ndim::Span<uint8_t *>, uint8_t) const;
        TerminalConstIter upperBound(ndim::Span<uint8_t *>) const;

    public:
        TerminalRange endsWith(suffix::SmallRollingNuclSeq) const;
        TerminalRange retrieve(suffix::SmallRollingNuclSeq) const;
    };

    // Operators

    struct TerminalsLessThan
    {
    public:
        TerminalsLessThan(const TerminalBuffer &buffer_)
            : _buffer(buffer_)
            , _keyWidth(buffer_.keyBytes())
            , _valueOffset(buffer_.valueOffset())
        {
        }

    protected:
        TerminalBuffer const &_buffer;
        uint32_t
            _keyWidth,
            _valueOffset;

    public:
        bool operator()(size_t i_, size_t j_) const;

    protected:
        static int _rng_cmp(const uint8_t *l_, const uint8_t *r_, uint32_t width);
    };

    struct TerminalIsSuffix
    {
        TerminalIsSuffix(uint8_t k_, uint8_t s_, long offset_ = 0)
            : _trg_width(_mkimem_details::key_size(k_) + offset_)
            , _sfx_shift((k_ - s_) % 4)
            , _sfx_width((s_ + _sfx_shift + 3) / 4)
            , _offset(offset_)
        {
            assert(s_ <= k_);
        }

    protected:
        long
            _trg_width; // number of bytes per search target (k-mer)
        size_t
            _sfx_shift; // shift each input suffix to align with k-mer byte orientation
        long
            _sfx_width, // number of bytes per (shifted) search query (terminal)
            _offset;

    public:
        bool operator()(uint8_t const *trg_, uint64_t sfx_) const;

        bool operator()(TerminalConstIter it, suffix::SmallRollingNuclSeq const &sfx) const;
    };

    /**
     * Compare two sequences whose bytes first need to be aligned (i.e. compare suffix with a terminal)
     */
    template <bool mask_lsb>
    struct UnalignedBytesLessThan
    {
        /**
         * `k_` defines the template orientation, and we need to modify `s_`
         * to fit the template.
         */
        UnalignedBytesLessThan(uint8_t k_, uint8_t s_, uint8_t offset_ = 0)
            : _trg_width(_mkimem_details::key_size(k_) + offset_)
            , _sfx_shift((k_ - s_) % 4)
            , _sfx_width((s_ + _sfx_shift + 3) / 4)
            , _trg_lsb_msk((0xFFu << (2 * _sfx_shift)) & 0xFFu)
        {
            assert(s_ <= k_);
        }
 
    protected:
        long _trg_width;        // number of bytes per search target (k-mer)
        size_t _sfx_shift;      // shift each input suffix to align with k-mer byte orientation
        long _sfx_width;        // number of bytes per (shifted) search query (terminal)
        uint8_t _trg_lsb_msk;   // apply mask to least-significant-byte of search target

    public:
        template <typename Op>
        bool operator()(ndim::Span<uint8_t const*> trg_, uint64_t sfx_, Op cmp) const
        {
            uint8_t const *trg_aln = trg_.data() + _trg_width - _sfx_width;
            sfx_ <<= 2 * _sfx_shift;

            for (long i = _sfx_width - 1; i > -1; --i)
            {
                uint8_t
                    trg_v = *(trg_aln + i),
                    sfx_v = (sfx_ >> (8 * i)) & 0xFFu;

                if constexpr (mask_lsb)
                {
                    if (i == 0)
                        trg_v &= _trg_lsb_msk;
                }

                if (cmp(trg_v, sfx_v))
                    return true;
                else if (cmp(sfx_v, trg_v))
                    return false;
            }

            return false;
        }

        inline bool operator()(ndim::Span<const uint8_t *> trg_, uint64_t sfx_) const {
            return operator()(trg_, sfx_, std::less{});
        }

        inline bool operator()(uint64_t sfx_, ndim::Span<const uint8_t *> trg_) const {
            return operator()(trg_, sfx_, std::greater{});
        }
    };

    enum mask_parity
    {
        LHS,
        RHS
    };

    /**
     * Compare two sequences whose bytes are aligned, according to a given k-mer size.
     * Used for comparing k-mers and terminals, where the terminal size may be greater
     * than the k-mer size (and so its most significant bytes should be masked). The
     * constructor parameter `k_` should be the smaller.
     */
    template <mask_parity P>
    struct MaskedAlignedBytesLessThan
    {
        MaskedAlignedBytesLessThan(uint8_t k_, long offset_ = 0)
            : _width(_mkimem_details::key_size(k_))
            , _msb_mask(k_ % 4 == 0 ? 0xFFu : (0xFFu >> (2 * (4 - (k_ % 4)))))
            , _offset(offset_)
        {
        }

    protected:
        size_t _width;
        uint8_t _msb_mask;
        long _offset;

    public:
        /**
         * Compare bytes of two arrays, where `l_` and `r_` point to the least significant
         * bytes. A mask is applied to the most significant byte, the target of which is
         * defined by the template parameter `P`.
         */
        bool operator()(const uint8_t *l_, const uint8_t *r_) const
        {
            uint8_t msk = _msb_mask;

            for (const uint8_t
                    *lhs = l_ + _width - 1,
                    *rhs = r_ + _width - 1;
                lhs != l_ - 1;
                (void)--lhs, (void)--rhs)
            {
                uint8_t
                    l_val = *lhs,
                    r_val = *rhs;

                if constexpr (P == LHS)
                    l_val &= msk;
                else
                    r_val &= msk;

                if (l_val < r_val)
                    return true;
                else if (l_val > r_val)
                    return false;

                msk = 0xFFu;
            }
            return false;
        }

        template <typename X, typename Y>
        inline bool operator()(ndim::Span<X *> lhs_, ndim::Span<Y *> rhs_) const {
            if constexpr (P == LHS)
                return operator()(lhs_.data() + _offset, rhs_.data()          );
            else
                return operator()(lhs_.data()          , rhs_.data() + _offset);
        }
    };

    template <typename Unwind_>
    TerminalBuffer extractTerminals(std::vector<range_t<Unwind_>> const &data, Unwind_ unwind, uint8_t k_)
    {
        auto blocks = accumulatePartials(data, [&](range_t<Unwind_> const &rng)->size_t
            { return rng.numTerminals(k_); });

        TerminalBuffer terminals(blocks.back(), k_, TerminalBuffer::autofit_tag);
        terminals.insertTerminals(data, unwind, blocks);
        
        return terminals;
    }

    template <typename Unwind_>
    inline TerminalBuffer extractTerminalsAndSort(std::vector<range_t<Unwind_>> const &data, Unwind_ unwind, uint8_t k_)
    {
        return extractTerminals(data, unwind, k_).OOPsort();
    }

    using mers::KmerDiff;
    using mers::KmerDiffClass;
    using mers::KmerOverlapVector;

    struct TerminalDiff : public KmerDiff
    {
        TerminalDiff(uint32_t lengthBytes, uint8_t k_eff)
            : KmerDiff(_mkimem_details::key_size(k_eff))
            , _lengthBytes(lengthBytes)
            , _k_eff(k_eff)
        {
        }

        // `_lhs` and `_rhs` point to the least significant bytes of each k-mer
        KmerDiffClass operator()(const uint8_t *_lhs, const uint8_t *_rhs) const;

    protected:
        uint32_t _lengthBytes;
        size_t _k_eff;
    };

} // namespace terminals

mers::KmerOverlapVector adjacentDifference(terminals::TerminalBuffer &);

mers::KmerOverlapVector adjacentDifference(terminals::TerminalRange);

template <mers::KmerDiffClass Xchr>
auto indexedAdjacentDifference(terminals::TerminalBuffer &buffer_) {
    return mers::indexedAdjacentDifference<Xchr>(
        buffer_.constBegin(),
        buffer_.constEnd(),
        terminals::TerminalDiff{ buffer_.lengthBytes(), buffer_.getEffK() });
}

template <mers::KmerDiffClass Xchr>
auto indexedAdjacentDifference(terminals::TerminalRange buffer_) {
    return mers::indexedAdjacentDifference<Xchr>(
        buffer_.constBegin(),
        buffer_.constEnd(),
        terminals::TerminalDiff{ buffer_.lengthBytes(), buffer_.getEffK() });
}
