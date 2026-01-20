#pragma once
#include <atomic>
#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_pipeline.h>
#include <sdsl/vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>
#include <maki/matrix.hpp>
#include <maki/mers.hpp>
#include <maki/terminals.hpp>
#include <maki/vectorbuffer.hpp>

using oneapi::tbb::filter_mode::parallel;
using oneapi::tbb::filter_mode::serial_in_order;
using mers::KmerDiffClass::BW_0_K;
using mers::KmerDiffClass::IS_0;
using mers::KmerDiffClass::IS_K;

namespace flush
{
    using KmerOverlapIter = mers::KmerOverlapVector::iterator;
    using TerminalConstIter = terminals::TerminalBuffer::const_iterator;

    using IndexVector = std::vector<size_t>;
    using IndexVectorConstIter = IndexVector::const_iterator;
    using AtomicCounter = std::array<std::atomic_size_t, 5>;

    namespace detail
    {

        template <typename Iter_>
        concept Sized = requires(Iter_ it)
        {
            { it.readSize() } -> std::same_as<size_t>;
        };

        template <typename Iter_>
        concept Unsized = !Sized<Iter_>;

        template <typename Cursor_>
        concept Coloured = requires(Cursor_ crs)
        {
            { crs.flushColours() } -> std::same_as<void>;
        };

    } // namespace detail

    template <typename Cursor_, typename Iter_>
    void insert(Cursor_ &crs, Iter_ it, Iter_ end, KmerOverlapIter dv, uint8_t msb) {
        while (it != end) {
            if (*dv == BW_0_K) {
                crs.setToCurrent();
            }
            
            do {
                crs.emplace(it.readValue());
                ++it;
                ++dv;
            } while ((it != end) && (*dv == IS_0));
            
            crs.writeNode(msb);
        }
    }

    template <typename Cursor_, typename Iter_>
    void insert(Cursor_ &crs, Iter_ it, Iter_ end, KmerOverlapIter dv) {
        uint8_t _msb;

        while (it != end) {
            if (*dv == BW_0_K) {
                crs.setToCurrent();
            }
            
            if constexpr (detail::Sized<Iter_>) {
                _msb = (it.readSize() == 0)
                    ? 0
                    : mers::dna4ToDna5(it.keyMSB() >> (2 * ((crs.k() - 1) % 4)));
            } else {
                _msb = mers::dna4ToDna5(it.keyMSB() >> (2 * ((crs.k() - 1) % 4)));
            }
            
            do {
                crs.emplace(it.readValue());
                ++it;
                ++dv;
            } while ((it != end) && (*dv == IS_0));
            
            crs.writeNode(_msb);
        }
    }

    template <typename Cursor_, typename Iter_, class ...Args>
    void interleave(Cursor_ &crs,
                    /* SORTED */ Iter_ km_it_,
                    /* SORTED */ Iter_ km_end_,
                    KmerOverlapIter bk_,
                    /* SORTED */ TerminalConstIter tm_it_,
                    /* SORTED */ TerminalConstIter tm_end_,
                    KmerOverlapIter bt_,
                    uint8_t k_,
                    uint8_t kEff_,
                    Args &&...args)
    {
        size_t
            prevTerminalSize = 0,
            kMinus1 = k_ - 1;
        assert(kMinus1 > prevTerminalSize);

        mers::MaskedBytesDiff Diff(kEff_);
        terminals::MaskedAlignedBytesLessThan<terminals::RHS> LessThan(kEff_, tm_it_.lengthBytes());

        while ((km_it_ != km_end_) || (tm_it_ != tm_end_))
        {
            // insert terminals (they always start a new block)

            auto tm_pivot = km_it_ == km_end_
                ? tm_end_
                : std::upper_bound(tm_it_, tm_end_, *km_it_, LessThan);

            if (tm_it_ != tm_pivot)
            {
                insert(crs, tm_it_, tm_pivot, bt_);

                size_t insertedTerminals = tm_pivot - tm_it_;
                assert(insertedTerminals > 0);
                bt_ += insertedTerminals;

                auto tm_penul = tm_it_ + insertedTerminals - 1;
                prevTerminalSize = tm_penul.readSize();
                tm_it_ = tm_pivot;
            }

            // insert k-mers

            auto km_pivot = tm_it_ == tm_end_
                ? km_end_
                : std::lower_bound(km_it_, km_end_, *tm_it_, LessThan);

            if (km_it_ != km_pivot)
            {
                // B correction

                if (prevTerminalSize == kMinus1)
                {
                    // if their bytes are equal and the terminal size is `k-1`, then the only difference is the final `$` char

                    auto dtk = Diff(std::to_address(/* guaranteed */ tm_it_ - 1) + tm_it_.lengthBytes(),
                                    std::to_address(km_it_));
                    *bk_ = dtk == IS_0 ? IS_K : dtk;
                }

                // insert
                insert(crs, km_it_, km_pivot, bk_, args...);

                size_t insertedKmers = km_pivot - km_it_;
                assert(insertedKmers > 0);
                bk_ += insertedKmers;

                km_it_ = km_pivot;
            }
        }
    }

    namespace ranges
    {

        namespace detail {

            size_t calculateLeftChunks(size_t width, size_t chunksize, oneapi::tbb::proportional_split p);

        }

        namespace counting {

            template <typename Iter>
            struct BlockRange {
                Iter begin, end;
                KmerOverlapIter b;
                size_t chunksize;

            public:
                inline size_t size() const
                { return end - begin; }

                inline bool empty() const
                { return begin == end; }

                inline bool is_divisible() const
                { return findNextChunkIndex(1) < size(); }

            private:
                size_t findNextChunkIndex(size_t chunks) const {                
                    size_t limit = end - begin, 
                           rsize = std::min(limit, chunksize * chunks);
                    KmerOverlapIter bp = b + rsize;

                    while ((rsize < limit) && (*bp == IS_0))
                    {
                        ++rsize;
                        ++bp;
                    }

                    return rsize;
                }
            
            public:
                BlockRange(Iter begin_, Iter end_, KmerOverlapIter b_, size_t chunksize_)
                    : begin(begin_)
                    , end(end_)
                    , b(b_)
                    , chunksize(chunksize_)
                {}

                BlockRange(BlockRange const&) =default;

                BlockRange(BlockRange &rhs, size_t chunks) : BlockRange(rhs) {
                    size_t lhsSize = findNextChunkIndex(chunks);
                    rhs.begin += lhsSize;
                    rhs.b     += lhsSize;
                    end = rhs.begin;
                }

                BlockRange(BlockRange &rhs, oneapi::tbb::proportional_split p)
                    : BlockRange(rhs, detail::calculateLeftChunks(rhs.size(), rhs.chunksize, p))
                {}

                BlockRange(BlockRange &rhs, oneapi::tbb::split)
                    : BlockRange(
                        rhs,
                        detail::calculateLeftChunks(
                            rhs.size(),
                            rhs.chunksize,
                            oneapi::tbb::proportional_split(1, 1)))
                {}
            };

            struct CountUnique {
                size_t count;

            protected:
                static constexpr uint8_t terminalBit = 1u << mers::KmerBuffer::terminalEdge;

            public:
                CountUnique() : count(0) {}

                CountUnique([[maybe_unused]] CountUnique &rhs,
                            [[maybe_unused]] oneapi::tbb::split)
                    : count(0)
                {
                }

                template <typename Iter>
                void operator()(Iter it, Iter end, KmerOverlapIter b) {
                    size_t localCount = 0;

                    uint8_t edgeCount;
                    while (it != end)
                    {
                        // collect edges of this node

                        edgeCount = 0;
                        do
                        {
                            edgeCount |= 1u << it.readEdge();
                            ++it;
                            ++b;
                        } while ((it != end) && (*b == IS_0));
                        
                        // if `$` appears with other edges, it is not counted

                        if ((edgeCount & terminalBit) && (edgeCount & ~terminalBit))
                            edgeCount ^= terminalBit;

                        localCount += std::popcount(edgeCount);
                    }

                    count += localCount;
                }

                template <typename Iter>
                inline void operator()(BlockRange<Iter> const &r) {
                    operator()(r.begin, r.end, r.b);
                }

                inline void join(CountUnique &rhs) { count += rhs.count; }
            };

            template <typename Iter>
            size_t uniqueEdges(Iter begin_, Iter end_, KmerOverlapIter b_, size_t chunksize_) {
                CountUnique Sum;
                oneapi::tbb::parallel_reduce(BlockRange<Iter>(begin_, end_, b_, chunksize_), Sum);
                return Sum.count;
            }

        }

        namespace detail
        {

            template <typename Iter_>
            class ForwardRange {
            protected:
                Iter_ it;
                size_t beginIndex, endIndex;

            public:
                ForwardRange(Iter_ it_, size_t begin_, size_t end_)
                    : it(it_)
                    , beginIndex(begin_)
                    , endIndex(end_)
                {
                }

                ForwardRange(Iter_ begin, Iter_ end)
                    : ForwardRange(begin, 0, end - begin)
                {}

                inline Iter_ begin() const noexcept { return it + beginIndex; }
                
                inline Iter_ end() const noexcept { return it + endIndex; }
                
                inline size_t size() const noexcept { return endIndex - beginIndex; }

                inline bool empty() const noexcept { return endIndex == beginIndex; }

                inline void stepForward(size_t di_) { beginIndex += di_; }
            };

            /**
             * Monotone increasing range. Step between adjacent elements stored in `_b`.
             */
            template <typename Iter_>
            class MonotoneRange : public ForwardRange<Iter_>
            {
            protected:
                KmerOverlapIter _b;

            public:
                MonotoneRange(Iter_ it, size_t begin_, size_t end_, KmerOverlapIter b_)
                    : ForwardRange<Iter_>::ForwardRange(it, begin_, end_)
                    , _b(b_)
                {
                }

                MonotoneRange(Iter_ begin, Iter_ end, KmerOverlapIter b_)
                    : ForwardRange<Iter_>::ForwardRange(begin, end)
                    , _b(b_)
                {
                }

                inline KmerOverlapIter B() const noexcept {
                    return _b + this->beginIndex;
                }
            };

            template <typename RankType, typename SelectType>
            class Granular {
            protected:
                size_t              _chunksize;
                RankType const      &rank;
                SelectType const    &select;

            public:
                Granular(size_t chunksize, RankType const &rank_, SelectType const &select_)
                    : _chunksize(chunksize)
                    , rank(rank_)
                    , select(select_)
                {
                }

                Granular(Granular const&) =default;

            public:
                inline size_t chunksize() const noexcept
                { return _chunksize; }

                /**
                 * Number of whole blocks contained in the range.
                 */
                inline size_t blocks(size_t begin, size_t end) const {
                    return rank(end) - rank(begin);
                }

                /**
                 * Number of whole chunks contained in the range.
                 */
                inline size_t chunks(size_t begin, size_t end) const {
                    return blocks(begin, end) / _chunksize;
                }

                /**
                 * Return true if the range contains at least 1 chunk.
                 */
                inline bool is_divisible(size_t begin, size_t end) const {
                    return chunks(begin, end) >= 1;
                }

                /**
                 * Return position corresponding to `chunks` chunks after the beginning, but not beyond `limit`.
                 */
                size_t findChunkPosition(size_t begin, size_t end, size_t chunks) const {
                    size_t requestedRank = rank(begin + 1) + (chunks * _chunksize),
                           endRank = rank(end);
                    return requestedRank <= endRank ? select(requestedRank) : end;
                }
            };

            static constexpr struct pivot_tag_t{} pivot_tag{};
            static constexpr struct chunk_tag_t{} chunk_tag{};

            template <typename Iter_, typename RankType, typename SelectType>
            struct DivisibleRange : public MonotoneRange<Iter_>,
                                    public Granular<RankType,SelectType>
            {
                DivisibleRange(Iter_ it, size_t begin_, size_t end_, KmerOverlapIter b_,
                               size_t chunksize_, RankType const &rank_, SelectType const &select_)
                    : MonotoneRange<Iter_>::MonotoneRange(it, begin_, end_, b_)
                    , Granular<RankType,SelectType>::Granular(chunksize_, rank_, select_)
                {
                }

                DivisibleRange(Iter_ begin, Iter_ end, KmerOverlapIter b_,
                               size_t chunksize_, RankType const &rank_, SelectType const &select_)
                    : MonotoneRange<Iter_>::MonotoneRange(begin, end, b_)
                    , Granular<RankType,SelectType>::Granular(chunksize_, rank_, select_)
                {
                }

                DivisibleRange(DivisibleRange const&) =default;

                inline bool is_divisible() const {
                    return Granular<RankType,SelectType>::is_divisible(this->beginIndex, this->endIndex);
                }

                inline size_t blocks() const {
                    return Granular<RankType,SelectType>::blocks(this->beginIndex, this->endIndex);
                }

                inline size_t chunks() const {
                    return Granular<RankType,SelectType>::chunks(this->beginIndex, this->endIndex);
                }

                DivisibleRange(DivisibleRange &rhs, size_t rhsBegin, pivot_tag_t)
                    : MonotoneRange<Iter_>::MonotoneRange(rhs.it, rhs.beginIndex, rhsBegin, rhs._b)
                    , Granular<RankType,SelectType>::Granular(rhs._chunksize, rhs.rank, rhs.select)
                {
                    rhs.beginIndex = rhsBegin;
                }

                DivisibleRange(DivisibleRange &rhs, size_t chunks, chunk_tag_t)
                    : DivisibleRange(
                        rhs,
                        rhs.findChunkPosition(rhs.beginIndex, rhs.endIndex, chunks),
                        pivot_tag)
                {
                }

                DivisibleRange(DivisibleRange &rhs, oneapi::tbb::proportional_split p)
                    : DivisibleRange(
                        rhs,
                        calculateLeftChunks(rhs.blocks(), rhs.chunksize(), p),
                        chunk_tag)
                {
                }

                DivisibleRange(DivisibleRange &rhs, oneapi::tbb::split)
                    : DivisibleRange(
                        rhs,
                        calculateLeftChunks(rhs.blocks(), rhs.chunksize(), oneapi::tbb::proportional_split(1, 1)),
                        chunk_tag)
                {
                }
            };

            class EdgeOutput {
            protected:
                size_t _edge;

            public:
                EdgeOutput(size_t edge_) : _edge(edge_) {}

                inline size_t edgePosition() const noexcept { return _edge; }

                inline void moveEdgeForward(size_t de_) { _edge += de_; }
            };

            template <typename Iter_, typename RankType, typename SelectType>
            struct IORange : public DivisibleRange<Iter_,RankType,SelectType>,
                             public EdgeOutput
            {
                IORange(Iter_ it, size_t begin_, size_t end_, KmerOverlapIter b_,
                        size_t chunksize_, RankType const &rank_, SelectType const &select_,
                        size_t edge_)
                    : DivisibleRange<Iter_,RankType,SelectType>(
                        it, begin_, end_, b_,
                        chunksize_, rank_, select_)
                    , EdgeOutput(edge_)
                {
                }

                template <class ...Args>
                IORange(IORange &rhs, Args&& ...args)
                    : DivisibleRange<Iter_,RankType,SelectType>(rhs, std::forward<Args>(args)...)
                    , EdgeOutput(rhs.edgePosition())
                {
                    rhs.moveEdgeForward(
                        counting::uniqueEdges<Iter_>(
                            this->begin(),
                            this->end(),
                            this->B(),
                            1000000));
                }
            };

            template <typename Iter, typename RankType, typename SelectType>
            class InterleaveRanges : public EdgeOutput
            {
            protected:
                DivisibleRange<Iter,RankType,SelectType>
                                         _rng;
                terminals::TerminalRange _tml;
                KmerOverlapIter          _tmlB;
                uint8_t                  _k, _kEff;

            public:
                inline uint8_t k()    const noexcept { return _k;    }
                inline uint8_t kEff() const noexcept { return _kEff; }

                // Data attributes

                inline auto   dataBegin() const noexcept { return _rng.begin(); }
                inline auto   dataEnd()   const noexcept { return _rng.end();   }
                inline auto   dataB()     const noexcept { return _rng.B();     }
                inline size_t dataSize()  const noexcept { return _rng.size();  }
                inline bool   dataEmpty() const noexcept { return _rng.empty(); }

                // Terminal attributes

                inline auto   terminalBegin() const noexcept { return _tml.constBegin(); }
                inline auto   terminalEnd() const noexcept   { return _tml.constEnd();   }
                inline auto   terminalB() const noexcept     { return _tmlB;             }
                inline size_t terminalSize() const noexcept  { return _tml.size();       }
                inline bool   terminalEmpty() const noexcept { return _tml.empty();      }

                // combined attributes

                inline size_t size() const noexcept  { return dataSize() + terminalSize();    }
                inline bool   empty() const noexcept { return dataEmpty() && terminalEmpty(); }

            public:
                InterleaveRanges(Iter it, size_t begin_, size_t end_, KmerOverlapIter b_,
                                 size_t chunksize_, RankType const &rank_, SelectType const &select_,
                                 terminals::TerminalRange tmnls_, KmerOverlapIter tmlB_,
                                 uint8_t k_, uint8_t kEff_,
                                 size_t edge_)
                    : EdgeOutput(edge_)
                    , _rng(
                        it, begin_, end_, b_,
                        chunksize_, rank_, select_)
                    , _tml(tmnls_)
                    , _tmlB(tmlB_)
                    , _k(k_)
                    , _kEff(kEff_)
                {
                }

                InterleaveRanges(Iter begin_, Iter end_, KmerOverlapIter b_,
                                 size_t chunksize_, RankType const &rank_, SelectType const &select_,
                                 terminals::TerminalRange tmnls_, KmerOverlapIter tmlB_,
                                 uint8_t k_, uint8_t kEff_,
                                 size_t edge_)
                    : EdgeOutput(edge_)
                    , _rng(
                        begin_, end_, b_,
                        chunksize_, rank_, select_)
                    , _tml(tmnls_)
                    , _tmlB(tmlB_)
                    , _k(k_)
                    , _kEff(kEff_)
                {
                }

                /** Split constructor.
                 */
                template <class ...Args>
                InterleaveRanges(InterleaveRanges &rhs, Args&& ...args) 
                    : EdgeOutput(rhs.edgePosition())
                    , _rng(rhs._rng, std::forward<Args>(args)...)
                    , _tml(
                        rhs._tml,
                        rhs.dataEmpty()
                            ? rhs.terminalEnd()
                            : rhs._tml.upperBound(*(dataEnd() - 1), rhs.kEff())
                        )
                    , _tmlB(rhs._tmlB)
                    , _k(rhs.k())
                    , _kEff(rhs.kEff())
                {
                    rhs._tmlB += terminalSize();
                    rhs.moveEdgeForward(
                        counting::uniqueEdges(dataBegin(), dataEnd(), dataB(), 1000000)
                            + counting::uniqueEdges(terminalBegin(), terminalEnd(), terminalB(), 1000000)
                        );
                }
            };

            struct TbbFree {
                template <typename T>
                inline constexpr void operator()(T *obj) const {
                    oneapi::tbb::tbb_allocator<char>().deallocate(reinterpret_cast<char*>(obj), sizeof(T));
                } 
            };

            template <typename T, class ...Args>
            T* alloc(Args&& ...args)
            {
                char *raw = oneapi::tbb::tbb_allocator<char>().allocate(sizeof(T));
                
                if (raw == nullptr)
                {
                    throw std::runtime_error(std::string("Could not allocate type ") + std::string(typeid(T).name()));
                }
                
                // in-place construction with pre-allocated memory
                T* obj = new (raw) T(std::forward<Args>(args)...);
                return obj;
            }

            template <typename Range_>
            struct SplitRange {
                Range_ *_r;

                SplitRange(Range_ *r_) : _r(r_) {}

                inline Range_* operator()(oneapi::tbb::flow_control &fc) const {
                    if (_r->empty())
                    {
                        fc.stop();
                        return nullptr;
                    }
                    else return nextRange();
                }

                inline Range_* nextRange() const {
                    Range_ *p = alloc<Range_>(*_r, size_t(1), chunk_tag);
                    return p;
                }
            };

            template <typename Cursor_, class ...Args>
            class _FlushOperator {
            protected:
                Cursor_             *_c;
                AtomicCounter       *_F,
                                    *_C;
                TbbFree             rel;
                std::tuple<Args...> _args;

            public:
                _FlushOperator(Cursor_ *c_, AtomicCounter *F_, AtomicCounter *C_, Args& ...args)
                    : _c(c_)
                    , _F(F_)
                    , _C(C_)
                    , rel()
                    , _args(args...)
                {
                }

                _FlushOperator(_FlushOperator const&) =default;
                _FlushOperator(_FlushOperator&&) =default;
                _FlushOperator& operator=(_FlushOperator const&) =default;
                _FlushOperator& operator=(_FlushOperator&&) =default;
            
            protected:
                template <class ...CArgs>
                inline Cursor_ _moveCursorTo(CArgs&& ...cargs) const
                { return _c->moveTo(std::forward<CArgs>(cargs)...); }

                inline void _flushCounters(Cursor_ const &crs) const {
                    (*_F)[0] += crs.F[0];
                    (*_F)[1] += crs.F[1];
                    (*_F)[2] += crs.F[2];
                    (*_F)[3] += crs.F[3];
                    (*_F)[4] += crs.F[4];
                    (*_C)[0] += crs.C[0];
                    (*_C)[1] += crs.C[1];
                    (*_C)[2] += crs.C[2];
                    (*_C)[3] += crs.C[3];
                    (*_C)[4] += crs.C[4];
                }

                template <typename Iter_, typename RankType, typename SelectType, std::size_t... I>
                void _flushRange(Cursor_ &crs, InterleaveRanges<Iter_,RankType,SelectType> *rng, std::index_sequence<I...>) const {
                    interleave(crs,
                               rng->dataBegin(),
                               rng->dataEnd(),
                               rng->dataB(),
                               rng->terminalBegin(),
                               rng->terminalEnd(),
                               rng->terminalB(),
                               rng->k(),
                               rng->kEff(),
                               std::get<I>(_args)...);
                }

                template <typename Iter_, typename RankType, typename SelectType, std::size_t... I>
                void _flushRange(Cursor_ &crs, IORange<Iter_,RankType,SelectType> *rng, std::index_sequence<I...>) const {
                    insert(crs,
                           rng->begin(),
                           rng->end(),
                           rng->B(),
                           std::get<I>(_args)...);
                }

                template <typename Range_>
                void _flush(Cursor_ &crs, Range_ *rng) const {
                    _flushRange(crs, rng, std::make_index_sequence<sizeof...(Args)>());
                    _flushCounters(crs);
                    rel(rng);
                }
            };

            template <typename Cursor_, class ...Args>
            struct Flush : public _FlushOperator<Cursor_,Args...>
            {
                Flush(Cursor_ *c_, AtomicCounter *F_, AtomicCounter *C_, Args& ...args)
                    : _FlushOperator<Cursor_,Args...>(c_, F_, C_, args...)
                {
                }

                Flush(Flush const&) =default;
                Flush(Flush&&) =default;
                Flush& operator=(Flush const&) =default;
                Flush& operator=(Flush&&) =default;

                template <typename Range_>
                void operator()(Range_ *p) const {
                    auto crs = this->_moveCursorTo(p->edgePosition());
                    this->_flush(crs, p);
                }
            };

            template <typename Cursor_, class ...Args>
                requires flush::detail::Coloured<Cursor_>
            struct Flush<Cursor_,Args...> : public _FlushOperator<Cursor_,Args...>
            {
            protected:
                std::shared_ptr<tbb::enumerable_thread_specific<ZstdCompressor>> cxp;
 
            public:
                Flush(Cursor_ *c_, AtomicCounter *F_, AtomicCounter *C_, Args& ...args)
                    : _FlushOperator<Cursor_,Args...>(c_, F_, C_, args...),
                      cxp(std::make_shared<tbb::enumerable_thread_specific<ZstdCompressor>>())
                {
                }

                Flush(Flush const&) =default;
                Flush(Flush&&) =default;
                Flush& operator=(Flush const&) =default;
                Flush& operator=(Flush&&) =default;

                template <typename Range_>
                void operator()(Range_ *p) const {
                    auto crs = this->_moveCursorTo(p->edgePosition(), p->size(), &cxp->local());
                    this->_flush(crs, p);
                }
            };

        } // namespace detail

    } // namespace ranges

    template <typename Cursor, typename OutputRange, class ...Args>
    void flushOutputRange(Cursor &c_, OutputRange *Range, Args ...args) {
        AtomicCounter
            _F = {},
            _C = {};

        ranges::detail::SplitRange            __Split(Range);
        ranges::detail::Flush<Cursor,Args...> __Flush(&c_, &_F, &_C, args...);

        oneapi::tbb::parallel_pipeline(
            PARALLEL_FLUSH_TOKENS,
            oneapi::tbb::make_filter<void,         OutputRange*>(serial_in_order, __Split),
            oneapi::tbb::make_filter<OutputRange*, void        >(parallel,        __Flush));

        // increment output cursor

        c_.F[0] += _F[0];
        c_.F[1] += _F[1];
        c_.F[2] += _F[2];
        c_.F[3] += _F[3];
        c_.F[4] += _F[4];
        c_.C[0] += _C[0];
        c_.C[1] += _C[1];
        c_.C[2] += _C[2];
        c_.C[3] += _C[3];
        c_.C[4] += _C[4];

        c_ += _F[0] + _F[1] + _F[2] + _F[3] + _F[4];
    }

    template <typename Cursor, typename InputRange, class ...Args>
    void flush(Cursor &c_, InputRange &ir, size_t grainmod, Args&& ...args) {
        LOG(INFO) << "calculating adjacent overlaps";
        std::pair<sdsl::int_vector<2>,sdsl::bit_vector> B = indexedAdjacentDifference<BW_0_K>(ir);

        // initialise rank-select support
        LOG(INFO) << "initialising rank-select support";
        sdsl::rank_support_v5<1,1>    rank;
        sdsl::select_support_mcl<1,1> select;
        oneapi::tbb::parallel_invoke(
            [&]{ sdsl::util::init_support(rank,   &std::get<1>(B)); },
            [&]{ sdsl::util::init_support(select, &std::get<1>(B)); }
        );

        size_t grainsize = std::max(
            (size_t)1,
            (rank(ir.size()) * grainmod) / (1000 * grainmod));

        LOG(INFO) << "flushing ranges (num. blocks=" << rank(ir.size()) << ") with grainsize=" << grainsize;
        ranges::detail::IORange<typename InputRange::iterator,sdsl::rank_support_v5<1,1>,sdsl::select_support_mcl<1,1>> Range(
            ir.begin(), 0, ir.size(), std::get<0>(B).begin(),
            grainsize, rank, select,
            c_.position());

        flushOutputRange(c_, &Range, std::forward<Args>(args)...);
    }

    template <typename Cursor, typename InputRange, class ...Args>
    void flush(Cursor &c_, InputRange &ir1, terminals::TerminalRange &ir2, size_t grainmod, Args&& ...args) {
        LOG(INFO) << "calculating adjacent overlaps B2";
        sdsl::int_vector<2> B2 = adjacentDifference(ir2);
        LOG(INFO) << "calculating adjacent overlaps B1";
        std::pair<sdsl::int_vector<2>,sdsl::bit_vector> B1 = indexedAdjacentDifference<BW_0_K>(ir1);

        // initialise rank-select support
        LOG(INFO) << "initialising rank-select support";
        sdsl::rank_support_v5<1,1>    rank;
        sdsl::select_support_mcl<1,1> select;
        oneapi::tbb::parallel_invoke(
            [&]{ sdsl::util::init_support(rank,   &std::get<1>(B1)); },
            [&]{ sdsl::util::init_support(select, &std::get<1>(B1)); }
        );

        size_t grainsize = std::max(
            (size_t)1,
            (rank(ir1.size()) * grainmod) / (1000 * grainmod));

        LOG(INFO) << "flushing ranges (num. blocks=" << rank(ir1.size()) << ") with grainsize=" << grainsize;
        ranges::detail::InterleaveRanges<typename InputRange::iterator,sdsl::rank_support_v5<1,1>,sdsl::select_support_mcl<1,1>> Range(
            ir1.begin(), 0, ir1.size(), std::get<0>(B1).begin(),
            grainsize, rank, select,
            ir2, B2.begin(),
            ir1.getK(), ir1.getEffK(),
            c_.position());

        flushOutputRange(c_, &Range, std::forward<Args>(args)...);
    }

} // namespace flush
