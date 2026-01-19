#pragma once
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <glog/logging.h>
#include <sdsl/vectors.hpp>
#include <gtl/phmap.hpp>
#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/blocked_range2d.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <indicators/block_progress_bar.hpp>
#include <oneapi/tbb/concurrent_hash_map.h>
#include <indicators/progress_bar.hpp>
#include <MurmurHash/MurmurHash3.h>
#include <maki/bitvectors.hpp>
#include <maki/fasta.hpp>
#include <maki/graph.hpp>
#include <maki/utils.hpp>
#include <maki/diskMap.hpp>

namespace dist
{

    using graph::ColourBuffer;
    using graph::ColourBufferConstIterator;
    using graph::ColourBufferIterator;
    using pmap::BlockIndex;
    using pmap::LazyGrid;

    using ColourType = uint64_t;
    using CountInt = uint32_t;
    using ColourVector = std::vector<ColourType, oneapi::tbb::tbb_allocator<ColourType>>;

    using ColourBufferRange = typename graph::TraversableGraph::ColourBufferRange;

    struct bitmask
    {
        bitmask(uint8_t numBits) : _w(numBits), _msk((1u << _w) - 1u) {}

        inline static constexpr size_t numHashes(uint8_t w_)
        { return 1u << w_; }

        template <typename T>
        inline constexpr T operator()(T const &x) const noexcept
        { return x & _msk; }

        template <typename T>
        inline constexpr T operator()(copyable_pair<T> const &x) const noexcept
        { return x.first & _msk; }

        inline uint8_t width() const noexcept
        { return _w; }

    protected:
        uint8_t _w;
        uint64_t _msk;
    };

    /** Lower Triangular Square Matrix
     *
     * Index structure
     * ---------------
     *
     *           0  1  2  3  4  5  6  7  8  9
     *         - - - - - - - - - - - - - - - -
     *     0  |  -
     *     1  |  0  -
     *     2  |  1  2  -
     *     3  |  3  4  5  -
     *     4  |  6  7  8  9  -
     *    ... |  ...
     *
     */
    template <typename Int>
    struct LTSquareMatrix
    {
        using IntType = copyable_atomic<Int>;

        LTSquareMatrix(size_t n_);
        
        static constexpr size_t sumUpToN(size_t n_);

        IntType& operator()(size_t r_, size_t c_);

        const IntType& operator()(size_t r_, size_t c_) const;

        size_t size() const { return _n; }

    protected:
        size_t _n;
        std::vector<IntType> _data;
    };

    using ThreadsZstd = oneapi::tbb::enumerable_thread_specific<ZstdDecompressor>;
    using ThreadsColourBuffers = oneapi::tbb::enumerable_thread_specific<graph::ColourBuffer>;

    template <class Functional>
    void ApplyAcrossColours(graph::ColourBufferRegistry const &colourBuffers, size_t grainsize, Functional &&Func);

    template <class Functional>
    void ApplyAcrossColours(graph::ColourBufferRegistry const &colourBuffers, size_t grainsize, Functional Func, const std::vector<std::pair<const size_t, graph::BufferIdPair>> &registerKeys);

    namespace pairwise
    {

        // counting shared kmers (dict-based) --------------------------------------------------------------------------- //

        template <uint8_t blockBitWidth, typename Block>
        struct SparseCountColourPairs
        {
            LazyGrid<blockBitWidth, Block> *obj;

            SparseCountColourPairs(LazyGrid<blockBitWidth, Block> *obj_) : obj(obj_) {}
            SparseCountColourPairs(SparseCountColourPairs &rhs, oneapi::tbb::split) : SparseCountColourPairs(rhs) {}

            void operator()(const ColourBufferRange &r) const;

        protected:
            void flushColours(ColourVector &colours, BlockIndex &block, Block *&mat) const;
        };

        template <uint8_t blockBitWidth, typename Int>
        pmap::DiskSparseTable<blockBitWidth,Int> countPairsSparse(const graph::TraversableGraph &g, std::string const &tempBase) {
            size_t numColoursBound = g.getMaxColour() + 1u;
            pmap::DiskSparseTable<blockBitWidth,Int> spmat(numColoursBound, tempBase);
            ApplyAcrossColours(g.colourBuffers, 1, SparseCountColourPairs(&spmat));
            return spmat;
        }

        // counting shared kmers (matrix-based) --------------------------------------------------------------------------- //

        using Matrix = LTSquareMatrix<uint32_t>;

        struct MatrixCountColourPairs
        {
            Matrix *m;
            bitmask proj;

        public:
            MatrixCountColourPairs(Matrix *m_, bitmask proj_);
            MatrixCountColourPairs(MatrixCountColourPairs &rhs, oneapi::tbb::split) : MatrixCountColourPairs(rhs) {}

            void operator()(const ColourBufferRange &r) const;

        protected:
            void flushColours(ColourVector &colours) const;
        };

        Matrix countPairsDense(graph::TraversableGraph const &g);

        // sample shared k-mers (matrix-based) --------------------------------------------------------------------------- //

        template <typename K, typename V>
        class PairSampler {
        protected:
            struct LessThan {
                bool operator()(const K &a_, const K &b_, size_t seed) const;
            private:
                uint32_t _hash(const K &x_, size_t seed) const;
            };
            size_t _mw;
        public:
            struct KeyValuePair { K key; V value; };
        protected:
            std::vector<KeyValuePair> data;
            bv::threaded::SpinRegionManager mtx;
        public:
            class Accessor {
                size_t _mw;
                std::vector<KeyValuePair> &data;
                bv::threaded::SpinRegionManager::Accessor lock;
                LessThan cmp;
            public:
                Accessor(PairSampler &obj_);
                void increment(size_t x_, size_t y_, const K &k_);
            };
            PairSampler(size_t n_);
            size_t width() const;
            const KeyValuePair& get(const K &x_, const K &y_) const;
        };

        struct ColourSampler : public PairSampler<uint64_t,uint32_t> {
            template <class Fn> void forEach(Fn f) const;
            template <class Op> void pforEach(size_t grainsize, oneapi::tbb::enumerable_thread_specific<Op> &tlocal) const;
            std::unordered_set<size_t> gatherKeys(uint8_t baseWidth);
        };

        struct SampleColourPairs {
            ColourSampler *obj;
            uint8_t keyWidth;
            dist::bitmask msk;

            SampleColourPairs(ColourSampler *obj_, uint8_t w_);
            SampleColourPairs(SampleColourPairs &rhs, oneapi::tbb::split);
            void operator()(const ColourBufferRange &r) const;
        };

        ColourSampler samplePairsDense(const graph::TraversableGraph &g);

        // results IO --------------------------------------------------------------------------- //

        template <uint8_t blockBitWidth, typename Int>
        void writeToCompressed(pmap::DiskSparseTable<blockBitWidth,Int> const &m, std::string filepath) {
            zstr::ofstream fh(filepath);
            m.forEach([&](ColourType c1, ColourType c2, CountInt v)
                      { fh << c1 << '\t' << c2 << '\t' << v << '\n'; });
            fh.close();
            LOG(INFO) << "pairwise counts written to " << filepath;
        }

        template <uint8_t blockBitWidth, typename Int>
        void writeToTsv(pmap::DiskSparseTable<blockBitWidth,Int> const &m, std::string filepath) {
            std::ofstream fh(filepath);
            m.forEach([&](ColourType c1, ColourType c2, CountInt v)
                      { fh << c1 << '\t' << c2 << '\t' << v << '\n'; });
            fh.close();
            LOG(INFO) << "pairwise counts written to " << filepath;
        }

        void writeToCompressed(const Matrix &m, std::string filepath);
        void writeToTsv(const Matrix &m, std::string filepath);
        template <typename InputMap> void writeToDisk(InputMap &&m, std::string filepath);
        template <typename InputMap> void writeToDisk(const InputMap &m, std::string filepath);

    } // namespace pairwise

    namespace query
    {
        // counting shared k-mers (dict) ------------------------------------------------------------------- //

        struct QuerySet {
            std::unordered_map<uint64_t,int64_t> keys;
            size_t size() const noexcept;
            void insert(uint64_t);
            int64_t index(uint64_t) const;
        };

        struct ConcurrentCounter {
            ConcurrentCounter() =default;
            oneapi::tbb::concurrent_unordered_map<
                copyable_pair<uint64_t>,
                std::atomic_uint32_t,
                copyable_pair_hash_64<uint64_t>> data;
            void incrementRange(uint64_t lhs,
                                sdsl::int_vector<0>::const_iterator rhsBegin,
                                sdsl::int_vector<0>::const_iterator rhsEnd);
        };

        struct ConcurrentBuffers {
            std::vector<ConcurrentCounter> data;
            ConcurrentBuffers() =default;
            ConcurrentBuffers(size_t n_);
            std::vector<ConcurrentCounter>::const_iterator begin() const;
            std::vector<ConcurrentCounter>::const_iterator end() const;
            size_t size() const noexcept;
            template <typename ...Args> auto& operator[](Args&& ...args) { return data.operator[](std::forward<Args>(args)...); }
            template <typename ...Args> const auto& operator[](Args&& ...args) const { return data.operator[](std::forward<Args>(args)...); }
            std::unordered_set<uint64_t> gatherKeys() const;
        };

        struct CountQueries {
            const QuerySet &queries;
            ConcurrentBuffers *buffers;
            bitmask msk;

            CountQueries(const QuerySet&, ConcurrentBuffers*, uint8_t);
            CountQueries(const CountQueries&) =default;
            CountQueries& operator=(const CountQueries&) =default;
            void operator()(const ColourBufferRange &r) const;
        protected:
            void searchQueries(const sdsl::int_vector<0>::const_iterator&,
                               const sdsl::int_vector<0>::const_iterator&) const;
        };

        ConcurrentBuffers searchGenomes(const QuerySet &queries, const graph::TraversableGraph &g);
        
    } // namespace query

    namespace unique
    {

        // counting unique kmers --------------------------------------------------------------------------- //

        using Colour = uint64_t;
        using Map = std::vector<copyable_atomic_uint32_t, oneapi::tbb::tbb_allocator<copyable_atomic_uint32_t>>;

        template <class Proj = std::identity>
        struct CountColours {
            Map *map;
            Proj proj;
            CountColours(Map *map_, Proj proj_) : map(map_), proj(proj_) {}
            CountColours(CountColours &rhs, oneapi::tbb::split) : CountColours(rhs) {}
            void operator()(ColourBufferRange &r) const;
            void flushColours(ColourVector &tmp) const;
        };

        template <class Proj, class Set>
        struct SubsampleCount {
            Proj                   proj;
            std::atomic_uint64_t   *map;
            std::vector<Set> const &trgs;
            SubsampleCount(Proj proj_, std::atomic_uint64_t *map_, std::vector<Set> const &trgs_) : proj(proj_), map(map_), trgs(trgs_) {}
            SubsampleCount(SubsampleCount &rhs, oneapi::tbb::split) : SubsampleCount(rhs) {}
            void operator()(ColourBufferRange &r) const;
            void checkEach(ColourVector &colours) const;
        };

        template <class Proj = std::identity>
        Map countUnique(graph::TraversableGraph const &g, Proj proj, size_t maxColour, size_t grainsize);

        template <typename Set>
        auto countKmerSubsamples(graph::TraversableGraph const &g, std::vector<Set> const &idSets, size_t grainsize);

        void writeToDisk(const Map &m, std::string filepath);
        void writeToCompressed(const Map &m, std::string filepath);
        void writeToTsv(const Map &m, std::string filepath);

    } // namespace unique

    namespace cooc
    {
        
        class Xlist {
            struct Hash {
                Xlist &arr;
                size_t operator()(const size_t &index) const;
            };

            struct Compare {
                Xlist &arr;
                bool operator()(const size_t &lhs, const size_t &rhs) const;
            };

            size_t _items, _itemWidth;
            sdsl::int_vector<0> _keys;
            Hash _hash;
            Compare _compare;
            gtl::flat_hash_map<size_t,uint32_t,Hash,Compare> _map;

        protected:
            // push a key
            template <typename InputIt>
            void push(InputIt begin_, InputIt end_);

            // increment the key that was just pushed, returning `true`
            // if the key was inserted (`false` if it already existed)
            bool update();

        public:
            Xlist(size_t bitWidth_, size_t itemWidth_);
            Xlist(size_t bitWidth_, size_t itemWidth_, size_t reserve_);

            template <typename InputIt>
            void increment(InputIt begin_, InputIt end_);

            template <typename Fn>
            void forEach(Fn f_) const;
        };

        class RankMap {
            using Key = std::vector<uint8_t, oneapi::tbb::tbb_allocator<uint8_t>>;

            struct ByteVectorHash {
                using is_transparent = void;
                size_t operator()(const Key& v) const noexcept;
                size_t operator()(std::string_view sv) const noexcept;
            };

            struct ByteVectorEq {
                using is_transparent = void;
                bool operator()(const Key& a, const Key& b) const noexcept;
                bool operator()(const Key& a, std::string_view b) const noexcept;
            };

            using Map = gtl::parallel_flat_hash_map<
                Key,                                 // Key
                uint32_t,                            // Mapped value
                ByteVectorHash,                      // Hash functor
                ByteVectorEq,                        // Equality comparator
                oneapi::tbb::tbb_allocator<std::pair<const Key, uint32_t>>,
                /*N=*/7,                             // number of submaps is 2^N => 128
                std::mutex                           // per-submap mutex for thread safety
            >;

            size_t _bypk;  // bytes per key
            Map _data;

        protected:
            Key compressU64(const ColourVector& src) const;
            ColourVector decompressU64(const Key& bytes) const;

        public:
            RankMap(size_t bitWidth_);
            void increment(const ColourVector &v);
            void set(const ColourVector &v, uint32_t z);
            Map::const_iterator find(const ColourVector &v) const;
            Map::const_iterator end() const;
            void flush(std::ostream&) const;
            template <typename Fn> void forEach(Fn f_) const;
            template <typename Fn> void pforEach(size_t grainsize, Fn f_) const;
            bool operator==(const RankMap &rhs) const;
        };

        struct BaseColourCounter {
            RankMap *m;
            BaseColourCounter(RankMap *m_);
            virtual void operator()(const ColourBufferRange &r) const =0;
        };

        struct CountColourLists : public BaseColourCounter {
            bitmask msk;
        public:
            CountColourLists(RankMap *m_, uint8_t w_);
            virtual void operator()(const ColourBufferRange &r) const override;
        };

        struct CountAllLists : public BaseColourCounter {
            using BaseColourCounter::BaseColourCounter;
            virtual void operator()(const ColourBufferRange &r) const override;
        };

        // count co-occurring annotated colours
        RankMap countPatterns(graph::TraversableGraph const &g, size_t grainsize);

        // count all co-occurring colours
        RankMap countRaw(graph::TraversableGraph const &g, size_t grainsize);

        struct RetrieveColouredKmers {
            const graph::TraversableGraph *g;
            const RankMap *qry;
            bitmask msk;
            std::ostream *os;
        public:
            RetrieveColouredKmers(const graph::TraversableGraph *g_, const RankMap *m_, uint8_t w_, std::ostream *os_);
            RetrieveColouredKmers(RetrieveColouredKmers &rhs, oneapi::tbb::split) : RetrieveColouredKmers(rhs) {}
            void operator()(const ColourBufferRange &r) const;
        protected:
            void printKmer(size_t edgeIndex, uint32_t value) const;
        };

        void printColours(graph::TraversableGraph const &g, const RankMap &qry, std::ostream &os);

    } // namespace cooc

    // Template Definitions =============================
    // ==================================================

    template<class Int>
    LTSquareMatrix<Int>::LTSquareMatrix(size_t n_)
        : _n(n_)
        , _data(sumUpToN(_n), IntType(Int(0u)))
    {
        LOG(INFO) << "allocating lower-triangular matrix with " << _n << " rows";
    }

    template<class Int>
    constexpr size_t LTSquareMatrix<Int>::sumUpToN(size_t n_) {
        if (n_ & 1u)
            return n_ * (n_ >> 1);
        else
            return (n_ >> 1) * (n_ - 1u);
    }

    template<class Int>
    LTSquareMatrix<Int>::IntType& LTSquareMatrix<Int>::operator()(size_t r_, size_t c_) {
        assert(r_ > c_);
        r_ = sumUpToN(r_);
        return _data[r_ + c_];
    }

    template<class Int>
    const LTSquareMatrix<Int>::IntType& LTSquareMatrix<Int>::operator()(size_t r_, size_t c_) const {
        assert(r_ > c_);
        r_ = sumUpToN(r_);
        return _data[r_ + c_];
    }

    template <class Functional>
    void ApplyAcrossColours(graph::ColourBufferRegistry const &colourBuffers, size_t grainsize, Functional &&Func) {
        std::vector registerKeys(colourBuffers.cbegin(), colourBuffers.cend());
        ApplyAcrossColours(colourBuffers, grainsize, std::move(Func), registerKeys);
    }

    template <class Functional>
    void ApplyAcrossColours(graph::ColourBufferRegistry const &colourBuffers, size_t grainsize, Functional Func, const std::vector<std::pair<const size_t, graph::BufferIdPair>> &registerKeys) {
        ThreadsZstd zstds;
        ThreadsColourBuffers buffers;

        // progress bar

        indicators::ProgressBar pbar{
            indicators::option::BarWidth{50},
            indicators::option::Start{"["},
            indicators::option::Fill{"="},
            indicators::option::Lead{">"},
            indicators::option::Remainder{" "},
            indicators::option::End{"]"},
            indicators::option::PrefixText{" Processing buffers "},
            indicators::option::ForegroundColor{indicators::Color::green},
            indicators::option::ShowElapsedTime{true},
            indicators::option::ShowRemainingTime{true},
            indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
            indicators::option::MaxProgress{registerKeys.size()}};

        // processing

        using c_iter = std::vector<std::pair<const size_t, graph::BufferIdPair>>::const_iterator;
        oneapi::tbb::parallel_for(
            oneapi::tbb::blocked_range<c_iter>(registerKeys.cbegin(), registerKeys.cend(), grainsize),
            [Func, &zstds, &buffers, &colourBuffers, &pbar](const oneapi::tbb::blocked_range<c_iter> &r)
            {
                graph::ColourBuffer &colours = buffers.local();
                ZstdDecompressor &zstd = zstds.local();

                for (auto const &[edgeIndex, kp] : r)
                {
                    // read buffer from disk

                    colours.resize(0);
                    colourBuffers.readColours(colours, kp, zstd);

                    // initialise rank/select support on edge boundaries

                    graph::TraversableGraph::ColourBufferRange Input(colours, 1u /* dummy */, edgeIndex);
                    Func(Input);

                    pbar.tick();
                }
            });

        pbar.mark_as_completed();
    }

    namespace pairwise
    {
        
        template<uint8_t blockBitWidth, typename Block>
        void SparseCountColourPairs<blockBitWidth, Block>::operator()(const ColourBufferRange &r) const {
            BlockIndex   matCoord;
            Block        *mat = nullptr;
            ColourVector tmp;

            auto it = r.begin(),
                 end = r.end();
            while (it != end) {
                do {
                    tmp.push_back(*it.colour);
                    ++it;
                } while ((it != end) && (*it.boundary == 0));
                flushColours(tmp, matCoord, mat);
            }
        }

        template<uint8_t blockBitWidth, typename Block>
        void SparseCountColourPairs<blockBitWidth, Block>::flushColours(ColourVector &colours, BlockIndex &block, Block *&mat) const {
            for (size_t i = 1; i < colours.size(); ++i) {
                for (size_t j = 0; j < i; ++j) {
                    typename LazyGrid<blockBitWidth,Block>::Coordinate coord
                        = obj->coord(colours[j], colours[i]);
                    
                    if ((!mat) || (coord.block != block)) {
                        std::pair<Block*, bool> alloc = obj->get(coord.block);
                        mat = alloc.first;
                        block = coord.block;
                    }

                    mat->increment(coord.cell);
                }
            }
            colours.clear();
        }

        template<class K, class V>
        bool PairSampler<K,V>::LessThan::operator()(const K &a_, const K &b_, size_t seed) const
        { return _hash(a_, seed) < _hash(b_, seed); }

        template<class K, class V>
        uint32_t PairSampler<K,V>::LessThan::_hash(const K &x_, size_t seed) const {
            uint8_t _data[sizeof(K)];
            std::memcpy(_data, &x_, sizeof(K));
            uint32_t _out;
            MurmurHash3_x86_32(_data, sizeof(K), seed, &_out);
            return _out;
        }

        template<class K, class V>
        PairSampler<K,V>::Accessor::Accessor(PairSampler<K,V> &obj_)
            : _mw(obj_._mw)
            , data(obj_.data)
            , lock(obj_.mtx)
            , cmp()
        {}

        template<class K, class V>
        void PairSampler<K,V>::Accessor::increment(size_t x_, size_t y_, const K &k_) {
            size_t index = x_ * _mw + y_;
            lock.access(index);
            KeyValuePair &item = data[index];
            if (item.value == 0) {
                item.key = k_;
                ++item.value;
            } else if (item.key == k_) {
                ++item.value;
            } else if (cmp(k_, item.key, /*seed*/index)) {
                item.key = k_;
                item.value = 1;
            }
        }

        template<class K, class V>
        PairSampler<K,V>::PairSampler(size_t n_)
            : _mw(n_)
            , data(n_ * n_)
            , mtx(n_ * n_, ceil_log2(n_))
        {}

        template<class K, class V>
        size_t PairSampler<K,V>::width() const { return _mw; }

        template<class K, class V>
        const PairSampler<K,V>::KeyValuePair& PairSampler<K,V>::get(const K &x_, const K &y_) const
        { return data[x_ * _mw + y_]; }

        template<class Fn>
        void ColourSampler::forEach(Fn f) const {
            for (size_t r = 0; r < _mw; ++r)
                for (size_t c = 0; c < _mw; ++c)
                    f(r+1, c+1, get(r, c));
        }

        template<class Op>
        void ColourSampler::pforEach(size_t grainsize, oneapi::tbb::enumerable_thread_specific<Op> &tlocal) const {
            oneapi::tbb::parallel_for(
                oneapi::tbb::blocked_range2d<size_t>(/*rows*/0, _mw, 1, /*cols*/0, _mw, grainsize),
                [&](const oneapi::tbb::blocked_range2d<size_t> &rng2d)->void {
                    Op &op = tlocal.local();
                    for (size_t r = rng2d.rows().begin(); r != rng2d.rows().end(); ++r)
                        for (size_t c = rng2d.cols().begin(); c != rng2d.cols().end(); ++c)
                            op(r+1, c+1, get(r, c));
                }
            );
        }

        template<class InputMap>
        void writeToDisk(InputMap &&m, std::string filepath) {
            if (filepath.ends_with(".gz"))
                writeToCompressed(std::move(m), filepath);
            else
                writeToTsv(std::move(m), filepath);
        }

        template <typename InputMap>
        void writeToDisk(const InputMap &m, std::string filepath)
        {
            if (filepath.ends_with(".gz"))
                writeToCompressed(m, filepath);
            else
                writeToTsv(m, filepath);
        }

    } // namespace pairwise

    namespace unique
    {

        template <class Proj>
        inline void CountColours<Proj>::operator()(ColourBufferRange &r) const
        {
            ColourVector tmp;

            auto it = r.begin(),
                 end = r.end();
            while (it != end) {
                do {
                    tmp.push_back(proj(*it.colour));
                    ++it;
                } while ((it != end) && (*it.boundary == 0));

                std::sort(tmp.begin(), tmp.end());
                auto end = std::unique(tmp.begin(), tmp.end());
                tmp.resize(end - tmp.begin());

                flushColours(tmp);
            }
        }

        template <>
        inline void CountColours<std::identity>::operator()(ColourBufferRange &r) const
        {
            auto
                it = r.begin(),
                end = r.end();
            while (it != end)
            {
                ++(*map)[*it.colour];
                ++it;
            }
        }

        template<class Proj>
        void CountColours<Proj>::flushColours(ColourVector &tmp) const {
            for (const Colour &key : tmp)
                ++(*map)[key];

            tmp.clear();
        }

        template <class Proj, class Set>
        void SubsampleCount<Proj,Set>::operator()(ColourBufferRange &r) const {
            ColourVector tmp;

            auto it  = r.begin(),
                 end = r.end();
            while (it != end) {
                do {
                    tmp.push_back(proj(*it.colour));
                    ++it;
                } while ((it != end) && (*it.boundary == 0));
                checkEach(tmp);
                tmp.clear();
            }
        }

        template <class Proj, class Set>
        void SubsampleCount<Proj,Set>::checkEach(ColourVector &colours) const {
            for (size_t i = 0; i < trgs.size(); ++i) {
                for (const Colour &c : colours) {
                    if (trgs[i].contains(c)) {
                        ++*(map + i);
                        break;
                    }
                }
            }
        }

        template <class Proj>
        Map countUnique(graph::TraversableGraph const &g, Proj proj, size_t maxColour, size_t grainsize)
        {
            Map counts(maxColour + 1u, copyable_atomic_uint32_t(0u));
            CountColours Count(&counts, proj);
            ApplyAcrossColours(g.colourBuffers, grainsize, Count);
            return counts;
        }

        template <typename Set>
        auto countKmerSubsamples(graph::TraversableGraph const &g, std::vector<Set> const &idSets, size_t grainsize) {
            auto p = std::make_unique<std::atomic_uint64_t[]>(idSets.size());
            SubsampleCount Count(bitmask(g.getBaseColourBitWidth()), p.get(), idSets);
            ApplyAcrossColours(g.colourBuffers, grainsize, Count);
            return p;
        }

    } // namespace unique

    namespace cooc
    {
        
        template <typename InputIt>
        void Xlist::push(InputIt begin_, InputIt end_) {
            // prepare space for another key
            size_t requiredSpace = (_items + 1) * _itemWidth;
            if (requiredSpace > _keys.size())
                _keys.insert(_keys.end(), requiredSpace - _keys.size(), 0);

            // write key values
            for (size_t i = 0;
                 (i < _itemWidth) && (begin_ != end_);
                 (void)++i, (void)++begin_)
            {
                _keys[(_items * _itemWidth) + i] = *begin_;
            }
        }

        template <typename InputIt>
        void Xlist::increment(InputIt begin_, InputIt end_) {
            push(begin_, end_);
            if (update()) ++_items;
        }

        template <typename Fn>
        void Xlist::forEach(Fn f_) const {
            auto begin_ = _keys.cbegin();
            for (const auto &[keyIndex, count] : _map)
                f_(begin_ + (keyIndex * _itemWidth),
                   begin_ + ((keyIndex + 1) * _itemWidth),
                   count);
        }

        template <typename Fn>
        void RankMap::forEach(Fn f_) const {
            for (const auto &[key, count] : _data) {
                auto bigKey = decompressU64(key);
                f_(bigKey.cbegin(), bigKey.cend(), count);
            }
        }

        template <typename Fn>
        void RankMap::pforEach(size_t grainsize, Fn f_) const {
            oneapi::tbb::parallel_for(
                oneapi::tbb::blocked_range(_data.cbegin(), _data.cend(), grainsize),
                [&](const auto &rng)->void {
                    for (const auto &[key, count] : rng) {
                        auto bigKey = decompressU64(key);
                        f_(bigKey.cbegin(), bigKey.cend(), count);
                    }
                }
            );
        }

    } // namespace cooc
    
} // namespace dist
