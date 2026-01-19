#pragma once

#include <atomic>
#include <concepts>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <foonathan/memory/allocator_traits.hpp>
#include <foonathan/memory/memory_pool_collection.hpp>
#include <foonathan/memory/segregator.hpp>
#include <foonathan/memory/std_allocator.hpp>
#include <foonathan/memory/allocator_storage.hpp>
#include <foonathan/memory/namespace_alias.hpp>
#include <mio/mio.hpp>
#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/concurrent_hash_map.h>
#include <oneapi/tbb/parallel_for.h>
#include <maki/bitvectors.hpp>

namespace pmap
{
    
    namespace alloc
    {

        class TbbBlockAllocator
        {
        public:
            TbbBlockAllocator(size_t blockSize);
            memory::memory_block allocate_block();
            void deallocate_block(memory::memory_block b);
            size_t next_block_size() const noexcept;

        private:
            size_t _bsize;
        };

        using TbbMemoryPool = memory::memory_pool_collection<
            memory::node_pool, memory::log2_buckets, TbbBlockAllocator>;

        class MMapBlockAllocator
        {
            struct File
            {
                std::string file;
                mio::mmap_sink mmap;
                size_t bytes;
                static std::string &&initialiseFile(std::string &&file_, size_t bytes_);
                File(std::string &&file_, size_t bytes_);
                void sync();
                ~File();
                File(File &&) = default;
                File(File const &) = delete;
                inline File &operator=(File &&) = default;
                inline File &operator=(File const &) = delete;
                inline char *data() { return reinterpret_cast<char *>(mmap.data()); }
            };

            using FilePointer = std::unique_ptr<File>;
            using Heap = oneapi::tbb::concurrent_hash_map<void *, FilePointer>;

        public:
            MMapBlockAllocator() =delete;
            MMapBlockAllocator(size_t blockSize, std::string_view base);
            MMapBlockAllocator(MMapBlockAllocator const &) = delete;
            MMapBlockAllocator(MMapBlockAllocator &&) = default;
            MMapBlockAllocator &operator=(MMapBlockAllocator const &) = delete;
            MMapBlockAllocator &operator=(MMapBlockAllocator &&) = default;
            memory::memory_block allocate_block();
            void deallocate_block(memory::memory_block b);
            size_t next_block_size() const noexcept;

        private:
            size_t _bsize;
            std::string _base;
            Heap _heap;
            copyable_atomic_uint64_t _id;
        };

        using MMapMemoryPool = memory::memory_pool_collection<
            memory::node_pool, memory::log2_buckets, MMapBlockAllocator>;

        template <typename RawAllocator>
        class CappedAllocator
        {
        public:
            using allocator_type = typename memory::allocator_traits<RawAllocator>::allocator_type;

            CappedAllocator(RawAllocator &&alloc, size_t limit) noexcept
                : _alloc(std::move(alloc))
                , _limit(limit)
            {
            }

            inline allocator_type& get_allocator() noexcept
            { return _alloc; }

            inline allocator_type const& get_allocator() const noexcept
            { return _alloc; }

            inline bool use_allocate_node(size_t size, [[maybe_unused]] size_t alignment) noexcept
            { return size <= _limit; }

            inline bool use_allocate_array(size_t count, size_t size, [[maybe_unused]] size_t alignment) noexcept
            { return (count * size) <= _limit; }

        private:
            allocator_type _alloc;
            size_t _limit;
        };
    
        using DynamicAllocator = memory::binary_segregator<
            CappedAllocator<TbbMemoryPool>,
            memory::thread_safe_allocator<MMapMemoryPool, std::mutex>>;

        DynamicAllocator makeDynamicAllocator(size_t maxRAMAlloc, size_t maxDiskAllocL2, std::string_view baseDir);

    } // namespace alloc

    template <typename T>
    struct AtomicFixedSizeVector
    {
        AtomicFixedSizeVector(size_t capacity_)
            : _data(capacity_)
            , _size(0u)
            , _capacity(capacity_)
        {
        }

        AtomicFixedSizeVector(size_t capacity_, size_t start_)
            : _data(capacity_)
            , _size(start_)
            , _capacity(capacity_)
        {
        }

        AtomicFixedSizeVector(AtomicFixedSizeVector const &) = delete;
        AtomicFixedSizeVector(AtomicFixedSizeVector &&) = default;
        inline AtomicFixedSizeVector &operator=(AtomicFixedSizeVector const &) = delete;
        inline AtomicFixedSizeVector &operator=(AtomicFixedSizeVector &&) = default;

        inline T &allocNext()
        { return _data[_size++]; }

        inline uint64_t allocNextPosition()
        { return _size++; }

        inline void shrink()
        { _data.resize(static_cast<uint64_t>(_size)); }

        inline T &operator[](size_t i)
        { return _data[i]; }

        inline T const &operator[](size_t i) const
        { return _data[i]; }

        inline std::vector<T> &memoryview() noexcept
        { return _data; }

        inline uint64_t size() const
        { return static_cast<uint64_t>(_size); }

    protected:
        std::vector<T> _data;
        copyable_atomic_uint64_t _size;
        size_t _capacity;
    };

    namespace detail
    {
        
        struct SharedGate
        {
            enum Status : bool { OPEN = false, CLOSED = true };
            std::atomic<Status> status{OPEN};
            std::shared_mutex mutex{};
            bool tryObtain();
            void release();
            void wait();
        };

        template <uint8_t N>
        struct uint_leastN {
            using type = typename std::conditional<
                N <= 8, std::uint_least8_t,
                typename std::conditional<
                    N <= 16, std::uint_least16_t,
                    typename std::conditional<
                        N <= 32, std::uint_least32_t,
                        std::uint_least64_t
                    >::type
                >::type
            >::type;
        };

        template <uint8_t N>
        using uint_leastN_t = typename uint_leastN<N>::type;

    } // namespace detail

    /** BlockIndex
     * 
     * `(x,y)` coordinate pair holding block position in grid.
     */
    using BlockIndex = copyable_pair<size_t>;

    /** CellIndex
     * 
     * `(w,z)` coordinate pair holding cell position within a block.
     * Template parameter `c_width` specifies bit width of the integer,
     * up to a max of `31`. The cell is lockable via the standard 
     * mutex interface, i.e. `try_lock`, `lock` and `unlock`.
     */
    template <uint8_t c_width>
    class CellIndex {
        static_assert(c_width >= 1 && c_width <= 31, "Bit width must be between 1 and 31");
    public:
        using IndexType = detail::uint_leastN_t<2 * (c_width + 1)>;
        static_assert(sizeof(std::atomic<IndexType>)==sizeof(IndexType));

    protected:
        static constexpr IndexType maxValue = (1ULL << c_width) - 1,
                                   xMask = maxValue,
                                   yMask = maxValue << c_width,
                                   lockBit = 1ULL << (2 * c_width);
    
        std::atomic<IndexType> _index;

    public:
        constexpr CellIndex() : _index(0) {}
        constexpr CellIndex(size_t x_, size_t y_) : _index(x_ | (y_ << c_width)) {}
        constexpr CellIndex(BlockIndex c_) : CellIndex(c_.first, c_.second) {}

        inline constexpr IndexType data() const
        { return _index.load(std::memory_order_relaxed) & ~lockBit; }

        CellIndex& operator=(const CellIndex &rhs);

        inline constexpr bool operator==(const CellIndex &rhs) const
        { return data() == rhs.data(); }

        inline constexpr IndexType x() const noexcept
        { return _index.load(std::memory_order_relaxed) & xMask; }

        inline constexpr IndexType y() const noexcept
        { return (_index.load(std::memory_order_relaxed) & yMask) >> c_width; }

        constexpr void x(size_t x_) noexcept;

        constexpr void y(size_t y_) noexcept;

        template <typename S = IndexType>
        inline constexpr copyable_pair<S> asPair() const
        { return copyable_pair<S>(static_cast<S>(x()), static_cast<S>(y())); }

        // lock API

        bool try_lock() noexcept;

        void lock() noexcept;

        void unlock() noexcept;
    };

} // namespace pmap

namespace std
{

    template <uint8_t c_width>
    struct hash<pmap::CellIndex<c_width>>
    {
        size_t operator()(const pmap::CellIndex<c_width> &x) const
        { return std::hash<typename pmap::CellIndex<c_width>::IndexType>{}(x.data()); }
    };

} // namespace std

namespace pmap
{

    /**
     * DynamicTable
     * 
     * Uses a `RawAllocator` (usually some form of `memory::RawAllocator`) as the
     * basis for a hash map.
     */
    template <typename Key,         // Key type
              typename Int,         // Value integer type
              class Hash,           // Hash for key type
              class RawAllocator>   // Allocator for table buffer
    class DynamicTable : public memory::allocator_reference<RawAllocator>
    {
    public:
        struct Item {
            Key key;
            Int value;

            // initialise a record that has been claimed to value 1
            void initialise(const Key &key_);

            // initialise a record that has been claimed to a specific value
            void initialise(const Key &key_, Int value_);

            // increment count
            void update();

            // set value
            void update(Int value_);
        };

        using AllocatorRef = memory::allocator_reference<RawAllocator>;

    protected:
        std::atomic_uint32_t _load;
        uint8_t              _allocatedBytesL2;
        uint32_t             _allocated,
                             _capacity;
        Item                 *_data;
        Hash                 _hash;
        detail::SharedGate   _gate;

    public:
        static constexpr size_t growthFactor = 3,
                                initialAllocationL2 = 7,
                                initialAllocationBytes = (size_t)1u << initialAllocationL2;

        // constructor API ---------------------------------------------------------------

        DynamicTable(AllocatorRef alloc_);
        ~DynamicTable();

        DynamicTable(const DynamicTable&) =delete;
        DynamicTable& operator=(const DynamicTable&) =delete;

        // log2 of bytes allocation required for a capacity that holds `maxItems`
        // (assuming a log2 bucket scheme)
        static constexpr uint8_t pageSizeLog2(uint32_t maxItems);
        
        // number of items that can fit in a given byte allocation
        static constexpr uint32_t items(size_t bytes);

        // capacity for a given allocation size
        static constexpr uint32_t capacity(uint32_t allocated);

        // required allocation for a given capacity
        static constexpr uint32_t allocate(uint32_t capacity);

        // allocation --------------------------------------------------------------------
    protected:
        // increment allocation by growth factor
        void reAlloc();

        // increment allocation to a specific item capacity and copy data
        void reAlloc(uint32_t newItemCapacity);

        // allocate and initialise `numItems` items using the `RawAllocator`
        Item* allocateItems(uint32_t numItems);

    public:
        // current number of items in the hash map
        uint32_t load() const
        { return _load.load(); }

        // current capacity of the hash map
        uint32_t capacity() const noexcept
        { return _capacity; }

        // insertion, updating -----------------------------------------------------------
    protected:
        // atomic update item whose key matches `key`, attempting to insert
        // and default-initialise if it does not exist (fail if there is
        // no capacity)
        template <class ...Args>
        bool tryUpdate(const Key &key, Args&& ...args);

        // atomic increment item whose key matches `key`, attempting to insert
        // and default-initialise if it does not exist
        bool tryIncrement(const Key &key)
        { return tryUpdate(key); }

        // atomic increment load by 1, failing if it reaches capacity
        bool increaseLoad();
    
    public:
        // atomic increment value
        void increment(const Key &key);

        // atomic increment value
        template <typename K>
        void increment(const K &key)
        { increment(Key(key)); }

        // apply read-only function across hash map (re-allocation, and therefore
        // inserting, results in UB)
        template <class Fn>
        void forEach(Fn f_) const;
    };

    template <uint8_t blockBitWidth, typename Int>
    using DynamicCellTable = DynamicTable<
        CellIndex<blockBitWidth>, Int,
        std::hash<CellIndex<blockBitWidth>>,
        alloc::DynamicAllocator >;

    /** LazyMap
     * 
     * Container that defers allocation of data until it is requested.
     * The container is a template parameter. Pointers to allocated
     * containers are obtained under lock protection.
     */
    template <
        typename Key,
        typename Container,
        typename HashCompare = oneapi::tbb::tbb_hash_compare<Key> >
    struct LazyMap {
        using Value = std::unique_ptr<Container>;
        using Map   = oneapi::tbb::concurrent_hash_map<
            Key, Value,
            HashCompare,
            oneapi::tbb::tbb_allocator<std::pair<const Key, Value>> >;
    
    protected:
        Map grid;

        // element access -------------------------------------------------------
    public:
        /**
         * Access element at coordinate `key` under lock protection.
         * If it doesn't exist, allocate it using forwarded arguments.
         * Returns pointer to element, and `bool` indicating whether
         * an allocation took place.
         */
        template <class ...Args>
        std::pair<Container*, bool> get(const Key &key, Args&& ...args);

        /**
         * Apply function across const-qualified reference to containers.
         */
        template <class Fn>
        void forEach(Fn f_) const;
    };

    /** RankClustering
     * 
     * Approximate clustering of integers based on rank occurrence. Requires
     * upper bound on number of unique objects.
     */
    class RankClustering {
        bv::threaded::ConcurrentBitVector     coloursObserved;
        AtomicFixedSizeVector<size_t>         colours;
        std::vector<copyable_atomic_uint64_t> ranks;

    public:
        RankClustering(size_t numColours);
        
        // get rank of colour, or assign one if it hasn't been seen before
        size_t rank(size_t colour_);

        // get colour that occurred at a specific rank
        size_t colour(size_t rank_) const;
    };

    /** LazyGrid
     * 
     * defer allocation of grid blocks until they are requested, where blocks
     * are of width `2 ^ blockBitWidth`, and the container is a template
     * parameter
     */
    template <uint8_t blockBitWidth, typename Block>
    class LazyGridBase {
    public:
        using BlockType = Block;
        using Grid = LazyMap<BlockIndex,Block,copyable_pair_hash_64<size_t>>;

        struct Coordinate {
            BlockIndex block, cell;
            
            Coordinate(size_t x_, size_t y_);
        };

        static constexpr size_t blockSize = (size_t)1 << blockBitWidth;
    protected:
        static constexpr size_t cellMask = blockSize - 1u;

        Grid           grid;
        RankClustering axis;

    public:
        // initialise grid that holds `size_` different elements
        LazyGridBase(size_t size_);

        // get grid coordinates
        Coordinate coord(size_t r_, size_t c_);
    };

    template <uint8_t blockBitWidth, typename Block>
    class LazyGrid : public LazyGridBase<blockBitWidth,Block> {
    public:
        LazyGrid(size_t size_);

        // fetch container for a given grid coordinate
        template <class ...Args>
        std::pair<Block*,bool> get(const BlockIndex &coord, Args&& ...args);
    };

    template <uint8_t blockBitWidth, typename Int>
    class LazyGrid<blockBitWidth,DynamicCellTable<blockBitWidth,Int>>
        : public LazyGridBase<blockBitWidth,DynamicCellTable<blockBitWidth,Int>>
    {
        alloc::DynamicAllocator disk;

    public:
        LazyGrid(size_t size_, const std::string &baseDir_);

        std::pair<DynamicCellTable<blockBitWidth,Int>*,bool> get(const BlockIndex &coord);

        // apply a function to all const-qualified blocks
        template <typename Fn>
        void forEach(Fn f_) const;
    };

    template <uint8_t blockBitWidth, typename Int>
    using DiskSparseTable = LazyGrid<
        blockBitWidth,
        DynamicCellTable<blockBitWidth,Int> >;

    // Template Definitions =============================
    // ==================================================

    template<uint8_t c_width>
    CellIndex<c_width>& CellIndex<c_width>::operator=(const CellIndex &rhs) {
        IndexType currentValue = _index.load(std::memory_order_relaxed);
        _index.store((currentValue & lockBit) | rhs.data(),
                      std::memory_order_relaxed);
        return *this;
    }

    template<uint8_t c_width>
    constexpr void CellIndex<c_width>::x(size_t x_) noexcept {
        IndexType currentValue = _index.load(std::memory_order_relaxed);
        _index.store((currentValue & ~xMask) | (x_ & xMask),
                      std::memory_order_relaxed);
    }

    template<uint8_t c_width>
    constexpr void CellIndex<c_width>::y(size_t y_) noexcept {
        IndexType currentValue = _index.load(std::memory_order_relaxed);
        _index.store((currentValue & ~yMask) | ((y_ << c_width) & yMask),
                      std::memory_order_relaxed);
    }

    template <uint8_t c_width>
    bool CellIndex<c_width>::try_lock() noexcept {
        IndexType expected = _index.load(std::memory_order_relaxed);
        if (expected & lockBit)
            return false;
        return _index.compare_exchange_strong(expected, expected | lockBit, std::memory_order_acquire);
    }

    template <uint8_t c_width>
    void CellIndex<c_width>::lock() noexcept {
        IndexType expected;
        do {
            do {
                expected = _index.load(std::memory_order_relaxed);
            } while (expected & lockBit);
        } while (!_index.compare_exchange_weak(expected, expected | lockBit, std::memory_order_acquire));
    }

    template <uint8_t c_width>
    void CellIndex<c_width>::unlock() noexcept {
        _index.fetch_and(~lockBit, std::memory_order_release);
    }

    template<class Key, class Int, class Hash, class RawAllocator>
    void DynamicTable<Key,Int,Hash,RawAllocator>::Item::initialise(const Key &key_) {
        key = key_;
        value = 1;
    }

    template<class Key, class Int, class Hash, class RawAllocator>
    void DynamicTable<Key,Int,Hash,RawAllocator>::Item::initialise(const Key &key_, Int value_) {
        key = key_;
        value = value_;
    }

    template<class Key, class Int, class Hash, class RawAllocator>
    void DynamicTable<Key,Int,Hash,RawAllocator>::Item::update()
    { ++value; }

    template<class Key, class Int, class Hash, class RawAllocator>
    void DynamicTable<Key,Int,Hash,RawAllocator>::Item::update(Int value_)
    { value = value_; }

    template<class Key, class Int, class Hash, class RawAllocator>
    DynamicTable<Key,Int,Hash,RawAllocator>::DynamicTable(AllocatorRef alloc_)
        : AllocatorRef(alloc_)
        , _load(0)
        , _allocatedBytesL2(initialAllocationL2)
        , _allocated(items(initialAllocationBytes))
        , _capacity(capacity(_allocated))
        , _data(allocateItems(_allocated))
        , _hash()
        , _gate()
    {
    }

    template<class Key, class Int, class Hash, class RawAllocator>
    DynamicTable<Key,Int,Hash,RawAllocator>::~DynamicTable() {
        this->deallocate_node(_data, _allocated * sizeof(Item), alignof(Item));
    }

    template<class Key, class Int, class Hash, class RawAllocator>
    constexpr uint8_t DynamicTable<Key,Int,Hash,RawAllocator>::pageSizeLog2(uint32_t maxItems) {
        size_t bytes = (size_t)sizeof(Item) * (size_t)allocate(maxItems);
            uint8_t bl2 = ceil_log2(bytes);
            if (bl2 < initialAllocationL2)
                return initialAllocationL2;
            else
            {
                uint8_t rem3 = (bl2 - initialAllocationL2 + 2) / 3;
                return initialAllocationL2 + (3 * rem3);
            }
    }
    
    template<class Key, class Int, class Hash, class RawAllocator>
    constexpr uint32_t DynamicTable<Key,Int,Hash,RawAllocator>::items(size_t bytes)
    { return static_cast<uint32_t>(bytes / sizeof(Item)); }

    template<class Key, class Int, class Hash, class RawAllocator>
    constexpr uint32_t DynamicTable<Key,Int,Hash,RawAllocator>::capacity(uint32_t allocated)
    { return ((7 * allocated) + 9) / 10; }

    template<class Key, class Int, class Hash, class RawAllocator>
    constexpr uint32_t DynamicTable<Key,Int,Hash,RawAllocator>::allocate(uint32_t capacity)
    { return (10 * capacity) / 7; }

    template<class Key, class Int, class Hash, class RawAllocator>
    void DynamicTable<Key,Int,Hash,RawAllocator>::reAlloc() {
        _allocatedBytesL2 += growthFactor;
        try {
            reAlloc(items((size_t)1u << _allocatedBytesL2));
        } catch(...) {
            _allocatedBytesL2 -= growthFactor;
            throw;
        }
    }

    template<class Key, class Int, class Hash, class RawAllocator>
    void DynamicTable<Key,Int,Hash,RawAllocator>::reAlloc(uint32_t newItemCapacity) {
        Item *trg = allocateItems(newItemCapacity);

        std::swap(_data, trg);
        std::swap(_allocated, newItemCapacity);
        _capacity = capacity(_allocated);

        // insert old values
        _load = 0;
        oneapi::tbb::parallel_for(
            oneapi::tbb::blocked_range<Item *>(trg, trg + newItemCapacity, 100),
            [&](oneapi::tbb::blocked_range<Item *> const &r)
            {
                for (Item *p = r.begin(); p != r.end(); ++p)
                    if (p->value)
                        tryUpdate(p->key, p->value);
            });

        // free old allocation
        this->deallocate_node(trg, newItemCapacity * sizeof(Item), alignof(Item));
    }

    template<class Key, class Int, class Hash, class RawAllocator>
    DynamicTable<Key,Int,Hash,RawAllocator>::Item*
    DynamicTable<Key,Int,Hash,RawAllocator>::allocateItems(uint32_t numItems) {
        Item *trg = reinterpret_cast<Item *>(
            this->allocate_node(numItems * sizeof(Item), alignof(Item))
            );

        try {
            for (Item *p = trg; p != trg + numItems; ++p)
                ::new (p) Item();
        } catch (...) {
            this->deallocate_node(trg, numItems * sizeof(Item), alignof(Item));
            throw;
        }

        return trg;
    }

    template<class Key, class Int, class Hash, class RawAllocator>
    template<class ...Args>
    bool DynamicTable<Key,Int,Hash,RawAllocator>::tryUpdate(const Key &key, Args&& ...args) {
        uint32_t index = _hash(key) % _allocated,
                 step = 0;
        while (step < _allocated) {
            Item &x = _data[index];
            {
                std::lock_guard guard{ x.key };
                if (x.value == 0) {
                    if (increaseLoad()) {
                        x.initialise(key, std::forward<Args>(args)...);
                        return true;
                    } else return false;
                } else if (x.key == key) {
                    x.update(std::forward<Args>(args)...);
                    return true;
                }
            }
            index = (index + (++step)) % _allocated;
        }
        return false;
    }

    template<class Key, class Int, class Hash, class RawAllocator>
    bool DynamicTable<Key,Int,Hash,RawAllocator>::increaseLoad() {
        uint32_t currentLoad = _load.load();
        while (currentLoad < _capacity) {
            if (_load.compare_exchange_weak(
                    currentLoad, currentLoad + 1,
                    std::memory_order_acquire,
                    std::memory_order_relaxed))
            {
                return true;
            }
        }
        return false;
    }

    template<class Key, class Int, class Hash, class RawAllocator>
    void DynamicTable<Key,Int,Hash,RawAllocator>::increment(const Key &key) {
        for (;;) {
            // try to update counter
            {
                std::shared_lock olock(_gate.mutex);
                if (tryIncrement(key))
                    return;
            }

            // update failed, need to reallocate
            if (_gate.tryObtain()) {
                std::unique_lock xlock(_gate.mutex);
                if (_load == _capacity)
                    reAlloc();
                _gate.release();
            }
            else
                _gate.wait();
        }
    }

    template<class Key, class Int, class Hash, class RawAllocator>
    template <class Fn>
    void DynamicTable<Key,Int,Hash,RawAllocator>::forEach(Fn f_) const {
        for (Item const *p = _data; p != _data + _allocated; ++p)
            if (p->value)
                f_(p->key, p->value);
    }

    template<class Key, class Container, class HashCompare>
    template<class... Args>
    std::pair<Container*, bool> LazyMap<Key,Container,HashCompare>::get(const Key &key, Args &&...args) {
        typename Map::accessor p;
        if (!grid.insert(p, key))
            return std::make_pair(p->second.get(), false);
        else
        {
            p->second.reset(new Container(std::forward<Args>(args)...));
            return std::make_pair(p->second.get(), true);
        }
    }

    template<class Key, class Container, class HashCompare>
    template<class Fn>
    void LazyMap<Key,Container,HashCompare>::forEach(Fn f_) const {
        for (auto const &kvp : grid)
            f_(kvp.first, *kvp.second);
    }

    template<uint8_t blockBitWidth, class Block>
    LazyGridBase<blockBitWidth, Block>::Coordinate::Coordinate(size_t x_, size_t y_)
        : block(x_ >> blockBitWidth, y_ >> blockBitWidth)
        , cell(x_ & cellMask, y_ & cellMask)
    {}

    template<uint8_t blockBitWidth, class Block>
    LazyGridBase<blockBitWidth, Block>::LazyGridBase(size_t size_)
        : grid()
        , axis(size_)
    {}

    template<uint8_t blockBitWidth, class Block>
    LazyGridBase<blockBitWidth, Block>::Coordinate
    LazyGridBase<blockBitWidth, Block>::coord(size_t r_, size_t c_)
    { return Coordinate(axis.rank(r_), axis.rank(c_)); }

    template<uint8_t blockBitWidth, class Block>
    LazyGrid<blockBitWidth, Block>::LazyGrid(size_t size_)
        : LazyGridBase<blockBitWidth,Block>::LazyGridBase(size_)
    {}

    template<uint8_t blockBitWidth, class Block>
    template <class ...Args>
    std::pair<Block *, bool> LazyGrid<blockBitWidth, Block>::get(const BlockIndex &coord, Args&& ...args)
    { return this->grid.get(coord, std::forward<Args>(args)...); }

    template<uint8_t blockBitWidth, class Int> 
    LazyGrid<blockBitWidth, DynamicCellTable<blockBitWidth, Int>>::LazyGrid(size_t size_, const std::string &baseDir_)
        : LazyGridBase<blockBitWidth,DynamicCellTable<blockBitWidth,Int>>(size_)
        , disk(alloc::makeDynamicAllocator(
            2048,
            DynamicCellTable<blockBitWidth, Int>::pageSizeLog2(this->blockSize * this->blockSize),
            baseDir_))
    {}

    template<uint8_t blockBitWidth, class Int>
    std::pair<DynamicCellTable<blockBitWidth,Int>*, bool>
    LazyGrid<blockBitWidth, DynamicCellTable<blockBitWidth, Int>>::get(const BlockIndex &coord)
    { return this->grid.get(coord, disk); }

    template<uint8_t blockBitWidth, class Int>
    template<class Fn>
    void LazyGrid<blockBitWidth, DynamicCellTable<blockBitWidth, Int>>::forEach(Fn f_) const {
        this->grid.forEach([&](const BlockIndex &block_, const DynamicCellTable<blockBitWidth, Int> &data_)->void {
            data_.forEach([&](const CellIndex<blockBitWidth> &cell_, Int count)->void {
                f_(
                    this->axis.colour(block_.first  * this->blockSize + cell_.x()),
                    this->axis.colour(block_.second * this->blockSize + cell_.y()),
                    count
                );
            });
        });
    }

} // namespace pmap
