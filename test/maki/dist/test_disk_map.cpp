#include <atomic>
#include <chrono>
#include <filesystem>
#include <memory>
#include <thread>
#include <tuple>
#include <string>
#include <unordered_map>
#include <vector>
#include <gtest/gtest.h>
#include <maki/diskMap.hpp>
#include <maki/utils.hpp>

class MockTbbBlockAllocator : public pmap::alloc::TbbBlockAllocator {
public:
    MockTbbBlockAllocator(size_t blockSize) : pmap::alloc::TbbBlockAllocator(blockSize) {}
};

class MockMMapBlockAllocator {
    std::string tempbase;
    pmap::alloc::MMapBlockAllocator alloc;

public:
    MockMMapBlockAllocator(size_t blockSize)
        : tempbase(tempio::create_temporary_directory())
        , alloc(blockSize, tempbase)
    {}

    ~MockMMapBlockAllocator() {
        std::filesystem::remove_all(tempbase);
    }

    memory::memory_block allocate_block() { return alloc.allocate_block(); }

    void deallocate_block(memory::memory_block b) { return alloc.deallocate_block(b); }
};

using BlockAllocatorTypes = ::testing::Types<MockTbbBlockAllocator, MockMMapBlockAllocator>;

template <typename AllocatorType>
class MemoryAllocatorTest : public ::testing::Test {
protected:
    static constexpr size_t kBlockSize = 128;
    std::unique_ptr<AllocatorType> allocator;

    void SetUp() override {
        allocator = std::make_unique<AllocatorType>(kBlockSize);
    }
};

TYPED_TEST_SUITE_P(MemoryAllocatorTest);

// Test 1: Allocation returns non-null pointer
TYPED_TEST_P(MemoryAllocatorTest, AllocateBlockReturnsNonNull) {
    auto block = this->allocator->allocate_block();
    EXPECT_NE(block.memory, nullptr);
    this->allocator->deallocate_block(block);
}

// Test 2: Allocation returns correct block size
TYPED_TEST_P(MemoryAllocatorTest, AllocateBlockReturnsCorrectSize) {
    auto block = this->allocator->allocate_block();
    EXPECT_EQ(block.size, this->kBlockSize);
    this->allocator->deallocate_block(block);
}

// Test 3: Multiple allocations return distinct blocks
TYPED_TEST_P(MemoryAllocatorTest, MultipleAllocationsAreDistinct) {
    auto block1 = this->allocator->allocate_block();
    auto block2 = this->allocator->allocate_block();
    EXPECT_NE(block1.memory, block2.memory);
    this->allocator->deallocate_block(block1);
    this->allocator->deallocate_block(block2);
}

REGISTER_TYPED_TEST_SUITE_P(
    MemoryAllocatorTest,
    AllocateBlockReturnsNonNull,
    AllocateBlockReturnsCorrectSize,
    MultipleAllocationsAreDistinct
);

INSTANTIATE_TYPED_TEST_SUITE_P(MyAllocators, MemoryAllocatorTest, BlockAllocatorTypes);

// Template fixture: parameterized by AllocatorType and IntegerType
template <typename AllocatorType, typename IntegerType>
class MemoryAllocatorConcurrencyTest : public ::testing::Test {
protected:
    static constexpr size_t kArrayLength = 1024;
    static constexpr size_t kBlockSize = sizeof(IntegerType) * kArrayLength;
    static constexpr size_t kThreadCount = 8;

    std::unique_ptr<AllocatorType> allocator;

    void SetUp() override {
        allocator = std::make_unique<AllocatorType>(kBlockSize);
    }
};

// Define integer types to test
using IntegerTypes = ::testing::Types<int8_t, int16_t, int32_t, int64_t, uint32_t, uint64_t>;

// Macro to generate typed test suite for a given allocator
#define DEFINE_CONCURRENCY_TESTS_FOR_ALLOCATOR(AllocatorType, TestSuiteName) \
template <typename T> \
using TestSuiteName = MemoryAllocatorConcurrencyTest<AllocatorType, T>; \
TYPED_TEST_SUITE(TestSuiteName, IntegerTypes); \
TYPED_TEST(TestSuiteName, CanWriteAndReadIntegerArray) { \
    auto block = this->allocator->allocate_block(); \
    ASSERT_EQ(block.size, this->kBlockSize); \
    ASSERT_NE(block.memory, nullptr); \
    TypeParam* array = static_cast<TypeParam*>(block.memory); \
    for (size_t i = 0; i < this->kArrayLength; ++i) { \
        array[i] = static_cast<TypeParam>(i * 2); \
    } \
    for (size_t i = 0; i < this->kArrayLength; ++i) { \
        EXPECT_EQ(array[i], static_cast<TypeParam>(i * 2)); \
    } \
    this->allocator->deallocate_block(block); \
}\
TYPED_TEST(TestSuiteName, ConcurrentAtomicIncrements) { \
    auto block = this->allocator->allocate_block(); \
    ASSERT_NE(block.memory, nullptr); \
    ASSERT_EQ(block.size, TestSuiteName<TypeParam>::kBlockSize); \
    using AtomicType = std::atomic<TypeParam>; \
    AtomicType* array = reinterpret_cast<AtomicType*>(block.memory); \
    for (size_t i = 0; i < this->kArrayLength; ++i) { \
        ::new (array + i) AtomicType(); \
    } \
    for (size_t i = 0; i < TestSuiteName<TypeParam>::kArrayLength; ++i) { \
        array[i].store(0, std::memory_order_relaxed); \
    } \
    std::vector<std::thread> threads; \
    for (size_t t = 0; t < TestSuiteName<TypeParam>::kThreadCount; ++t) { \
        threads.emplace_back([&] { \
            for (size_t i = 0; i < TestSuiteName<TypeParam>::kArrayLength; ++i) { \
                array[i].fetch_add(1, std::memory_order_relaxed); \
            } \
        }); \
    } \
    for (auto& t : threads) t.join(); \
    for (size_t i = 0; i < TestSuiteName<TypeParam>::kArrayLength; ++i) { \
        EXPECT_EQ(array[i].load(std::memory_order_relaxed), TestSuiteName<TypeParam>::kThreadCount); \
    } \
    this->allocator->deallocate_block(block); \
}

DEFINE_CONCURRENCY_TESTS_FOR_ALLOCATOR(MockTbbBlockAllocator, TbbBlockConcurrencyTests)
DEFINE_CONCURRENCY_TESTS_FOR_ALLOCATOR(MockMMapBlockAllocator, MMapBlockConcurrencyTests)

class DynamicAllocatorTests : public ::testing::TestWithParam<std::tuple<size_t,size_t>> {
protected:
    std::string tempbase;
    std::unique_ptr<pmap::alloc::DynamicAllocator> allocator;

    using Item = std::atomic<uint64_t>;
    static constexpr size_t ItemSize = sizeof(Item), kThreadCount = 8;

    void SetUp() override {
        tempbase = tempio::create_temporary_directory();

        auto [ramLimitItems, maxLimitItems] = GetParam();
        allocator = std::make_unique<pmap::alloc::DynamicAllocator>(
            pmap::alloc::makeDynamicAllocator(
                ItemSize * ramLimitItems,
                ceil_log2(ItemSize * maxLimitItems),
                tempbase));
    }

    void TearDown() override {
        std::filesystem::remove_all(tempbase);
    }

    void InitialiseMemory(Item *data, size_t numItems) {
        for (Item *p = data; p != data + numItems; ++p)
            ::new (p) Item();
    }
};

TEST_P(DynamicAllocatorTests, RAMAllocationAlignment) {
    size_t items = std::get<0>(GetParam()) - 1;
    void *p = this->allocator->allocate_node(items * sizeof(Item), alignof(Item));
    EXPECT_EQ(reinterpret_cast<std::uintptr_t>(p) % alignof(Item), 0);
    this->allocator->deallocate_node(p, items * sizeof(Item), alignof(Item));
}

TEST_P(DynamicAllocatorTests, DiskAllocationAlignment) {
    size_t items = std::get<0>(GetParam()) + 1;
    void *p = this->allocator->allocate_node(items * sizeof(Item), alignof(Item));
    EXPECT_EQ(reinterpret_cast<std::uintptr_t>(p) % alignof(Item), 0);
    this->allocator->deallocate_node(p, items * sizeof(Item), alignof(Item));
}

TEST_P(DynamicAllocatorTests, RAMAllocationUse) {
    size_t items = std::get<0>(GetParam()) - 1;
    Item *data = reinterpret_cast<Item*>(
        this->allocator->allocate_node(items * sizeof(Item), alignof(Item)));
    InitialiseMemory(data, items);
    
    for (Item *p = data; p < data + items; ++p) {
        EXPECT_EQ(p->load(), (size_t)0);
        p->store(1u);
    }
    for (Item *p = data; p < data + items; ++p)
        EXPECT_EQ(p->load(), (size_t)1);
    
    this->allocator->deallocate_node(data, items * sizeof(Item), alignof(Item));
}

TEST_P(DynamicAllocatorTests, DiskAllocationUse) {
    size_t items = std::get<0>(GetParam()) + 1;
    Item *data = reinterpret_cast<Item*>(
        this->allocator->allocate_node(items * sizeof(Item), alignof(Item)));
    InitialiseMemory(data, items);
    
    for (Item *p = data; p < data + items; ++p) {
        EXPECT_EQ(p->load(), (size_t)0);
        p->store(1u);
    }
    for (Item *p = data; p < data + items; ++p)
        EXPECT_EQ(p->load(), (size_t)1);
    
    this->allocator->deallocate_node(data, items * sizeof(Item), alignof(Item));
}

TEST_P(DynamicAllocatorTests, MaxAllocationConcurrentUse) {
    size_t items = std::get<1>(GetParam());
    Item *data = reinterpret_cast<Item*>(
        this->allocator->allocate_node(items * sizeof(Item), alignof(Item)));
    InitialiseMemory(data, items);
    
    std::vector<std::thread> threads;
    for (size_t t = 0; t < kThreadCount; ++t) {
        threads.emplace_back([&] { \
            for (size_t i = 0; i < items; ++i) {
                data[i].fetch_add(1, std::memory_order_relaxed);
            }
        });
    }

    for (auto& t : threads) t.join();
    for (size_t i = 0; i < items; ++i) {
        EXPECT_EQ(data[i].load(std::memory_order_relaxed), kThreadCount);
    }
    
    this->allocator->deallocate_node(data, items * sizeof(Item), alignof(Item));
}

INSTANTIATE_TEST_SUITE_P(
    DynamicAllocatorParameters,
    DynamicAllocatorTests,
    ::testing::Combine(
        ::testing::Values(2, 10),
        ::testing::Values(100, 500)
    )
);

class AtomicFixedSizeVectorTests : public testing::Test
{
protected:
    AtomicFixedSizeVectorTests()
    {
    }

    void RunCompetitiveAllocTest(size_t size_, size_t grainsize_) {
        std::vector<size_t> ids(size_);
        std::iota(ids.begin(), ids.end(), 0u);

        pmap::AtomicFixedSizeVector<size_t> arr(size_);
        oneapi::tbb::parallel_for(
            oneapi::tbb::blocked_range(ids.begin(), ids.end(), grainsize_),
            [&](oneapi::tbb::blocked_range<typename std::vector<size_t>::iterator> &r)
            {
                for (size_t val_ : r)
                {
                    size_t &dest = arr.allocNext();
                    dest = val_;
                }
            });

        std::vector<size_t> &out = arr.memoryview();
        std::sort(out.begin(), out.end());

        ASSERT_EQ(out, ids);
    }
};

TEST_F(AtomicFixedSizeVectorTests, Size10Grainsize1) { RunCompetitiveAllocTest(10, 1); }
TEST_F(AtomicFixedSizeVectorTests, Size10Grainsize10) { RunCompetitiveAllocTest(10, 10); }
TEST_F(AtomicFixedSizeVectorTests, Size1Grainsize1) { RunCompetitiveAllocTest(1, 1); }
TEST_F(AtomicFixedSizeVectorTests, Size1Grainsize10) { RunCompetitiveAllocTest(1, 10); }
TEST_F(AtomicFixedSizeVectorTests, Size10000Grainsize1) { RunCompetitiveAllocTest(10000, 1); }
TEST_F(AtomicFixedSizeVectorTests, Size10000Grainsize100) { RunCompetitiveAllocTest(10000, 100); }

class SharedGateTest : public ::testing::Test {
protected:
    pmap::detail::SharedGate gate;
    std::atomic<int> sharedCounter{0};
    std::atomic<int> stateChangerCount{0};
};

// Test that multiple threads can acquire shared_lock concurrently
TEST_F(SharedGateTest, MultipleSharedLocks) {
    const int threadCount = 10;
    std::vector<std::thread> threads;
    std::atomic<int> concurrentReaders{0};
    std::atomic<int> maxConcurrentReaders{0};

    for (int i = 0; i < threadCount; ++i) {
        threads.emplace_back([&] {
            std::shared_lock lock(gate.mutex);
            int current = ++concurrentReaders;
            maxConcurrentReaders.store(std::max(maxConcurrentReaders.load(), current));
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
            --concurrentReaders;
        });
    }

    for (auto& t : threads) t.join();

    EXPECT_EQ(maxConcurrentReaders.load(), threadCount);
}

// Test that only one thread performs the state change
TEST_F(SharedGateTest, SingleStateChanger) {
    const int threadCount = 20;
    std::vector<std::thread> threads;

    for (int i = 0; i < threadCount; ++i) {
        threads.emplace_back([&] {
            if (gate.tryObtain()) {
                std::unique_lock lock(gate.mutex);
                ++stateChangerCount;
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                gate.release();
            } else {
                gate.wait();
            }
        });
    }

    for (auto& t : threads) t.join();

    EXPECT_EQ(stateChangerCount.load(), 1);
}

// Test high concurrency with shared access and a single state change
TEST_F(SharedGateTest, HighConcurrencyWithStateChange) {
    const int readerThreads = 50;
    std::vector<std::thread> threads;

    // Start readers
    for (int i = 0; i < readerThreads; ++i) {
        threads.emplace_back([&] {
            std::shared_lock lock(gate.mutex);
            ++sharedCounter;
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            --sharedCounter;
        });
    }

    // Wait a bit to ensure readers are active
    std::this_thread::sleep_for(std::chrono::milliseconds(20));

    // Start state changer
    std::thread changer([&] {
        if (gate.tryObtain()) {
            std::unique_lock lock(gate.mutex);
            ++stateChangerCount;
            EXPECT_EQ(sharedCounter.load(), 0); // Ensure no readers
            gate.release();
        }
    });

    for (auto& t : threads) t.join();
    changer.join();

    EXPECT_EQ(stateChangerCount.load(), 1);
}

// Test multiple state change attempts under load
TEST_F(SharedGateTest, MultipleStateChangeAttempts) {
    const int changer_threads = 10;
    std::vector<std::thread> threads;

    for (int i = 0; i < changer_threads; ++i) {
        threads.emplace_back([&] {
            if (gate.tryObtain()) {
                std::unique_lock lock(gate.mutex);
                ++stateChangerCount;
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
                gate.release();
            } else {
                gate.wait();
            }
        });
    }

    for (auto& t : threads) t.join();

    EXPECT_EQ(stateChangerCount.load(), 1);
}

TEST_F(SharedGateTest, RapidStateChangesUnderHighLoad) {
    const int threadCount = 100;
    const int iterations = 50;
    std::atomic<int> successfulChanges{0};
    std::atomic<int> concurrentChanges{0};

    auto worker = [&] {
        for (int i = 0; i < iterations; ++i) {
            {
                std::shared_lock shlock(gate.mutex);
                std::this_thread::sleep_for(std::chrono::milliseconds(1)); // simulate read work
            }

            if (gate.tryObtain()) {
                std::unique_lock unique_lock(gate.mutex);
                ++concurrentChanges;
                ASSERT_EQ(concurrentChanges.load(), 1) << "More than one thread entered critical section!";
                ++successfulChanges;
                std::this_thread::sleep_for(std::chrono::milliseconds(2)); // simulate write work
                --concurrentChanges;
                gate.release();
            } else {
                gate.wait();
            }
        }
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < threadCount; ++i) threads.emplace_back(worker);
    for (auto& t : threads) t.join();
}

#define DEFINE_CELL_INDEX_TESTS_FOR_BIT_WIDTH(BitWidth) \
class CellIndexTests##BitWidth : public ::testing::Test { \
protected: \
    pmap::CellIndex<BitWidth> cell; \
}; \
TEST_F(CellIndexTests##BitWidth, SingleThreadedReadWrite) { \
    cell.lock(); \
    cell.x(BitWidth); \
    cell.y(BitWidth - 1); \
    size_t x = cell.x(); \
    size_t y = cell.y(); \
    cell.unlock(); \
    EXPECT_EQ(x, BitWidth); \
    EXPECT_EQ(y, BitWidth - 1); \
} \
TEST_F(CellIndexTests##BitWidth, ConcurrentAccessWithLocking) { \
    const int threadCount = 20; \
    std::vector<std::thread> threads; \
    std::atomic<int> successfulReads{0}; \
    std::atomic<int> successfulWrites{0}; \
    for (int i = 0; i < threadCount; ++i) { \
        threads.emplace_back([&, i] { \
            if (i % 2 == 0) { \
                cell.lock(); \
                cell.x(i); \
                cell.y(i + 1); \
                cell.unlock(); \
                ++successfulWrites; \
            } else { \
                cell.lock(); \
                volatile size_t x = cell.x(); \
                volatile size_t y = cell.y(); \
                (void)x; (void)y; \
                cell.unlock(); \
                ++successfulReads; \
            } \
        }); \
    } \
    for (auto& t : threads) t.join(); \
    EXPECT_EQ(successfulReads + successfulWrites, threadCount); \
} \
TEST_F(CellIndexTests##BitWidth, TryLockBehavior) { \
    std::atomic<int> lockSuccessCount{0}; \
    const int threadCount = 50; \
    std::vector<std::thread> threads; \
    for (int i = 0; i < threadCount; ++i) { \
        threads.emplace_back([&] { \
            if (cell.try_lock()) { \
                ++lockSuccessCount; \
                std::this_thread::sleep_for(std::chrono::milliseconds(2)); \
                cell.unlock(); \
            } \
        }); \
    } \
    for (auto& t : threads) t.join(); \
    EXPECT_GE(lockSuccessCount.load(), 1); \
    EXPECT_LE(lockSuccessCount.load(), threadCount); \
} \
TEST_F(CellIndexTests##BitWidth, ContendedAccessTest) { \
    std::atomic<int> state{0}; \
    auto firstThread = std::thread([&] { \
        cell.lock(); \
        state.store(1); \
        state.notify_one(); \
        std::this_thread::sleep_for(std::chrono::milliseconds(100)); \
        cell.x((1u << BitWidth) - 1); \
        cell.y((1u << BitWidth) - 1); \
        cell.unlock(); \
    }); \
    auto secondThread = std::thread([&] { \
        state.wait(0); \
        cell.lock(); \
        EXPECT_EQ(cell.x(), (1u << BitWidth) - 1); \
        EXPECT_EQ(cell.y(), (1u << BitWidth) - 1); \
        cell.unlock(); \
    }); \
    secondThread.join(); \
    firstThread.join(); \
} \
TEST_F(CellIndexTests##BitWidth, HighConcurrencyStressTest) { \
    const int threadCount = 100; \
    std::vector<std::thread> threads; \
    for (int i = 0; i < threadCount; ++i) { \
        threads.emplace_back([&] { \
            for (size_t j = 0; j < BitWidth; ++j) { \
                cell.lock(); \
                cell.x((1u << j) - 1); \
                cell.y((1u << j) - 1); \
                EXPECT_EQ(cell.x(), (1u << j) - 1); \
                EXPECT_EQ(cell.y(), (1u << j) - 1); \
                cell.unlock(); \
            } \
        }); \
    } \
    for (auto& t : threads) t.join(); \
}

DEFINE_CELL_INDEX_TESTS_FOR_BIT_WIDTH(1);
DEFINE_CELL_INDEX_TESTS_FOR_BIT_WIDTH(4);
DEFINE_CELL_INDEX_TESTS_FOR_BIT_WIDTH(8);
DEFINE_CELL_INDEX_TESTS_FOR_BIT_WIDTH(16);
DEFINE_CELL_INDEX_TESTS_FOR_BIT_WIDTH(24);
DEFINE_CELL_INDEX_TESTS_FOR_BIT_WIDTH(25);
DEFINE_CELL_INDEX_TESTS_FOR_BIT_WIDTH(31);

template <typename Key>
struct LockableKey {
    std::mutex mtx;
    Key data;

    LockableKey() =default;
    LockableKey(const Key &data_) : mtx(), data(data_) {}
    LockableKey(const LockableKey &rhs) : mtx(), data(rhs.data) {}
    template <typename K> LockableKey(const K &k_) : LockableKey(Key(k_)) {}

    LockableKey& operator=(const LockableKey &rhs) { data = rhs.data; return *this; }
    LockableKey& operator=(const Key &rhs) { data = rhs; return *this; }
    bool operator!=(const LockableKey &rhs) const { return data != rhs.data; }
    bool operator!=(const Key &rhs) const { return data != rhs; }
    bool operator==(const LockableKey &rhs) const { return data == rhs.data; }
    bool operator==(const Key &rhs) const { return data == rhs; }

    bool try_lock() noexcept { return mtx.try_lock(); }
    void lock() noexcept { return mtx.lock(); }
    void unlock() noexcept { return mtx.unlock(); }
};

namespace std
{
    template <typename Key>
    struct hash<LockableKey<Key>> {
        size_t operator()(const LockableKey<Key> &x) const
        { return std::hash<Key>{}(x.data); }
    };
} // namespace std

struct CustomKey {
    int id;
    std::string name;

    CustomKey() =default;
    CustomKey(int id_) : id(id_), name(std::to_string(id_)) {}
    CustomKey(int id_, const std::string &name_) : id(id_), name(name_) {}
    CustomKey(const CustomKey &rhs) =default;
    CustomKey(CustomKey &&rhs) =default;

    CustomKey& operator=(const CustomKey &rhs) =default;
    CustomKey& operator=(CustomKey&& rhs) =default;

    bool operator==(const CustomKey& other) const {
        return id == other.id && name == other.name;
    }
};

namespace std {
    template <>
    struct hash<CustomKey> {
        std::size_t operator()(const CustomKey& k) const {
            return std::hash<int>()(k.id) ^ (std::hash<std::string>()(k.name) << 1);
        }
    };
}

// Test fixture template
template <typename T>
class DynamicTableTest : public ::testing::Test {
protected:
    using Key = typename T::first_type;
    using Int = typename T::second_type;
    using Alloc = pmap::alloc::TbbMemoryPool;
    using Table = pmap::DynamicTable<Key, Int, std::hash<Key>, Alloc>;
    
    static constexpr uint32_t initialSpace = Table::items(Table::initialAllocationBytes),
                              initialCapacity = Table::capacity(initialSpace);

    size_t maxAllocL2;
    Alloc alloc;
    Table table;

    DynamicTableTest()
        : maxAllocL2(Table::pageSizeLog2(2 * std::max((uint32_t)6400, (uint32_t)32 * initialCapacity)))
        , alloc(1u << maxAllocL2, (1u << maxAllocL2) * maxAllocL2)
        , table(alloc)
    {}

    std::unordered_map<Key,Int> ExportData() const {
        std::unordered_map<Key,Int> data;
        table.forEach([&](const Key &key, const Int &val)->void { data[key] = val; });
        return data;
    }
};

// Define test types
using TestTypes = ::testing::Types<
    std::pair<LockableKey<int>, int>,
    std::pair<LockableKey<CustomKey>, uint64_t>
>;

TYPED_TEST_SUITE(DynamicTableTest, TestTypes);

TYPED_TEST(DynamicTableTest, InitialSpace) {
    std::cout << this->initialSpace << " " << this->initialCapacity << std::endl;
}

TYPED_TEST(DynamicTableTest, InsertAndGetAKey) {
    using Key = typename TypeParam::first_type;
    using Int = typename TypeParam::second_type;

    Key key = Key{ 1 };
    this->table.increment(key);

    auto data = this->ExportData();
    EXPECT_EQ(data[key], static_cast<Int>(1));
    EXPECT_EQ(this->table.load(), 1);
}

// Basic single-threaded increment test
TYPED_TEST(DynamicTableTest, SingleThreadedIncrement) {
    using Key = typename TypeParam::first_type;
    using Int = typename TypeParam::second_type;

    Key key1 = Key{ 1 };
    Key key2 = Key{ 2 };

    this->table.increment(key1);
    this->table.increment(key1);
    this->table.increment(key2);

    auto data = this->ExportData();
    EXPECT_EQ(data[key1], static_cast<Int>(2));
    EXPECT_EQ(data[key2], static_cast<Int>(1));
    EXPECT_EQ(this->table.load(), 2);
}

// Concurrent increments on same key
TYPED_TEST(DynamicTableTest, ConcurrentIncrementsSameKey) {
    using Key = typename TypeParam::first_type;
    using Int = typename TypeParam::second_type;

    Key key = Key{ 1 };
    const int threadCount = 32;
    const int incrementsPerThread = 1000;

    std::vector<std::thread> threads;
    for (int i = 0; i < threadCount; ++i) {
        threads.emplace_back([this, &key, incrementsPerThread] {
            for (int j = 0; j < incrementsPerThread; ++j) {
                this->table.increment(key);
            }
        });
    }

    for (auto& t : threads) t.join();

    auto data = this->ExportData();
    EXPECT_EQ(data[key], static_cast<Int>(threadCount * incrementsPerThread));
    EXPECT_EQ(this->table.load(), 1);
}

// Concurrent increments on many keys
TYPED_TEST(DynamicTableTest, ConcurrentIncrementsManyKeys) {
    using Key = typename TypeParam::first_type;
    using Int = typename TypeParam::second_type;

    const int threadCount = 64;
    const int keysPerThread = 100;

    std::vector<std::thread> threads;
    for (int i = 0; i < threadCount; ++i) {
        threads.emplace_back([this, i, keysPerThread] {
            for (int j = 0; j < keysPerThread; ++j) {
                Key key = Key{ i * keysPerThread + j };
                this->table.increment(key);
            }
        });
    }

    for (auto& t : threads) t.join();

    EXPECT_EQ(this->table.load(), threadCount * keysPerThread);

    auto data = this->ExportData();
    for (const auto& [key, count] : data) {
        EXPECT_EQ(count, static_cast<Int>(1));
    }
}

// Test reallocation under load
TYPED_TEST(DynamicTableTest, ReallocationConsistency) {
    using Key = typename TypeParam::first_type;
    using Int = typename TypeParam::second_type;

    const int threadCount = 32;
    const int keysPerThread = this->initialCapacity;

    std::vector<std::thread> threads;
    for (int i = 0; i < threadCount; ++i) {
        threads.emplace_back([this, i, keysPerThread] {
            for (int j = 0; j < keysPerThread; ++j) {
                Key key = Key{ i * keysPerThread + j };
                this->table.increment(key);
            }
        });
    }

    for (auto& t : threads) t.join();

    auto data = this->ExportData();
    EXPECT_EQ(data.size(), threadCount * keysPerThread);
    EXPECT_GE(this->table.capacity(), this->table.load());

    for (const auto& [key, count] : data) {
        EXPECT_EQ(count, static_cast<Int>(1));
    }
}

// LazyMap tests

// Sample Container class
struct HeavyObject {
    std::string name;
    int value;
    HeavyObject(std::string n, int v) : name(std::move(n)), value(v) {}
};

// Basic test with int key and HeavyObject container
TEST(LazyMapTest, BasicConstructionAndRetrieval) {
    pmap::LazyMap<int, HeavyObject> map;
    auto [ptr, created] = map.get(1, "Test", 42);
    ASSERT_TRUE(created);
    ASSERT_NE(ptr, nullptr);
    EXPECT_EQ(ptr->name, "Test");
    EXPECT_EQ(ptr->value, 42);

    auto [ptr2, created2] = map.get(1, "Ignored", 0);
    ASSERT_FALSE(created2);
    EXPECT_EQ(ptr2, ptr);
}

TEST(LazyMapTest, StringKeyVectorContainer) {
    pmap::LazyMap<std::string, std::vector<int>> map;
    auto [ptr, created] = map.get("key1", std::initializer_list<int>{1, 2, 3});
    ASSERT_TRUE(created);
    ASSERT_EQ(ptr->size(), 3);
    EXPECT_EQ((*ptr)[0], 1);
}

struct ComplexObject {
    std::string label;
    std::vector<int> data;
    ComplexObject(std::string l, std::vector<int> d) : label(std::move(l)), data(std::move(d)) {}
};

TEST(LazyMapTest, ForwardingArguments) {
    pmap::LazyMap<int, ComplexObject> map;
    std::vector<int> vec = {10, 20, 30};
    auto [ptr, created] = map.get(5, "Label", vec);
    ASSERT_TRUE(created);
    EXPECT_EQ(ptr->label, "Label");
    EXPECT_EQ(ptr->data.size(), 3);
}

TEST(LazyMapTest, ConcurrentAccess) {
    pmap::LazyMap<int, HeavyObject> map;
    const int key = 99;
    std::atomic<int> created_count{0};

    auto thread_func = [&]{
        auto [ptr, created] = map.get(key, "Concurrent", 100);
        if (created) created_count++;
        ASSERT_NE(ptr, nullptr);
        EXPECT_EQ(ptr->name, "Concurrent");
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < 10; ++i)
        threads.emplace_back(thread_func);

    for (auto& t : threads)
        t.join();

    EXPECT_EQ(created_count.load(), 1);
}

TEST(LazyMapTest, HighLoadMultipleKeys) {
    pmap::LazyMap<int, HeavyObject> map;
    const int num_keys = 1000;

    for (int i = 0; i < num_keys; ++i) {
        auto [ptr, created] = map.get(i, "Obj" + std::to_string(i), i);
        ASSERT_TRUE(created);
        EXPECT_EQ(ptr->value, i);
    }

    for (int i = 0; i < num_keys; ++i) {
        auto [ptr, created] = map.get(i, "ShouldNotCreate", -1);
        ASSERT_FALSE(created);
        EXPECT_EQ(ptr->value, i);
    }
}
