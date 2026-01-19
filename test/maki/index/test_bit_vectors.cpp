#include <atomic>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <thread>
#include <unordered_set>
#include <vector>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <sdsl/vectors.hpp>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_invoke.h>
#include <maki/bitvectors.hpp>

std::vector<size_t> generateRandomIntegers(size_t low_, size_t high_, size_t n_) {
    std::vector<size_t> _data(n_);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> unif(low_, high_);
    for (size_t i = 0; i < n_; ++i)
        _data[i] = unif(gen);

    return _data;
}

std::vector<size_t> generateShuffledIndex(size_t n_) {
    std::vector<size_t> _data(n_);
    std::iota(_data.begin(), _data.end(), (size_t)0u);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(_data.begin(), _data.end(), gen);

    return _data;
}

template <uint8_t t_width>
std::vector<size_t> exportValues(sdsl::int_vector<t_width> const &vec) {
    std::vector<size_t> result(vec.size());

    for (size_t i = 0; i < vec.size(); ++i)
        result[i] = vec[i];

    return result;
}

template <uint8_t t_width>
bool sawSDSLDataRace(size_t n, size_t gsize) {
    auto vals = generateRandomIntegers(0, (1<<t_width)-1, n);
    auto index = generateShuffledIndex(n);

    sdsl::int_vector<t_width> dest(n);
    oneapi::tbb::parallel_for(
        (size_t)0, n, gsize,
        [&](size_t i)->void {
            dest[index[i]] = vals[index[i]];
        }
    );
    auto result = exportValues(dest);

    return result != vals;
}

template <uint8_t t_width>
bool sawThreadedDataRace(size_t n, size_t samplerate, size_t gsize) {
    auto vals = generateRandomIntegers(0, (1<<t_width)-1, n);
    auto index = generateShuffledIndex(n);

    bv::threaded::IntVector<t_width> dest(n, samplerate);

    oneapi::tbb::parallel_for(
        (size_t)0, n, gsize,
        [&](size_t i)->void {
            typename bv::threaded::IntVector<t_width>::Accessor p(dest);
            p[index[i]] = vals[index[i]];
        }
    );
    typename bv::threaded::IntVector<t_width>::ConstAccessor cp(dest);
    auto result = exportValues<t_width>(cp);

    return result != vals;
}

class MutexManagerTests : public testing::Test {
};

TEST_F(MutexManagerTests, Allocate) {
    bv::detail::MutexManager mc(100);
    auto *p = mc.raw();
    for (size_t i = 0; i < (mc.capacity() / 8); ++i)
        ASSERT_EQ((p+i)->load(), 0);
}

TEST_F(MutexManagerTests, TryLock) {
    bv::detail::MutexManager lockMgr(1);
    lockMgr.lock(0);
    ASSERT_FALSE(lockMgr.tryLock(0));
    ASSERT_FALSE(lockMgr.tryLock(0));
    ASSERT_FALSE(lockMgr.tryLock(0));
    lockMgr.unlock(0);
    ASSERT_EQ(lockMgr.raw()->load(), 0);
}

#pragma GCC push_options
#pragma GCC optimize ("O0")
TEST_F(MutexManagerTests, AtomicIncrement) {
    bv::detail::MutexManager mc(100);
    size_t count = 0,
           upto = 10000,
           numthreads = 20;

    std::vector<std::thread> threads;
    for (size_t t = 0; t < numthreads; ++t)
        threads.emplace_back(
            [&]->void {
                for (size_t i = 0; i < upto; ++i) {
                    mc.lock(0);
                    ++count;
                    mc.unlock(0);
                }
            }
        );
    for (auto &thread : threads)
        thread.join();

    ASSERT_EQ(count, upto * numthreads);
}
#pragma GCC pop_options

TEST_F(MutexManagerTests, SerialAllPositions) {
    std::vector<size_t> data(100),
                        exp(100, 1);

    bv::detail::MutexManager mc(100);
    for (size_t i = 0; i < 100; ++i) {
        mc.lock(i);
        data[i] = 1;
        mc.unlock(i);
    }

    auto *p = mc.raw();
    for (size_t i = 0; i < (mc.capacity() / 8); ++i)
        ASSERT_EQ((p+i)->load(), 0);

    ASSERT_EQ(data, exp);
}

TEST_F(MutexManagerTests, ParallelAllPositions) {
    size_t n = 1000;
    std::vector<size_t> data(n),
                        exp(n, 1);

    bv::detail::MutexManager mc(n);
    oneapi::tbb::parallel_for(
        (size_t)0, n, (size_t)1,
        [&](size_t i)->void {
            mc.lock(i);
            data[i] = 1;
            mc.unlock(i);
        }
    );

    auto *p = mc.raw();
    for (size_t i = 0; i < (mc.capacity() / 8); ++i)
        ASSERT_EQ((p+i)->load(), 0);

    ASSERT_EQ(data, exp);
}

class IntVectorTests : public testing::Test {
protected:
    IntVectorTests()
    {
    }
};

TEST_F(IntVectorTests, OneByte) {
    size_t numelems = 8;
    bv::threaded::IntVector<1> v(numelems, 64);

    oneapi::tbb::parallel_for(
        (size_t)0, numelems, (size_t)1,
        [&](size_t i)->void {
            bv::threaded::IntVector<1>::Accessor p(v);
            p[i] = 1;
        }
    );

    bv::threaded::IntVector<1>::ConstAccessor p(v);
    for (size_t i = 0; i < numelems; ++i)
        ASSERT_EQ(p[i], 1);
}

TEST_F(IntVectorTests, ManyBytes) {
    size_t numelems = 20000;
    bv::threaded::IntVector<1> v(numelems, 64);

    oneapi::tbb::parallel_for(
        (size_t)0, numelems, (size_t)1,
        [&](size_t i)->void {
            bv::threaded::IntVector<1>::Accessor p(v);
            p[i] = 1;
        }
    );

    bv::threaded::IntVector<1>::ConstAccessor p(v);
    for (size_t i = 0; i < numelems; ++i)
        ASSERT_EQ(p[i], 1);
}

class BVDataRaceTests : public testing::Test {
protected:
    BVDataRaceTests() : tries(100) {}
    size_t tries;
};

TEST_F(BVDataRaceTests, SDSLw1bv100) {
    size_t errors = 0;
    for (size_t i = 0; i < tries; ++i)
        if (sawSDSLDataRace<1>(100, 1))
            ++errors;
    std::cout << "Saw " << errors << " data races in " << tries << " attempts" << std::endl;
}

TEST_F(BVDataRaceTests, Threadedw1bv100) {
    size_t errors = 0;
    for (size_t i = 0; i < tries; ++i)
        if (sawThreadedDataRace<1>(100, 64, 1))
            ++errors;
    ASSERT_EQ(errors, 0);
}

TEST_F(BVDataRaceTests, SDSLw4bv100) {
    size_t errors = 0;
    for (size_t i = 0; i < tries; ++i)
        if (sawSDSLDataRace<4>(100, 1))
            ++errors;
    std::cout << "Saw " << errors << " data races in " << tries << " attempts" << std::endl;
}

TEST_F(BVDataRaceTests, Threadedw4bv100) {
    size_t errors = 0;
    for (size_t i = 0; i < tries; ++i)
        if (sawThreadedDataRace<4>(100, 16, 1))
            ++errors;
    ASSERT_EQ(errors, 0);
}

TEST_F(BVDataRaceTests, SDSLw1bv10000) {
    size_t errors = 0;
    for (size_t i = 0; i < tries; ++i)
        if (sawSDSLDataRace<1>(10000, 1))
            ++errors;
    std::cout << "Saw " << errors << " data races in " << tries << " attempts" << std::endl;
}

TEST_F(BVDataRaceTests, Threadedw1bv10000) {
    size_t errors = 0;
    for (size_t i = 0; i < tries; ++i)
        if (sawThreadedDataRace<1>(10000, 64, 1))
            ++errors;
    ASSERT_EQ(errors, 0);
}

TEST_F(BVDataRaceTests, SDSLw4bv10000) {
    size_t errors = 0;
    for (size_t i = 0; i < tries; ++i)
        if (sawSDSLDataRace<4>(10000, 1))
            ++errors;
    std::cout << "Saw " << errors << " data races in " << tries << " attempts" << std::endl;
}

TEST_F(BVDataRaceTests, Threadedw4bv10000) {
    size_t errors = 0;
    for (size_t i = 0; i < tries; ++i)
        if (sawThreadedDataRace<4>(10000, 16, 1))
            ++errors;
    ASSERT_EQ(errors, 0);
}

// Concurrent Bit Vector Tests

// Helper to run multiple threads and collect results
void runConcurrentSet(bv::threaded::ConcurrentBitVector& bv, size_t index, std::atomic<int>& successCount, int threadCount) {
    std::vector<std::thread> threads;
    for (int i = 0; i < threadCount; ++i) {
        threads.emplace_back([&bv, index, &successCount] {
            if (bv.trySet(index)) {
                successCount++;
            }
        });
    }
    for (auto& t : threads) t.join();
}

// Test with small vector size
TEST(ConcurrentBitVectorTest, SmallVectorSingleThread) {
    bv::threaded::ConcurrentBitVector bv(8);
    EXPECT_TRUE(bv.trySet(3));
    EXPECT_FALSE(bv.trySet(3));
    EXPECT_TRUE(bv.trySet(7));
    EXPECT_FALSE(bv.trySet(7));
}

// Test with large vector size
TEST(ConcurrentBitVectorTest, LargeVectorSingleThread) {
    bv::threaded::ConcurrentBitVector bv(10000);
    EXPECT_TRUE(bv.trySet(9999));
    EXPECT_FALSE(bv.trySet(9999));
    EXPECT_TRUE(bv.trySet(0));
    EXPECT_FALSE(bv.trySet(0));
}

// Test concurrent access to the same bit
TEST(ConcurrentBitVectorTest, ConcurrentSetSameBit) {
    bv::threaded::ConcurrentBitVector bv(64);
    std::atomic<int> successCount(0);
    runConcurrentSet(bv, 42, successCount, 100);
    EXPECT_EQ(successCount.load(), 1);  // Only one thread should succeed
}

// Test concurrent access to different bits
TEST(ConcurrentBitVectorTest, ConcurrentSetDifferentBits) {
    const size_t size = 100;
    bv::threaded::ConcurrentBitVector bv(size);
    std::vector<std::thread> threads;
    std::atomic<int> successCount(0);

    for (size_t i = 0; i < size; ++i) {
        threads.emplace_back([&bv, i, &successCount] {
            if (bv.trySet(i)) {
                successCount++;
            }
        });
    }

    for (auto& t : threads) t.join();
    EXPECT_EQ(successCount.load(), size);  // All bits should be set successfully
}

// Stress test with high thread count
TEST(ConcurrentBitVectorTest, StressTestHighThreadCount) {
    const size_t size = 1000;
    bv::threaded::ConcurrentBitVector bv(size);
    std::atomic<int> successCount(0);
    std::vector<std::thread> threads;

    for (size_t i = 0; i < size * 10; ++i) {
        threads.emplace_back([&bv, i, &successCount, size] {
            size_t index = i % size;
            if (bv.trySet(index)) {
                successCount++;
            }
        });
    }

    for (auto& t : threads) t.join();
    EXPECT_EQ(successCount.load(), size);  // Each bit should be set only once
}

// Random access in single-threaded context
TEST(ConcurrentBitVectorTest, RandomAccessSingleThread) {
    const size_t size = 1000;
    bv::threaded::ConcurrentBitVector bv(size);
    std::mt19937 rng(42);  // Fixed seed for reproducibility
    std::uniform_int_distribution<size_t> dist(0, size - 1);

    std::unordered_set<size_t> accessed;

    for (int i = 0; i < 500; ++i) {
        size_t index = dist(rng);
        bool result = bv.trySet(index);
        if (accessed.count(index) == 0) {
            EXPECT_TRUE(result);
            accessed.insert(index);
        } else {
            EXPECT_FALSE(result);
        }
    }
}

// Random access in multi-threaded context
TEST(ConcurrentBitVectorTest, RandomAccessMultiThread) {
    const size_t size = 1000;
    bv::threaded::ConcurrentBitVector bv(size);
    std::mt19937 rng(123);  // Fixed seed
    std::uniform_int_distribution<size_t> dist(0, size - 1);

    std::atomic<int> successCount(0);
    std::vector<std::thread> threads;

    for (int i = 0; i < 1000; ++i) {
        threads.emplace_back([&bv, &dist, &rng, &successCount] {
            static thread_local std::mt19937 local_rng(rng());
            size_t index = dist(local_rng);
            if (bv.trySet(index)) {
                successCount++;
            }
        });
    }

    for (auto& t : threads) t.join();

    // Success count should be <= size since each bit can only be set once
    EXPECT_LE(successCount.load(), size);
}

class RandomConcurrentBitVectorTests : public testing::Test
{
protected:
    RandomConcurrentBitVectorTests()
    {
    }

    static std::vector<uint64_t> generateBitVector(size_t nbits_, std::vector<size_t> const &setBits_)
    {
        std::vector<uint64_t> _data((nbits_ + 63) / 64, 0u);

        for (size_t const &b_ : setBits_)
            _data[b_ / 64] |= (size_t)1u << (b_ % 64);

        return _data;
    }

    static inline std::vector<size_t> generateSetBits(size_t low_, size_t high_, size_t n_)
    {
        auto raw = generateRandomIntegers(low_, high_, n_);
        std::sort(raw.begin(), raw.end());
        auto end = std::unique(raw.begin(), raw.end());
        raw.erase(end, raw.end());
        return raw;
    }

    static std::vector<size_t> shuffleIota(size_t n_)
    {
        std::vector<size_t> _data(n_);
        std::iota(_data.begin(), _data.end(), (size_t)0u);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(_data.begin(), _data.end(), gen);

        return _data;
    }

    static std::vector<size_t> shuffleIndex(std::vector<size_t> const &src)
    {
        size_t n_ = src.size();
        std::vector<size_t> _dest(n_);

        std::vector<size_t> newInds = shuffleIota(n_);
        for (size_t i = 0; i < n_; ++i)
            _dest[newInds[i]] = src[i];

        return _dest;
    }
};

TEST_F(RandomConcurrentBitVectorTests, trySet)
{
    size_t bitvecSize = 1000, density = 900;

    // generate bits to set

    auto bitIndices = generateSetBits(0, bitvecSize - 1, density);
    std::atomic_uint64_t bitCounts = 0u;

    // truth

    auto expected = generateBitVector(bitvecSize, bitIndices);

    // experiment

    bv::threaded::ConcurrentBitVector bits(bitvecSize);

    auto f = [&] {
        auto takeOrder = shuffleIota(bitIndices.size());
        for (size_t &i : takeOrder)
            if (bits.trySet(bitIndices[i]))
                ++bitCounts;
    };

    oneapi::tbb::parallel_invoke(f, f, f, f, f, f, f, f, f, f);

    std::vector<uint64_t> actual;
    actual.reserve(bits.raw().size());
    for (copyable_atomic_uint64_t const &v : bits.raw())
        actual.push_back(static_cast<uint64_t>(v));

    // check

    ASSERT_EQ(actual, expected);
    ASSERT_EQ(bitCounts.load(), bitIndices.size());
}
