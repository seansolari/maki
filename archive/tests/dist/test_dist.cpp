#include <algorithm>
#include <atomic>
#include <chrono>
#include <iostream>
#include <filesystem>
#include <numeric>
#include <random>
#include <unordered_map>
#include <vector>
#include <glog/logging.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <sdsl/vectors.hpp>
#include <maki/fasta.hpp>
#include <maki/graph.hpp>
#include <maki/dist.hpp>
#include <maki/utils.hpp>

namespace fs = std::filesystem;

std::vector<size_t> generateRandomIntegers(size_t low_, size_t high_, size_t n_)
{
    std::vector<size_t> _data;
    _data.reserve(n_);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> unif(low_, high_);
    for (size_t i = 0; i < n_; ++i)
        _data.push_back(unif(gen));

    std::sort(_data.begin(), _data.end());
    auto _end = std::unique(_data.begin(), _data.end());
    _data.resize(_end - _data.begin());

    return _data;
}

graph::WriteableGraph makeGraph(const Dna4GenomeVector &data, uint8_t k, fs::path bufferPath) {
    graph::WriteableGraph tempgraph = graph::initialiseEmptyGraph(
        data,
        k,
        detail::calculateNumEdgesBound(data),
        bufferPath);
    graph::makeGraphBulk(tempgraph, data, 1);
    return tempgraph;
}

class RandomKmerPointersTests : public testing::Test
{
protected:
    RandomKmerPointersTests()
        : k(31),
          rd(),
          dna(fasta_utils::randomGenomes(5, rd(), k, 100)),
          kmers(k),
          tempTestDir(tempio::create_temporary_directory()),
          dbg(makeGraph(dna, k, tempio::create_temporary_file(tempTestDir, ".h5")))
    {
        for (auto const &g : dna)
        {
            for (auto const &ctg : g.contigs)
            {
                for (auto const &rec : ctg.sequences.first)
                {
                    kmers.push_back(rec.sequence);
                }
            }
        }

        // sort k-mers

        kmers.sort();
        kmers.uniq();
        kmers.remove_null_terminals();
    }

    ~RandomKmerPointersTests()
    {
        LOG(INFO) << "removing temporary directory " << tempTestDir;
        fs::remove_all(tempTestDir);
    }

    size_t k;

    std::random_device rd;
    Dna4GenomeVector dna;

    mers::detail::_Kmers kmers;

    fs::path tempTestDir;
    graph::TraversableGraph dbg;
};

TEST_F(RandomKmerPointersTests, IsSorted)
{
    auto
        it = kmers.begin(),
        next = it + 1,
        end = kmers.end();

    mers::detail::_Kmers_LessThan LT{};

    while (next != end)
    {
        ASSERT_FALSE(LT(*next, *it));
        ++it;
        ++next;
    }
}

TEST_F(RandomKmerPointersTests, EdgesSize)
{
    ASSERT_EQ(kmers.size(), dbg.last.size());
}

TEST_F(RandomKmerPointersTests, CheckFoward)
{
    mers::detail::_Kmers_FwdEq FE{};

    for (size_t i = 0; i < kmers.size(); ++i)
    {
        auto optj = dbg.fwd(i);
        if (optj)
        {
            ASSERT_TRUE(FE(kmers[i], kmers[*optj]));
        }
        else
        {
            ASSERT_EQ(dbg.edge(i), 0u);
        }
    }
}

TEST_F(RandomKmerPointersTests, LoadSave)
{
    fs::path tmp = tempio::create_temporary_directory();
    graph::saveToDisk(dbg, tmp);

    graph::TraversableGraph newdbg;
    graph::loadFromDisk(tmp, newdbg);

    // check 20 random k-mers

    mers::detail::_Kmers_FwdEq FE{};

    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> len(0, kmers.size() - 1);

    for (size_t i = 0; i < 20; ++i)
    {
        size_t j = len(gen);
        auto optj = dbg.fwd(j);
        if (optj)
        {
            ASSERT_TRUE(FE(kmers[j], kmers[*optj]));
        }
        else
        {
            ASSERT_EQ(dbg.edge(j), 0u);
        }
    }

    LOG(INFO) << "removing temporary directory " << tmp;
    fs::remove_all(tmp);
}

TEST_F(RandomKmerPointersTests, GetTerminals)
{
    kmers.uniq_nodes();
    sdsl::bit_vector arr = dbg.getTerminals();

    ASSERT_EQ(kmers.size(), arr.size());
}

TEST_F(RandomKmerPointersTests, TerminalRanges)
{
    auto gen = dbg.nonTerminalNodes();

    size_t i = 0;
    while (auto res = gen.next())
    {
        // up to this range should be terminals

        while (i < res->begin)
        {
            ASSERT_LT(kmers[i++].sequence.size(), k);
        }

        // this range should contain no terminals

        while (i < res->end)
        {
            ASSERT_EQ(kmers[i++].sequence.size(), k);
        }
    }

    // remaining should be terminals

    while (i < kmers.size())
    {
        ASSERT_LT(kmers[i++].sequence.size(), k);
    }
}

TEST_F(RandomKmerPointersTests, TerminalWholeSubrange)
{
    auto gen = dbg.nonTerminalNodes().subrange(0, kmers.size());

    size_t i = 0;
    while (auto res = gen.next())
    {
        // up to this range should be terminals

        while (i < res->begin)
        {
            ASSERT_LT(kmers[i++].sequence.size(), k);
        }

        // this range should contain no terminals

        while (i < res->end)
        {
            ASSERT_EQ(kmers[i++].sequence.size(), k);
        }
    }

    // remaining should be terminals

    while (i < kmers.size())
    {
        ASSERT_LT(kmers[i++].sequence.size(), k);
    }
}

TEST_F(RandomKmerPointersTests, TerminalSubrange)
{
    auto gen = dbg.nonTerminalNodes();

    std::mt19937 f(rd());
    std::uniform_int_distribution<> dist(0, kmers.size());

    size_t
        begin = dist(f),
        end = dist(f);

    while (begin == end)
        end = dist(f);
    if (begin > end)
        std::swap(begin, end);
    auto subgen = gen.subrange(begin, end);

    size_t i = begin;
    while (auto res = subgen.next())
    {
        // up to this range should be terminals

        while (i < res->begin)
        {
            ASSERT_LT(kmers[i++].sequence.size(), k);
        }

        // this range should contain no terminals

        while (i < res->end)
        {
            ASSERT_EQ(kmers[i++].sequence.size(), k);
        }
    }

    // remaining should be terminals

    while (i < end)
    {
        ASSERT_LT(kmers[i++].sequence.size(), k);
    }
}

class DynamicMapTests : public testing::Test {
protected:
    using Cell = pmap::CellIndex<15>;
    using Dmap = pmap::DynamicCellTable<15, uint32_t>;

    DynamicMapTests()
        : _tempBase(tempio::create_temporary_directory())
        , _alloc(pmap::alloc::makeDynamicAllocator(64 /* bytes */, 11 /* 2048 bytes max */, _tempBase))
        , _map(_alloc)
    {
    }

    ~DynamicMapTests()
    { std::filesystem::remove_all(_tempBase); }

    std::string                   _tempBase;
    pmap::alloc::DynamicAllocator _alloc;
    Dmap                          _map;
};

TEST_F(DynamicMapTests, minBuckets) {
    ASSERT_EQ(Dmap::pageSizeLog2(Dmap::items(Dmap::initialAllocationBytes/2)),
              Dmap::initialAllocationL2);
}

TEST_F(DynamicMapTests, RoundBuckets) {
    // 8 -> 10
    constexpr size_t bytes = (size_t)1u << 8;
    static_assert(Dmap::items(bytes>>1) < Dmap::items(bytes));
    ASSERT_EQ(Dmap::pageSizeLog2(Dmap::items(bytes)),
              10);
}

TEST_F(DynamicMapTests, InitialSize) {
    ASSERT_EQ(_map.load(),      0);
    ASSERT_EQ(_map.capacity(),  12);
}

TEST_F(DynamicMapTests, IncrementNoRealloc) {
    std::vector<pmap::BlockIndex> occs = {
        pmap::BlockIndex{ 2, 0 }, pmap::BlockIndex{ 5, 0 },
        pmap::BlockIndex{ 0, 0 }, pmap::BlockIndex{ 0, 4 },
        pmap::BlockIndex{ 1, 0 }, pmap::BlockIndex{ 1, 0 },
        pmap::BlockIndex{ 3, 2 }, pmap::BlockIndex{ 2, 0 },
        pmap::BlockIndex{ 0, 4 }, pmap::BlockIndex{ 2, 0 },
        pmap::BlockIndex{ 2, 2 }, pmap::BlockIndex{ 3, 3 },
        pmap::BlockIndex{ 1, 5 }, pmap::BlockIndex{ 2, 5 },
        pmap::BlockIndex{ 5, 0 }, pmap::BlockIndex{ 0, 5 },
        pmap::BlockIndex{ 2, 4 }
    };

    // truth
    std::unordered_map<copyable_pair<size_t>,uint32_t,copyable_pair_hash_32<size_t>> exp, actual;
    for (pmap::BlockIndex const &v : occs) ++exp[v];

    // actual
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<std::vector<pmap::BlockIndex>::const_iterator>(occs.cbegin(), occs.cend(), 10),
        [&](oneapi::tbb::blocked_range<std::vector<pmap::BlockIndex>::const_iterator> const &r)->void {
            for (auto it = r.begin(); it != r.end(); ++it)
                _map.increment(*it);
        }
    );

    _map.forEach(
        [&](const Cell &k, uint32_t z)->void {
            actual[k.asPair<size_t>()] = z;
        }
    );

    ASSERT_EQ(exp, actual);

    ASSERT_EQ(_map.load(),      12);
    ASSERT_EQ(_map.capacity(),  12);
}

TEST_F(DynamicMapTests, IncrementRealloc) {
    std::vector<pmap::BlockIndex> occs = {
        pmap::BlockIndex{ 2, 0 }, pmap::BlockIndex{ 5, 0 },
        pmap::BlockIndex{ 0, 0 }, pmap::BlockIndex{ 0, 4 },
        pmap::BlockIndex{ 1, 0 }, pmap::BlockIndex{ 1, 0 },
        pmap::BlockIndex{ 3, 2 }, pmap::BlockIndex{ 2, 0 },
        pmap::BlockIndex{ 0, 4 }, pmap::BlockIndex{ 2, 0 },
        pmap::BlockIndex{ 2, 2 }, pmap::BlockIndex{ 3, 3 },
        pmap::BlockIndex{ 1, 5 }, pmap::BlockIndex{ 2, 5 },
        pmap::BlockIndex{ 5, 0 }, pmap::BlockIndex{ 0, 5 },
        pmap::BlockIndex{ 2, 4 }, pmap::BlockIndex{ 32767, 32767 }
    };

    // truth
    std::unordered_map<copyable_pair<size_t>,uint32_t,copyable_pair_hash_32<size_t>> exp, actual;
    for (pmap::BlockIndex const &v : occs) ++exp[v];

    // actual
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<std::vector<pmap::BlockIndex>::const_iterator>(occs.cbegin(), occs.cend(), 10),
        [&](oneapi::tbb::blocked_range<std::vector<pmap::BlockIndex>::const_iterator> const &r)->void {
            for (auto it = r.begin(); it != r.end(); ++it)
                _map.increment(*it);
        }
    );

    _map.forEach(
        [&](const Cell &k, uint32_t z)->void {
            actual[k.asPair<size_t>()] = z;
        }
    );

    ASSERT_EQ(exp, actual);

    ASSERT_EQ(_map.load(),      13);
    ASSERT_EQ(_map.capacity(),  90);
}

// LazyGrid Tests

TEST(GridTests, CreateGrid) {
    size_t n = 1000;
    pmap::LazyGrid<2, copyable_atomic_uint64_t> grid(n);
}

TEST(GridTests, GridAllocate) {
    pmap::LazyGrid<2, copyable_atomic_uint64_t> grid(501);
    std::pair<copyable_atomic_uint64_t*,bool> qres = grid.get(grid.coord(500, 500).block);
    ASSERT_TRUE(qres.second);
    ASSERT_EQ(qres.first->load(), 0u);
}

TEST(GridTests, GridParallelIncrement) {
    pmap::LazyGrid<15, copyable_atomic_uint64_t> grid(501);
    std::atomic_uint64_t numAllocs = 0;

    auto f = [&] {
        std::pair<copyable_atomic_uint64_t*,bool> qres = grid.get(grid.coord(500, 500).block);
        if (qres.second)
            ++numAllocs;
        ++(*qres.first);
    };
    oneapi::tbb::parallel_invoke(f, f, f, f, f, f, f, f, f, f);

    ASSERT_EQ(numAllocs, 1);

    std::pair<copyable_atomic_uint64_t*,bool> qres = grid.get(grid.coord(500, 500).block);
    ASSERT_FALSE(qres.second);
    ASSERT_EQ(qres.first->load(), 10);
}

/* PRACTICAL TESTING */

class RandomKmerDistTests : public testing::Test {
protected:
    RandomKmerDistTests()
        : k(31)
        , numSequences(5)
        , rd()
        , templateSequence(fasta_utils::randomSequence(k, 1000, rd()))
        , dna(fasta_utils::randomMutations(templateSequence, numSequences, 0.01, rd()))
        , kmers()
        , tempTestDir(tempio::create_temporary_directory())
        , dbg(makeGraph(dna, k, tempio::create_temporary_file(tempTestDir, ".h5")))
    {
        for (auto const& g : dna)
        {
            for (auto const& ctg : g.contigs)
            {
                for (auto const& rec : ctg.sequences.first)
                {
                    // create k-mer store
                    auto &obj = kmers.emplace_back(k);

                    // insert sequence
                    obj.push_back(rec.sequence);

                    // sort k-mers
                    obj.force_remove_terminals();
                    obj.sort();
                    obj.uniq();
                }
            }
        }
    }

    ~RandomKmerDistTests()
    {
        LOG(INFO) << "removing temporary directory " << tempTestDir;
        fs::remove_all(tempTestDir);
    }

    size_t k, numSequences;

    std::random_device rd;
    Dna4Sequence templateSequence;
    Dna4GenomeVector dna;

    std::vector<mers::detail::_Kmers> kmers;

    fs::path tempTestDir;
    graph::TraversableGraph dbg;
};

TEST_F(RandomKmerDistTests, DensePairwiseDistances)
{
    // compute pairwise distances using graph
    auto m = dist::pairwise::countPairsDense(dbg);

    std::map<std::pair<size_t, size_t >, size_t > res;

    size_t N = m.size();
    for (size_t r = 2; r < N; ++r)
    {
        for (size_t c = 1; c < r; ++c)
        {
            auto v = (uint32_t)m(r, c);
            if (v > 0)
                res[std::make_pair(c, r)] = v;
        }
    }
    // compute pairwise distances using arrays

    auto z = mers::detail::_pairwise(kmers);

    ASSERT_EQ(res, z);
}

TEST_F(RandomKmerDistTests, SparsePairwiseDistances)
{
    // compute pairwise distances using graph

    auto m = dist::pairwise::countPairsSparse<2,uint32_t>(dbg, tempTestDir);

    std::map<std::pair<size_t,size_t>,size_t> res;
    m.forEach([&](size_t a, size_t b, uint32_t count) {
        res[std::make_pair(a, b)] = count;
    });

    // compute pairwise distances using arrays

    auto z = mers::detail::_pairwise(kmers);

    ASSERT_EQ(res, z);
}

// Xlist tests

class XlistTest : public ::testing::Test {
protected:
    std::mt19937 rng{42};

    std::vector<uint64_t> GenerateTuple(int size, int bit_width) {
        std::uniform_int_distribution<uint64_t> dist(0, (bit_width == 64) ? UINT64_MAX : ((1LL << bit_width) - 1));
        std::vector<uint64_t> tuple(size);
        for (uint64_t& val : tuple)
            val = dist(rng);
        return tuple;
    }

    std::vector<std::pair<std::vector<uint64_t>, uint32_t>> ExportData(dist::cooc::Xlist& xlist) const {
        std::vector<std::pair<std::vector<uint64_t>, uint32_t>> result;
        xlist.forEach(
            [&](sdsl::int_vector<0>::const_iterator begin,
                sdsl::int_vector<0>::const_iterator end,
                uint32_t count)->void
        {
            result.emplace_back(std::make_pair(std::vector(begin, end), count));
        });
        return result;
    }

    void runBasicInsertTest(size_t bit_width, size_t tuple_size) {
        dist::cooc::Xlist xlist(bit_width, tuple_size);

        std::map<std::vector<uint64_t>, uint32_t> expected_counts;
        // Insert 5 unique tuples, each with a random number of repetitions (1–3)
        for (int i = 0; i < 5; ++i) {
            auto tuple = GenerateTuple(tuple_size, bit_width);

            int repetitions = 1 + (rng() % 3);
            for (int j = 0; j < repetitions; ++j)
                xlist.increment(tuple.begin(), tuple.end());

            expected_counts[tuple] += repetitions;
        }

        auto data = ExportData(xlist);
        ASSERT_EQ(data.size(), expected_counts.size());

        for (const auto& [tuple, count] : data) {
            auto it = expected_counts.find(tuple);
            ASSERT_NE(it, expected_counts.end()) << "Unexpected tuple found in export.";
            EXPECT_EQ(count, it->second) << "Count mismatch for tuple.";
        }
    }

    void benchmarkInsert(dist::cooc::Xlist& xlist, int num_tuples, int tuple_size, int bit_width) {
        std::uniform_int_distribution<int> dist(0, (1 << bit_width) - 1);
        std::vector<int> tuple(tuple_size);

        std::chrono::duration<double> elapsed(0);
        std::chrono::system_clock::time_point start, end;

        for (int i = 0; i < num_tuples; ++i) {
            for (int& val : tuple) val = dist(rng);
            start = std::chrono::high_resolution_clock::now();
            xlist.increment(tuple.begin(), tuple.end());
            end = std::chrono::high_resolution_clock::now();
            elapsed += end - start;
        }

        std::cout << "Inserted " << num_tuples << " tuples of size " << tuple_size
                  << " (bit width " << bit_width << ") in " << elapsed.count() << " seconds.\n";
    }
};

// --- Basic Functionality Tests ---

TEST_F(XlistTest, CountSingleTuple) {
    dist::cooc::Xlist xlist(3, 3);
    std::vector<uint64_t> tuple = {1, 2, 3};
    xlist.increment(tuple.begin(), tuple.end());

    auto data = ExportData(xlist);
    ASSERT_EQ(data.size(), 1);
    EXPECT_EQ(data[0].first, tuple);
    EXPECT_EQ(data[0].second, 1);
}

TEST_F(XlistTest, CountMultipleSameTuples) {
    dist::cooc::Xlist xlist(3, 3);
    std::vector<uint64_t> tuple = {1, 2, 3};
    for (int i = 0; i < 5; ++i)
        xlist.increment(tuple.begin(), tuple.end());

    auto data = ExportData(xlist);
    ASSERT_EQ(data.size(), 1);
    EXPECT_EQ(data[0].first, tuple);
    EXPECT_EQ(data[0].second, 5);
}

TEST_F(XlistTest, CountDifferentTuples) {
    dist::cooc::Xlist xlist(3, 3);
    std::vector<uint64_t> t1 = {1, 2, 3};
    std::vector<uint64_t> t2 = {4, 5, 6};
    xlist.increment(t1.begin(), t1.end());
    xlist.increment(t2.begin(), t2.end());

    auto data = ExportData(xlist);
    ASSERT_EQ(data.size(), 2);
}

// Macro to generate tests for each combination
#define GENERATE_BASIC_TEST(bit_width, tuple_size) \
    TEST_F(XlistTest, BitWidth##bit_width##_TupleSize##tuple_size) { \
        runBasicInsertTest(bit_width, tuple_size); \
    }

// Generate tests
GENERATE_BASIC_TEST(1, 1)
GENERATE_BASIC_TEST(1, 3)
GENERATE_BASIC_TEST(1, 10)
GENERATE_BASIC_TEST(1, 50)
GENERATE_BASIC_TEST(1, 100)

GENERATE_BASIC_TEST(3, 3)
GENERATE_BASIC_TEST(3, 10)
GENERATE_BASIC_TEST(3, 100)

GENERATE_BASIC_TEST(7, 3)
GENERATE_BASIC_TEST(7, 10)
GENERATE_BASIC_TEST(7, 100)

GENERATE_BASIC_TEST(15, 3)
GENERATE_BASIC_TEST(15, 10)
GENERATE_BASIC_TEST(15, 100)

GENERATE_BASIC_TEST(31, 3)
GENERATE_BASIC_TEST(31, 10)
GENERATE_BASIC_TEST(31, 100)

GENERATE_BASIC_TEST(63, 3)
GENERATE_BASIC_TEST(63, 10)
GENERATE_BASIC_TEST(63, 100)

GENERATE_BASIC_TEST(64, 3)
GENERATE_BASIC_TEST(64, 10)
GENERATE_BASIC_TEST(64, 100)

// --- Stress Tests and Benchmarks ---

TEST_F(XlistTest, StressTest_BitWidth3_TupleSize3) {
    dist::cooc::Xlist xlist(3, 3);
    benchmarkInsert(xlist, 1000000, 3, 3);
    auto data = ExportData(xlist);
    EXPECT_GT(data.size(), 0);
}

TEST_F(XlistTest, StressTest_BitWidth15_TupleSize10) {
    dist::cooc::Xlist xlist(15, 10);
    benchmarkInsert(xlist, 500000, 10, 15);
    auto data = ExportData(xlist);
    EXPECT_GT(data.size(), 0);
}

TEST_F(XlistTest, StressTest_BitWidth31_TupleSize100) {
    dist::cooc::Xlist xlist(31, 100);
    benchmarkInsert(xlist, 100000, 100, 31);
    auto data = ExportData(xlist);
    EXPECT_GT(data.size(), 0);
}

// RankMap tests

class RankMapTest : public ::testing::Test {
protected:
    std::mt19937 rng{42};

    dist::ColourVector GenerateTuple(int size, int bit_width) {
        std::uniform_int_distribution<uint64_t> dist(0, (bit_width == 64) ? UINT64_MAX : ((1LL << bit_width) - 1));
        dist::ColourVector tuple(size);
        for (uint64_t& val : tuple)
            val = dist(rng);
        return tuple;
    }

    std::vector<std::pair<dist::ColourVector, uint32_t>> ExportData(dist::cooc::RankMap& map) const {
        std::vector<std::pair<dist::ColourVector, uint32_t>> result;
        map.forEach([&](auto begin, auto end, uint32_t count)->void
        {
            result.emplace_back(std::make_pair(dist::ColourVector(begin, end), count));
        });
        return result;
    }

    void InsertTuples(dist::cooc::RankMap& map,
                      std::map<dist::ColourVector, uint32_t>& expected,
                      std::mutex &mtx,
                      const std::vector<size_t>& sizes,
                      size_t bit_width,
                      size_t count_per_size)
    {
        std::uniform_int_distribution<uint64_t> dist(0, (bit_width == 64) ? UINT64_MAX : ((1LL << bit_width) - 1));

        for (size_t size : sizes) {
            dist::ColourVector tuple(size);
            for (size_t i = 0; i < count_per_size; ++i) {
                for (size_t& val : tuple)
                    val = dist(rng);

                map.increment(tuple);
                
                {
                    std::lock_guard<std::mutex> lock(mtx);
                    expected[tuple]++;
                }
            }
        }
    }
};

TEST_F(RankMapTest, MixedTupleSizesSingleThread_Correctness) {
    unsigned int bit_width = 15;
    dist::cooc::RankMap map(bit_width);

    std::map<dist::ColourVector, uint32_t> expected;
    std::mutex mtx;
    std::vector<size_t> sizes = {3, 10, 50, 100};

    InsertTuples(map, expected, mtx, sizes, bit_width, 5);

    auto data = ExportData(map);
    ASSERT_EQ(data.size(), expected.size());

    for (const auto& [tuple, count] : data) {
        auto it = expected.find(tuple);
        ASSERT_NE(it, expected.end()) << "Unexpected tuple found: size=" << tuple.size();
        EXPECT_EQ(count, it->second) << "Count mismatch for tuple of size " << tuple.size();
        EXPECT_EQ(tuple, it->first) << "Tuple content mismatch.";
    }
}

TEST_F(RankMapTest, ConcurrentInsertions_Correctness) {
    unsigned int bit_width = 31;
    dist::cooc::RankMap map(bit_width);

    std::map<dist::ColourVector, uint32_t> expected;
    std::mutex mtx;
    std::vector<size_t> sizes = {3, 10, 50};

    const size_t threads = 8;
    const size_t count_per_size = 1000;
    std::vector<std::thread> workers;

    for (size_t i = 0; i < threads; ++i) {
        workers.emplace_back([&] {
            InsertTuples(map, expected, mtx, sizes, bit_width, count_per_size);
        });
    }

    for (auto& t : workers) t.join();

    auto data = ExportData(map);
    ASSERT_EQ(data.size(), expected.size());

    for (const auto& [tuple, count] : data) {
        auto it = expected.find(tuple);
        ASSERT_NE(it, expected.end()) << "Unexpected tuple found.";
        EXPECT_EQ(count, it->second) << "Count mismatch for tuple.";
        EXPECT_EQ(tuple, it->first) << "Tuple content mismatch.";
    }
}
