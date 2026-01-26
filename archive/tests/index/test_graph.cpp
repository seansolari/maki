#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <map>
#include <sstream>
#include <vector>
#include <glog/logging.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <highfive/highfive.hpp>
#include <sdsl/vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>
#include <cereal/archives/binary.hpp>
#include <maki/maki.h>
#include <maki/fasta.hpp>
#include <maki/utils.hpp>
#include <maki/mers.hpp>
#include <maki/terminals.hpp>
#include <maki/graph.hpp>
#include <maki/flush.hpp>

#define USE_LOGS 0

namespace fs = std::filesystem;

// namespace declarations

using ::testing::ContainerEq;
using ::testing::ElementsAreArray;

class StringSerialiseTests : public testing::Test {
protected:
    StringSerialiseTests()
        : vec({ 1, 2, 3, 2, 2, 2, 2, 1, 0 })
        , tmpDir(tempio::create_temporary_directory())
    {}

    ~StringSerialiseTests()
    {
        LOG(INFO) << "removing temporary directory " << tmpDir;
        fs::remove_all(tmpDir);
    }

    sdsl::int_vector<0> vec;
    fs::path tmpDir;
};

TEST_F(StringSerialiseTests, IntVectorToFromString) {
    std::stringstream ss;
    {
        cereal::BinaryOutputArchive oarchive(ss);
        oarchive( vec );
    }

    sdsl::int_vector<0> out;
    {
        cereal::BinaryInputArchive iarchive(ss);
        iarchive( out );
    }

    ASSERT_EQ(vec, out);
}

TEST_F(StringSerialiseTests, EnvVectorToFromString) {
    sdsl::enc_vector<> enc(vec);

    std::stringstream ss;
    {
        cereal::BinaryOutputArchive oarchive(ss);
        oarchive( enc );
    }

    sdsl::enc_vector<> out;
    {
        cereal::BinaryInputArchive iarchive(ss);
        iarchive( out );
    }

    ASSERT_EQ(enc, out);
}

TEST_F(StringSerialiseTests, EnvVectorToFromMovedString) {
    sdsl::enc_vector<> enc(vec);

    std::ostringstream oss;
    {
        cereal::BinaryOutputArchive oarchive(oss);
        oarchive( enc );
    }

    std::string _data = std::move(oss).str();
    std::istringstream iss(std::move(_data));

    sdsl::enc_vector<> out;
    {
        cereal::BinaryInputArchive iarchive(iss);
        iarchive( out );
    }

    ASSERT_EQ(enc, out);
}

TEST_F(StringSerialiseTests, SaveLoadString) {
    // string test data

    sdsl::enc_vector<> enc(vec);

    std::ostringstream oss;
    {
        cereal::BinaryOutputArchive oarchive(oss);
        oarchive( enc );
    }

    std::string _data = std::move(oss).str();

    // save to disk

    {
        std::unique_ptr<HighFive::File> h5p = std::make_unique<HighFive::File>(
            tmpDir / "test.h5",
            HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate
            );
    
        // allocate datasets

        std::array<size_t, 1> dims = { _data.size() };
        HighFive::DataSet dset = h5p->createDataSet<std::remove_pointer_t<decltype(_data.data())>>("string_data", HighFive::DataSpace(dims));
        dset.write_raw(_data.data());
    }

    // load from disk

    std::string _in;
    _in.resize(_data.size());

    {
        std::unique_ptr<HighFive::File> h5p = std::make_unique<HighFive::File>(
            tmpDir / "test.h5",
            HighFive::File::ReadOnly
            );
    
        auto dset = h5p->getDataSet("string_data");
        dset.read_raw(_in.data());
    }

    ASSERT_EQ(_data, _in);
}

class ZstdTests : public testing::Test {
protected:
    ZstdTests()
        : buffer1(3, 10)
        , buffer2(3, 12)
        , buffer3(3, 15)
        , zcx()
        , zdx()
    {
        buffer1.colours.assign(   { 0, 1, 2, 3, 4, 5, 6, 7, 0, 1 });
        buffer1.boundaries.assign({ 1, 0, 0, 0, 0, 0, 0, 0, 1, 0 });

        buffer2.colours.assign(   { 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5 });
        buffer2.boundaries.assign({ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 });

        buffer3.colours.assign(   { 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4 });
        buffer3.boundaries.assign({ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 });
    }

    graph::ColourBuffer
        buffer1,
        buffer2,
        buffer3;
    ZstdCompressor zcx;
    ZstdDecompressor zdx;
};

TEST_F(ZstdTests, CompressDecompressOne) {
    graph::ColourBuffer t1(buffer1);

    // result buffers

    sdsl::int_vector<0> t1c;
    sdsl::bit_vector t1b;

    // compress-decompress

    zcx.compress(std::move(t1.colours));
    zdx.src = zcx.rsl;    // copy compressed data to decompressor
    zdx.decompress(t1c);

    zcx.compress(std::move(t1.boundaries));
    zdx.src = zcx.rsl;    // copy compressed data to decompressor
    zdx.decompress(t1b);

    // check

    ASSERT_EQ(t1c, buffer1.colours);
    ASSERT_EQ(t1b, buffer1.boundaries);
}

TEST_F(ZstdTests, CompressDecompressMany) {
    graph::ColourBuffer t1(buffer1),
                        t2(buffer2),
                        t3(buffer3);

    // result buffers

    sdsl::int_vector<0> t1c,
                        t2c,
                        t3c;
    sdsl::bit_vector
        t1b,
        t2b,
        t3b;

    // compress-decompress

    zcx.compress(std::move(t1.colours));
    zdx.src = zcx.rsl;    // copy compressed data to decompressor
    zdx.decompress(t1c);

    zcx.compress(std::move(t1.boundaries));
    zdx.src = zcx.rsl;    // copy compressed data to decompressor
    zdx.decompress(t1b);

    zcx.compress(std::move(t2.colours));
    zdx.src = zcx.rsl;    // copy compressed data to decompressor
    zdx.decompress(t2c);

    zcx.compress(std::move(t2.boundaries));
    zdx.src = zcx.rsl;    // copy compressed data to decompressor
    zdx.decompress(t2b);

    zcx.compress(std::move(t3.colours));
    zdx.src = zcx.rsl;    // copy compressed data to decompressor
    zdx.decompress(t3c);

    zcx.compress(std::move(t3.boundaries));
    zdx.src = zcx.rsl;    // copy compressed data to decompressor
    zdx.decompress(t3b);

    // check

    ASSERT_EQ(t1c, buffer1.colours);
    ASSERT_EQ(t1b, buffer1.boundaries);
    ASSERT_EQ(t2c, buffer2.colours);
    ASSERT_EQ(t2b, buffer2.boundaries);
    ASSERT_EQ(t3c, buffer3.colours);
    ASSERT_EQ(t3b, buffer3.boundaries);
}

static void CheckBufferEdgePositions(const graph::ColourBufferRegistry &buffers) {
    ZstdDecompressor zstd;
    // order of buffers
    std::vector<size_t> sortedKeys;
    for (auto const &rec : buffers) sortedKeys.push_back(rec.first);
    std::sort(sortedKeys.begin(), sortedKeys.end());
    // trace edges
    size_t edgeCount = 0;
    for (size_t const &key : sortedKeys) {
        graph::ColourBuffer src = buffers.readColours(buffers.getKey(key)->second, zstd);
        for (auto b : src.boundaries) edgeCount += b;
        ASSERT_EQ(edgeCount, key);
    }
}

class ColourBufferRegistryTests : public testing::Test {
protected:
    ColourBufferRegistryTests()
        : buffer1(3, 10)
        , buffer2(3, 12)
        , buffer3(3, 15)
        , outputPath(fs::temp_directory_path() / "__test_ColourBufferRegistryTests.h5")
        , colRegister(outputPath)
        , zcx()
        , zdx()
    {
        // create data

        buffer1.colours.assign(   { 0, 1, 2, 3, 4, 5, 6, 7, 0, 1 });
        buffer1.boundaries.assign({ 1, 0, 0, 0, 0, 0, 0, 0, 1, 0 });

        buffer2.colours.assign(   { 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5 });
        buffer2.boundaries.assign({ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 });

        buffer3.colours.assign(   { 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4 });
        buffer3.boundaries.assign({ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 });

        // save data

        colRegister.saveColours(0, std::move(buffer1), zcx);
        colRegister.saveColours(10, std::move(buffer2), zcx);
        colRegister.saveColours(22, std::move(buffer3), zcx);
    }

    ~ColourBufferRegistryTests()
    {
        fs::remove(outputPath);
    }

    graph::ColourBuffer
        buffer1,
        buffer2,
        buffer3;
    std::string outputPath;
    graph::ColourBufferRegistry colRegister;
    ZstdCompressor zcx;
    ZstdDecompressor zdx;
};

TEST_F(ColourBufferRegistryTests, ReadBuffer1) {
    const graph::BufferIdPair &key = colRegister.getKey(0)->second;
    graph::ColourBuffer res = colRegister.readColours(key, zdx);

    EXPECT_EQ(res.colours, buffer1.colours);
    EXPECT_EQ(res.boundaries, buffer1.boundaries);
}

TEST_F(ColourBufferRegistryTests, ReadBuffer2) {
    const graph::BufferIdPair &key = colRegister.getKey(10)->second;
    graph::ColourBuffer res = colRegister.readColours(key, zdx);

    EXPECT_EQ(res.colours, buffer2.colours);
    EXPECT_EQ(res.boundaries, buffer2.boundaries);
}

TEST_F(ColourBufferRegistryTests, ReadBuffer3) {
    const graph::BufferIdPair &key = colRegister.getKey(22)->second;
    graph::ColourBuffer res = colRegister.readColours(key, zdx);

    EXPECT_EQ(res.colours, buffer3.colours);
    EXPECT_EQ(res.boundaries, buffer3.boundaries);
}

TEST_F(ColourBufferRegistryTests, TestMaxColourWidth) {
    EXPECT_EQ(colRegister.getMaxColourWidth(), 3);
}

TEST_F(ColourBufferRegistryTests, TestTotalSize) {
    EXPECT_EQ(colRegister.getTotalSize(), 37);
}

TEST_F(ColourBufferRegistryTests, RecoverBuffers) {
    auto result = colRegister.combineAllBuffers();

    graph::ColourBuffer expected(3, 37);
    expected.colours.assign(   { 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4 });
    expected.boundaries.assign({ 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 });

    ASSERT_THAT(result.colours, ContainerEq(expected.colours));
    ASSERT_THAT(result.boundaries, ContainerEq(expected.boundaries));
}

//
// Serial inserting
// --------------------------------------------------

class GraphInsertTests : public testing::Test {
protected:
    GraphInsertTests()
        : k(7)
        , num_colours(11)
    {}

    uint8_t
        k,
        num_colours;
};

TEST_F(GraphInsertTests, WriteUncolouredEdges) {
    graph::BaseGraph graph(k, 1, 4);

    auto crs = graph.from(0);
    crs.emplace(BufferValue(0ull));
    crs.emplace(BufferValue(2ull));
    crs.emplace(BufferValue(3ull));
    crs.emplace(BufferValue(mers::KmerBuffer::terminalEdge));
    crs.emplace(BufferValue(0ull));
    crs.emplace(BufferValue(3ull));
    crs.emplace(BufferValue(1ull));
    crs.emplace(BufferValue(mers::KmerBuffer::terminalEdge));
    crs.emplace(BufferValue(1ull));
    crs.emplace(BufferValue(3ull));
    crs.writeNode(1u /* A */);

    // check graph structure

    sdsl::int_vector<4> edges = { 0b1001, 0b1010, 0b1011, 0b1100 };
    sdsl::bit_vector    wplus = { 0, 0, 0, 1 };
    
    EXPECT_EQ(crs.currentBlock, 0);
    EXPECT_EQ(crs.position(), 4);
    
    ASSERT_THAT(graph.edges_ref(), ContainerEq(edges));
    ASSERT_THAT(graph.wplus_ref(), ContainerEq(wplus));
}

TEST_F(GraphInsertTests, WriteColouredEdges) {
    graph::WriteableGraph graph(k, num_colours, num_colours, 4);
    
    auto crs = graph.from(0, 8, nullptr);

    crs.emplace(BufferValue((size_t)2ull, 0ull));
    crs.emplace(BufferValue((size_t)5ull, 2ull));
    crs.emplace(BufferValue((size_t)8ull, 3ull));
    crs.emplace(BufferValue((size_t)0ull, mers::KmerBuffer::terminalEdge));     // should be ignored
    crs.emplace(BufferValue((size_t)1ull, 0ull));
    crs.emplace(BufferValue((size_t)7ull, 3ull));
    crs.emplace(BufferValue((size_t)4ull, 1ull));
    crs.emplace(BufferValue((size_t)0ull, mers::KmerBuffer::terminalEdge));    // should be ignored
    crs.emplace(BufferValue((size_t)3ull, 1ull));
    crs.emplace(BufferValue((size_t)6ull, 3ull));

    crs.writeNode(1u /* A */);

    // check graph structure

    EXPECT_EQ(crs.currentBlock, 0);
    EXPECT_EQ(crs.position(), 4);
    EXPECT_EQ(crs.getColourPosition(), 8);

    sdsl::int_vector<4> edges = { 0b1001, 0b1010, 0b1011, 0b1100 };
    sdsl::bit_vector wplus = { 0, 0, 0, 1 };

    ASSERT_THAT(graph.edges_ref(), ContainerEq(edges));
    ASSERT_THAT(graph.wplus_ref(), ContainerEq(wplus));

    auto expColours = graph::ColourBuffer(graph.getColourBitWidth(), 8);
    expColours.colours.assign(   { 1, 2, 3, 4, 5, 6, 7, 8 });
    expColours.boundaries.assign({ 1, 0, 1, 0, 1, 1, 0, 0 });

    ASSERT_THAT(crs.getColours(), ContainerEq(expColours.colours));
    ASSERT_THAT(crs.getBoundaries(), ContainerEq(expColours.boundaries));
    crs.resetColourBuffer();
}

TEST_F(GraphInsertTests, InsertUncolouredKmers) {
    // create k-mers to insert
    mers::UncolouredKmerBuffer buffer(
        k,
        k,
        {
            { 0b01111011, 0b00001010, 0b00000001 }, // TGTCGGA -> C
            { 0b01111011, 0b00001010, 0b00000001 }, // TGTCGGA -> C
            { 0b10010000, 0b00011000, 0b00000011 }, // AACGAGC -> T
            { 0b10010000, 0b00011000, 0b00000011 }, // AACGAGC -> T
            { 0b10010001, 0b00011000, 0b00000000 }, // CACGAGC -> A
            { 0b10010001, 0b00011000, 0b00000000 }, // CACGAGC -> A
            { 0b10010011, 0b00011000, 0b00000000 }, // TACGAGC -> A
            { 0b10010011, 0b00011000, 0b00000000 }, // TACGAGC -> A
            { 0b01110000, 0b00111010, 0b00000001 }, // AATCGGT -> C
            { 0b01110011, 0b00111010, 0b00000001 }, // TATCGGT -> C
        }
    );
    mers::KmerOverlapVector B = {
        mers::BW_0_K, mers::IS_0, mers::BW_0_K, mers::IS_0,
        mers::IS_K, mers::IS_0, mers::IS_K, mers::IS_0,
        mers::BW_0_K, mers::IS_K
    };

    graph::BaseGraph graph(k, 1, 6);
    auto crs = graph.from(0);

    flush::insert(crs, buffer.begin(), buffer.end(), B.begin());

    // check graph structure

    sdsl::int_vector<4> edges = { 0b1010, 0b1100, 0b1001, 0b0001, 0b1010, 0b0010 };
    sdsl::bit_vector wplus = { 1, 1, 1, 1, 1, 1 };

    ASSERT_THAT(graph.edges_ref(), ContainerEq(edges));
    ASSERT_THAT(graph.wplus_ref(), ContainerEq(wplus));

    ASSERT_THAT(crs.C, ElementsAreArray({ 0, 1, 3, 0, 2 }));
    ASSERT_THAT(crs.F, ElementsAreArray({ 0, 1, 3, 0, 2 }));
    ASSERT_THAT(crs.position(), 6);
}

TEST_F(GraphInsertTests, InsertColouredKmers) {
    // create k-mers to insert
    mers::KmerBuffer buffer(
        _mkimem_details::value_size(ceil_log2(num_colours)),
        k,
        k,
        {
            { 0b01111011, 0b00001010, 0b00000001 }, // TGTCGGA -> C
            { 0b01111011, 0b00001010, 0b00001001 }, // TGTCGGA -> C
            { 0b10010000, 0b00011000, 0b00000011 }, // AACGAGC -> T
            { 0b10010000, 0b00011000, 0b00011011 }, // AACGAGC -> T
            { 0b10010001, 0b00011000, 0b00000000 }, // CACGAGC -> A
            { 0b10010001, 0b00011000, 0b00111000 }, // CACGAGC -> A
            { 0b10010011, 0b00011000, 0b00001000 }, // TACGAGC -> A
            { 0b10010011, 0b00011000, 0b00101000 }, // TACGAGC -> A
            { 0b01110000, 0b00111010, 0b00000001 }, // AATCGGT -> C
            { 0b01110011, 0b00111010, 0b00000001 }, // TATCGGT -> C
        }
    );
    mers::KmerOverlapVector B = {
        mers::BW_0_K, mers::IS_0, mers::BW_0_K, mers::IS_0,
        mers::IS_K, mers::IS_0, mers::IS_K, mers::IS_0,
        mers::BW_0_K, mers::IS_K
    };

    graph::WriteableGraph graph(k, num_colours, num_colours, 6);
    auto crs = graph.from(0, 10, nullptr);

    flush::insert(crs, buffer.begin(), buffer.end(), B.begin());

    // check graph structure

    sdsl::int_vector<4> edges = { 0b1010, 0b1100, 0b1001, 0b0001, 0b1010, 0b0010 };
    sdsl::bit_vector wplus = { 1, 1, 1, 1, 1, 1 };

    ASSERT_THAT(graph.edges_ref(), ContainerEq(edges));
    ASSERT_THAT(graph.wplus_ref(), ContainerEq(wplus));

    graph::ColourBuffer expColours(graph.getColourBitWidth(), 10);
    expColours.colours.assign(   { 0, 1, 0, 3, 0, 7, 1, 5, 0, 0 });
    expColours.boundaries.assign({ 1, 0, 1, 0, 1, 0, 1, 0, 1, 1 });

    ASSERT_THAT(crs.getColours(), ContainerEq(expColours.colours));
    ASSERT_THAT(crs.getBoundaries(), ContainerEq(expColours.boundaries));

    ASSERT_THAT(crs.C, ElementsAreArray({ 0, 1, 3, 0, 2 }));
    ASSERT_THAT(crs.F, ElementsAreArray({ 0, 1, 3, 0, 2 }));
    ASSERT_THAT(crs.position(), 6);
    ASSERT_THAT(crs.getColourPosition(), 10);
    crs.resetColourBuffer();
}

TEST_F(GraphInsertTests, InsertKmers_Suffix) {
    mers::KmerBuffer buffer(
        _mkimem_details::value_size(ceil_log2(num_colours)),
        k,
        k,
        {
            { 0b01111011, 0b00001010, 0b00000001 }, // TGTCGGA -> C
            { 0b01111011, 0b00001010, 0b00001001 }, // TGTCGGA -> C
            { 0b10010000, 0b00011000, 0b00000011 }, // AACGAGC -> T
            { 0b10010000, 0b00011000, 0b00011011 }, // AACGAGC -> T
            { 0b10010001, 0b00011000, 0b00000000 }, // CACGAGC -> A
            { 0b10010001, 0b00011000, 0b00111000 }, // CACGAGC -> A
            { 0b10010011, 0b00011000, 0b00001000 }, // TACGAGC -> A
            { 0b10010011, 0b00011000, 0b00101000 }, // TACGAGC -> A
            { 0b01110000, 0b00111010, 0b00000001 }, // AATCGGT -> C
            { 0b01110011, 0b00111010, 0b00000001 }, // TATCGGT -> C
        }
    );
    mers::KmerOverlapVector B = {
        mers::BW_0_K, mers::IS_0, mers::BW_0_K, mers::IS_0,
        mers::IS_K, mers::IS_0, mers::IS_K, mers::IS_0,
        mers::BW_0_K, mers::IS_K
    };

    graph::WriteableGraph graph(k, num_colours, num_colours, 6);
    auto crs = graph.from(0, 10, nullptr);

    flush::insert(crs, buffer.begin(), buffer.end(), B.begin(), (uint8_t)0b01u);

    // check graph structure

    sdsl::int_vector<4> edges = { 0b1010, 0b1100, 0b1001, 0b0001, 0b1010, 0b0010 };
    sdsl::bit_vector wplus = { 1, 1, 1, 1, 1, 1 };

    ASSERT_THAT(graph.edges_ref(), ContainerEq(edges));
    ASSERT_THAT(graph.wplus_ref(), ContainerEq(wplus));

    graph::ColourBuffer expColours(graph.getColourBitWidth(), 10);
    expColours.colours.assign(   { 0, 1, 0, 3, 0, 7, 1, 5, 0, 0 });
    expColours.boundaries.assign({ 1, 0, 1, 0, 1, 0, 1, 0, 1, 1 });

    ASSERT_THAT(crs.getColours(), ContainerEq(expColours.colours));
    ASSERT_THAT(crs.getBoundaries(), ContainerEq(expColours.boundaries));

    ASSERT_THAT(crs.C, ElementsAreArray({ 0, 6, 0, 0, 0 }));
    ASSERT_THAT(crs.F, ElementsAreArray({ 0, 6, 0, 0, 0 }));
    ASSERT_THAT(crs.position(), 6);
    ASSERT_THAT(crs.getColourPosition(), 10);
    crs.resetColourBuffer();
}

/**
 * $$$$ -$$$ -> A
 * $$$$ -$$$ -> C
 * $$$$ -AAA -> G
 * AAA$ -AAA -> G
 * AAC$ -AAA -> T
 * $$$$ -TAA -> A
 * $$$$ -TAA -> C
 * 
 */
TEST_F(GraphInsertTests, InsertTerminals) {
    terminals::TerminalBuffer terminals(
        k,
        {
            { /* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000000 },
            { /* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000001 },
            { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000010 },
            { /* size */ 0b00000110, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000010 },
            { /* size */ 0b00000110, /* kmer */ 0b00000100, 0b00000000, /* edge */ 0b00000011 },
            { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00110000, /* edge */ 0b00000000 },
            { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00110000, /* edge */ 0b00000001 }
        },
        terminals::TerminalBuffer::autofit_tag);
    mers::KmerOverlapVector B = {
        mers::BW_0_K, mers::IS_0, mers::BW_0_K,
        mers::BW_0_K, mers::BW_0_K, mers::BW_0_K,
        mers::IS_0
    };

    graph::WriteableGraph graph(k, num_colours, num_colours, 7);
    auto crs = graph.from(0, 7, nullptr);

    flush::insert(crs, terminals.constBegin(), terminals.constEnd(), B.begin());

    // check graph structure

    sdsl::int_vector<4> edges = { 0b1001, 0b1010, 0b1011, 0b1011, 0b1100, 0b1001, 0b1010 };
    sdsl::bit_vector wplus = { 0, 1, 1, 1, 1, 0, 1 };

    ASSERT_THAT(graph.edges_ref(), ContainerEq(edges));
    ASSERT_THAT(graph.wplus_ref(), ContainerEq(wplus));

    auto expColours = graph::ColourBuffer(graph.getColourBitWidth(), 10);
    expColours.colours.assign(   { 0, 0, 0, 0, 0, 0, 0 });
    expColours.boundaries.assign({ 1, 1, 1, 1, 1, 1, 1 });

    ASSERT_THAT(crs.getColours(), ContainerEq(expColours.colours));
    ASSERT_THAT(crs.getBoundaries(), ContainerEq(expColours.boundaries));

    ASSERT_THAT(crs.C, ElementsAreArray({ 1, 3, 0, 0, 1 }));
    ASSERT_THAT(crs.F, ElementsAreArray({ 2, 3, 0, 0, 2 }));
    ASSERT_THAT(crs.position(), 7);
    ASSERT_THAT(crs.getColourPosition(), 7);
    crs.resetColourBuffer();
}

/**
 * $$$$ -$$$ -> A
 * $$$$ -$$$ -> C
 * $$$$ -AAA -> G
 * AAA$ -AAA -> G
 * AAC$ -AAA -> T
 * AACT -AAA -> T
 * $$$$ -TAA -> A
 * $$$$ -TAA -> C
 * 
 */
TEST_F(GraphInsertTests, InsertTerminalKmers) {
    terminals::TerminalBuffer terminals(
        k,
        {
            { /* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000000 },
            { /* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000001 },
            { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000010 },
            { /* size */ 0b00000110, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000010 },
            { /* size */ 0b00000110, /* kmer */ 0b00000100, 0b00000000, /* edge */ 0b00000011 },
            { /* size */ 0b00000111, /* kmer */ 0b00000111, 0b00000000, /* edge */ 0b00000011 },
            { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00110000, /* edge */ 0b00000000 },
            { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00110000, /* edge */ 0b00000001 }
        },
        terminals::TerminalBuffer::autofit_tag);
    mers::KmerOverlapVector B = {
        mers::BW_0_K, mers::IS_0, mers::BW_0_K,
        mers::BW_0_K, mers::BW_0_K, mers::IS_K,
        mers::BW_0_K, mers::IS_0
    };

    graph::BaseGraph graph(k, 1, 8);
    auto crs = graph.from(0);

    flush::insert(crs, terminals.constBegin(), terminals.constEnd(), B.begin());

    // check graph structure

    sdsl::int_vector<4> edges = { 0b1001, 0b1010, 0b1011, 0b1011, 0b1100, 0b0100, 0b1001, 0b1010 };
    sdsl::bit_vector wplus = { 0, 1, 1, 1, 1, 1, 0, 1 };

    ASSERT_THAT(graph.edges_ref(), ContainerEq(edges));
    ASSERT_THAT(graph.wplus_ref(), ContainerEq(wplus));

    ASSERT_THAT(crs.C, ElementsAreArray({ 1, 4, 0, 0, 1 }));
    ASSERT_THAT(crs.F, ElementsAreArray({ 2, 4, 0, 0, 2 }));
    ASSERT_THAT(crs.position(), 8);
}

class FlushRangeTests : public testing::Test {
protected:
    FlushRangeTests()
        : k(7)
        , terminals(
            k,
            {
                { /* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000000 },
                { /* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000000 },
                { /* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000000 },
                { /* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000001 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000001 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000010 },
                { /* size */ 0b00000110, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000010 },
                { /* size */ 0b00000110, /* kmer */ 0b00000100, 0b00000000, /* edge */ 0b00000011 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00110000, /* edge */ 0b00000000 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00110000, /* edge */ 0b00000001 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00110000, /* edge */ 0b00000000 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00110000, /* edge */ 0b00000001 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00110000, /* edge */ 0b00000011 }
            },
            terminals::TerminalBuffer::autofit_tag)
        , Bt(adjacentDifference(terminals))
    {
    }

    uint8_t k;
    terminals::TerminalBuffer terminals;
    mers::KmerOverlapVector Bt;
};

TEST_F(FlushRangeTests, ForwardRangeInit) {
    flush::ranges::detail::ForwardRange Rng(terminals.begin(), terminals.end());
    ASSERT_EQ(Rng.size(), 13);
    ASSERT_FALSE(Rng.empty());
}

TEST_F(FlushRangeTests, ForwardRangeMove) {
    flush::ranges::detail::ForwardRange Rng(terminals.begin(), terminals.end());
    Rng.stepForward(13);
    ASSERT_EQ(Rng.size(), 0);
    ASSERT_TRUE(Rng.empty());
}

TEST_F(FlushRangeTests, MonotoneRangeInit) {
    flush::ranges::detail::MonotoneRange Rng(terminals.begin(), terminals.end(), Bt.begin());
    ASSERT_EQ(Rng.size(), 13);
    ASSERT_EQ(Rng.B(), Bt.begin());
    ASSERT_FALSE(Rng.empty());
}

TEST_F(FlushRangeTests, MonotoneRangeMove) {
    flush::ranges::detail::MonotoneRange Rng(terminals.begin(), terminals.end(), Bt.begin());
    Rng.stepForward(13);
    ASSERT_EQ(Rng.size(), 0);
    ASSERT_EQ(Rng.B(), Bt.end());
    ASSERT_TRUE(Rng.empty());
}

class GranularSplitTests : public testing::Test {
protected:
    GranularSplitTests()
        : k(7)
        , num_records(20)
        , kmers(num_records, k, terminals::TerminalBuffer::autofit_tag)
        , B({
            mers::BW_0_K, mers::IS_0, mers::IS_0, mers::IS_0, mers::IS_0,
            mers::BW_0_K, mers::IS_0, mers::IS_0, mers::IS_0, mers::IS_0,
            mers::BW_0_K, mers::IS_0, mers::IS_0, mers::IS_0, mers::IS_0,
            mers::BW_0_K, mers::IS_0, mers::IS_0, mers::IS_0, mers::IS_0
        })
        , Bblocks({
            1, 0, 0, 0, 0,
            1, 0, 0, 0, 0,
            1, 0, 0, 0, 0,
            1, 0, 0, 0, 0
        })
    {
        oneapi::tbb::parallel_invoke(
            [&]{ sdsl::util::init_support(rank,   &Bblocks); },
            [&]{ sdsl::util::init_support(select, &Bblocks); }
        );
    }

    uint8_t k;
    size_t num_records;
    terminals::TerminalBuffer kmers;
    sdsl::int_vector<2> B;
    sdsl::bit_vector Bblocks;
    sdsl::rank_support_v5<1,1> rank;
    sdsl::select_support_mcl<1,1> select;
};

TEST_F(GranularSplitTests, GranularCalcEvenChunk) {
    ASSERT_EQ(
        flush::ranges::detail::calculateLeftChunks(10, 5, { 1, 1 }),
        1);
}

TEST_F(GranularSplitTests, GranularCalcUnbalancedChunk) {
    ASSERT_EQ(
        flush::ranges::detail::calculateLeftChunks(100, 30, { 3, 1 }),
        3);
}

TEST_F(GranularSplitTests, GranularCalcUnbalancedSmallChunk) {
    ASSERT_EQ(
        flush::ranges::detail::calculateLeftChunks(10, 1, { 3, 2 }),
        6);
}

TEST_F(GranularSplitTests, TotalChunks) {
    flush::ranges::detail::Granular rng(1, rank, select);
    ASSERT_EQ(rng.chunks(0,  20), 4);
    ASSERT_EQ(rng.chunks(0,  5),  1);
    ASSERT_EQ(rng.chunks(5,  10), 1);
    ASSERT_EQ(rng.chunks(10, 15), 1);
    ASSERT_EQ(rng.chunks(15, 20), 1);
}

TEST_F(GranularSplitTests, GranularChunkPositions) {
    flush::ranges::detail::Granular rng(1, rank, select);
    ASSERT_EQ(rng.findChunkPosition(0, 20, 1), 5);
    ASSERT_EQ(rng.findChunkPosition(0, 20, 2), 10);
    ASSERT_EQ(rng.findChunkPosition(0, 20, 3), 15);
    ASSERT_EQ(rng.findChunkPosition(0, 20, 4), 20);
    ASSERT_EQ(rng.findChunkPosition(0, 20, 10), 20);
}

TEST_F(GranularSplitTests, PivotSplit) {
    flush::ranges::detail::DivisibleRange R1(kmers.begin(), kmers.end(), B.begin(), 1, rank, select);
    flush::ranges::detail::DivisibleRange R2(R1, 10, flush::ranges::detail::pivot_tag);
    ASSERT_EQ(R1.size(), 10);
    ASSERT_FALSE(R1.empty());
    ASSERT_TRUE(R1.is_divisible());
    ASSERT_EQ(R2.size(), 10);
    ASSERT_FALSE(R2.empty());
    ASSERT_TRUE(R2.is_divisible());
}

TEST_F(GranularSplitTests, ChunkSplit) {
    flush::ranges::detail::DivisibleRange R1(kmers.begin(), kmers.end(), B.begin(), 1, rank, select);
    flush::ranges::detail::DivisibleRange R2(R1, 1, flush::ranges::detail::chunk_tag);
    ASSERT_EQ(R1.size(), 15);
    ASSERT_FALSE(R1.empty());
    ASSERT_TRUE(R1.is_divisible());
    ASSERT_EQ(R2.size(), 5);
    ASSERT_FALSE(R2.empty());
    ASSERT_TRUE(R2.is_divisible());
}

TEST_F(GranularSplitTests, PropSplit) {
    flush::ranges::detail::DivisibleRange R1(kmers.begin(), kmers.end(), B.begin(), 1, rank, select);
    flush::ranges::detail::DivisibleRange R2(R1, { 5, 4 });
    ASSERT_EQ(R1.size(), 10);
    ASSERT_FALSE(R1.empty());
    ASSERT_TRUE(R1.is_divisible());
    ASSERT_EQ(R2.size(), 10);
    ASSERT_FALSE(R2.empty());
    ASSERT_TRUE(R2.is_divisible());
}

TEST_F(GranularSplitTests, PropTakeAll) {
    flush::ranges::detail::DivisibleRange R1(kmers.begin(), kmers.end(), B.begin(), 1, rank, select);
    flush::ranges::detail::DivisibleRange R2(R1, { 10, 1 });
    ASSERT_EQ(R1.size(), 5);
    ASSERT_FALSE(R1.empty());
    ASSERT_TRUE(R1.is_divisible());
    ASSERT_EQ(R2.size(), 15);
    ASSERT_FALSE(R2.empty());
    ASSERT_TRUE(R2.is_divisible());
}

TEST_F(GranularSplitTests, EvenSplit) {
    flush::ranges::detail::DivisibleRange R1(kmers.begin(), kmers.end(), B.begin(), 1, rank, select);
    flush::ranges::detail::DivisibleRange R2(R1, oneapi::tbb::split{});
    ASSERT_EQ(R1.size(), 10);
    ASSERT_FALSE(R1.empty());
    ASSERT_TRUE(R1.is_divisible());
    ASSERT_EQ(R2.size(), 10);
    ASSERT_FALSE(R2.empty());
    ASSERT_TRUE(R2.is_divisible());
}

TEST_F(GranularSplitTests, TerminalPivot) {
    terminals::TerminalRange T1 = kmers.asRange();
    terminals::TerminalRange T2(T1, T1.constBegin() + 15);
    ASSERT_EQ(T2.size(), 15);
    ASSERT_EQ(T2.constBegin(), kmers.constBegin());
    ASSERT_EQ(T2.constEnd(), kmers.constBegin() + 15);
    ASSERT_EQ(T1.size(),  5);
    ASSERT_EQ(T1.constBegin(), T2.constEnd());
    ASSERT_EQ(T1.constEnd(), kmers.constEnd());
}

class UniqueCountTests : public testing::Test {
protected:
    UniqueCountTests()
        : k(7)
        , num_colours(11)
        , kmers(
            _mkimem_details::value_size(ceil_log2(num_colours)),
            k,
            k,
            {
                { /* kmer */ 0b01111011, 0b00001010, /* edge */ 0b00000001 }, // TGTCGGA -> C
                { /* kmer */ 0b01111011, 0b00001010, /* edge */ 0b00001001 }, // TGTCGGA -> C
                { /* kmer */ 0b01111011, 0b00001010, /* edge */ 0b00001011 }, // TGTCGGA -> T
                { /* kmer */ 0b10010000, 0b00011000, /* edge */ 0b00000011 }, // AACGAGC -> T
                { /* kmer */ 0b10010000, 0b00011000, /* edge */ 0b00011011 }, // AACGAGC -> T
                { /* kmer */ 0b10010001, 0b00011000, /* edge */ 0b00000000 }, // CACGAGC -> A
                { /* kmer */ 0b10010001, 0b00011000, /* edge */ 0b00111000 }, // CACGAGC -> A
                { /* kmer */ 0b10010011, 0b00011000, /* edge */ 0b00001000 }, // TACGAGC -> A
                { /* kmer */ 0b10010011, 0b00011000, /* edge */ 0b00101000 }, // TACGAGC -> A
                { /* kmer */ 0b01110000, 0b00111010, /* edge */ 0b00000001 }, // AATCGGT -> C
                { /* kmer */ 0b01110011, 0b00111010, /* edge */ 0b00000001 }, // TATCGGT -> C
            })
        , Bk(adjacentDifference(kmers))
        , terminals(
            k,
            {
                { /* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000000 },
                { /* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000000 },
                { /* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000000 },
                { /* size */ 0b00000000, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000001 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000001 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000010 },
                { /* size */ 0b00000110, /* kmer */ 0b00000000, 0b00000000, /* edge */ 0b00000010 },
                { /* size */ 0b00000110, /* kmer */ 0b00000100, 0b00000000, /* edge */ 0b00000011 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00110000, /* edge */ 0b00000000 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00110000, /* edge */ 0b00000001 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00110000, /* edge */ 0b00000000 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00110000, /* edge */ 0b00000001 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00110000, /* edge */ 0b00000011 }
            },
            terminals::TerminalBuffer::autofit_tag)
        , Bt(adjacentDifference(terminals))
    {}

    uint8_t
        k,
        num_colours;
    mers::KmerBuffer kmers;
    mers::KmerOverlapVector Bk;
    terminals::TerminalBuffer terminals;
    mers::KmerOverlapVector Bt;
};

TEST_F(UniqueCountTests, CountKmersStart) {
    auto uniqueKmers = flush::ranges::counting::uniqueEdges(kmers.begin(), kmers.begin() + 3, Bk.begin(), 5);
    EXPECT_EQ(uniqueKmers, 2);
}

TEST_F(UniqueCountTests, CountKmersFull) {
    auto uniqueKmers = flush::ranges::counting::uniqueEdges(kmers.begin(), kmers.end(), Bk.begin(), 5);
    EXPECT_EQ(uniqueKmers, 7);
}

TEST_F(UniqueCountTests, CountTerminalsStart) {
    auto uniqueTerminals = flush::ranges::counting::uniqueEdges(terminals.constBegin(), terminals.constBegin() + 5, Bt.begin(), 5);
    EXPECT_EQ(uniqueTerminals, 3);
}

TEST_F(UniqueCountTests, CountTerminalsEnd) {
    auto uniqueTerminals = flush::ranges::counting::uniqueEdges(terminals.constBegin() + 8, terminals.constEnd(), Bt.begin() + 8, 5);
    EXPECT_EQ(uniqueTerminals, 3);
}

TEST_F(UniqueCountTests, CountTerminalsFull) {
    auto uniqueTerminals = flush::ranges::counting::uniqueEdges(terminals.constBegin(), terminals.constEnd(), Bt.begin(), 5);
    EXPECT_EQ(uniqueTerminals, 9);
}

TEST_F(UniqueCountTests, CountTerminalsFullLargeChunk) {
    auto uniqueTerminals = flush::ranges::counting::uniqueEdges(terminals.constBegin(), terminals.constEnd(), Bt.begin(), 1000000);
    EXPECT_EQ(uniqueTerminals, 9);
}

/**
 * Example 1.
 * 
 * Sorted Structure
 * ----------------
 *                                                        // $$$$$$$$ -> A
 *                                                        // $$$$$$$$ -> G
 *                                                        // $$$$$$$$ -> T
 *                                                        // $$$$$$$A -> A
 *                                                        // $$$$GAAA -> A
 * { 0b01111011, 0b00001010, 0, 0, 0, 0, 0, 0b00000001 }, // TGTCGGAA -> C
 * { 0b01111011, 0b00001010, 0, 0, 0, 0, 0, 0b00001001 }, // TGTCGGAA -> C
 *                                                        // $$$$$GCA -> A
 *                                                        // $$CGAGCA -> T
 * { 0b10010000, 0b00011000, 0, 0, 0, 0, 0, 0b00011011 }, // AACGAGCA -> T
 * { 0b10010001, 0b00011000, 0, 0, 0, 0, 0, 0b00000000 }, // CACGAGCA -> A
 * { 0b10010001, 0b00011000, 0, 0, 0, 0, 0, 0b00111000 }, // CACGAGCA -> A
 * { 0b10010011, 0b00011000, 0, 0, 0, 0, 0, 0b00001000 }, // TACGAGCA -> A
 * { 0b10010011, 0b00011000, 0, 0, 0, 0, 0, 0b00101000 }, // TACGAGCA -> A
 *                                                        // $ATCGGTA -> C
 * { 0b01110011, 0b00111010, 0, 0, 0, 0, 0, 0b00000001 }, // TATCGGTA -> C
 *                                                        // $$$$$$$C -> T
 * 
 */
class GraphFlushingTests : public testing::Test {
protected:
    GraphFlushingTests()
        : k(10)
        , k_eff(8)
        , num_colours(7)
        , kmers(
            _mkimem_details::value_size(ceil_log2(num_colours)),
            k,
            k_eff,
            {
                { /* kmer */ 0b01111011, 0b00001010, /* edge */ 0b00000001 },
                { /* kmer */ 0b01111011, 0b00001010, /* edge */ 0b00001001 },
                { /* kmer */ 0b10010000, 0b00011000, /* edge */ 0b00011011 },
                { /* kmer */ 0b10010001, 0b00011000, /* edge */ 0b00000000 },
                { /* kmer */ 0b10010001, 0b00011000, /* edge */ 0b00111000 },
                { /* kmer */ 0b10010011, 0b00011000, /* edge */ 0b00001000 },
                { /* kmer */ 0b10010011, 0b00011000, /* edge */ 0b00101000 },
                { /* kmer */ 0b01110011, 0b00111010, /* edge */ 0b00000001 }
            })
        , Bk(indexedAdjacentDifference<BW_0_K>(kmers))
        , terminals(
            k,
            {
                { /* size */ 0b00000010, /* kmer */ 0b00000000, 0b00000000, 0b00001111, /* edge */ 0b00000000 },
                { /* size */ 0b00000010, /* kmer */ 0b00000000, 0b00000000, 0b00001111, /* edge */ 0b00000010 },
                { /* size */ 0b00000010, /* kmer */ 0b00000000, 0b00000000, 0b00001111, /* edge */ 0b00000011 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b00000000, 0b00001111, /* edge */ 0b00000000 },
                { /* size */ 0b00000110, /* kmer */ 0b00000000, 0b00000010, 0b00001111, /* edge */ 0b00000000 },
                { /* size */ 0b00000101, /* kmer */ 0b00000000, 0b00011000, 0b00001111, /* edge */ 0b00000000 },
                { /* size */ 0b00001000, /* kmer */ 0b10010000, 0b00011000, 0b00001111, /* edge */ 0b00000011 },
                { /* size */ 0b00001001, /* kmer */ 0b01110000, 0b00111010, 0b00001111, /* edge */ 0b00000001 },
                { /* size */ 0b00000011, /* kmer */ 0b00000000, 0b01000000, 0b00001111, /* edge */ 0b00000011 }
            }, terminals::TerminalBuffer::autofit_tag)
        , Bt(adjacentDifference(terminals))
        , LessThan(k_eff, terminals.lengthBytes())
    {
        oneapi::tbb::parallel_invoke(
            [&]{ sdsl::util::init_support(rank,   &std::get<1>(Bk)); },
            [&]{ sdsl::util::init_support(select, &std::get<1>(Bk)); }
        );
    }

    uint8_t
        k,
        k_eff;
    size_t num_colours;
    mers::KmerBuffer kmers;
    std::pair<sdsl::int_vector<2>,sdsl::bit_vector> Bk;
    terminals::TerminalBuffer terminals;
    mers::KmerOverlapVector Bt;
    sdsl::rank_support_v5<1,1> rank;
    sdsl::select_support_mcl<1,1> select;
    terminals::MaskedAlignedBytesLessThan<terminals::RHS> LessThan;

    inline auto makeRange() {
        return flush::ranges::detail::InterleaveRanges(
            kmers.begin(), kmers.end(), std::get<0>(Bk).begin(),
            1, rank, select,
            terminals.asRange(), Bt.begin(),
            kmers.getK(), kmers.getEffK(),
            0);
    }
};

/**
 * MaskedBytesNotEqual used to compared terminals (LHS, masked) to k-mers (RHS, unmasked).
 * 
 */
TEST_F(GraphFlushingTests, MaskedBytesDiff_4k_bw0k) {
    mers::MaskedBytesDiff Diff { 8 };
    std::vector<uint8_t>
        tmn = { /* size */ 4, /* kmer */ 0b00000000, 0b00010111, /* edge */ 0b00000000 },
        kmr = {               /* kmer */ 0b11000000, 0b00010111, /* edge */ 0b00001000 };
    EXPECT_EQ(Diff(tmn.data() + 1, kmr.data()), mers::BW_0_K);
}

TEST_F(GraphFlushingTests, MaskedBytesDiff_4k_isk) {
    mers::MaskedBytesDiff Diff { 8 };
    std::vector<uint8_t>
        tmn = { /* size */ 4, /* kmer */ 0b00000000, 0b00010111, /* edge */ 0b00000000 },
        kmr = {               /* kmer */ 0b00000001, 0b00010111, /* edge */ 0b00001000 };
    EXPECT_EQ(Diff(tmn.data() + 1, kmr.data()), mers::IS_K);
}

TEST_F(GraphFlushingTests, MaskedBytesDiff_4k_is0) {
    mers::MaskedBytesDiff Diff { 8 };
    std::vector<uint8_t>
        tmn = { /* size */ 7, /* kmer */ 0b11100000, 0b01101100, /* edge */ 0b00000000 },
        kmr = {               /* kmer */ 0b11100000, 0b01101100, /* edge */ 0b00001000 };
    EXPECT_EQ(Diff(tmn.data() + 1, kmr.data()), mers::IS_0);
}

TEST_F(GraphFlushingTests, MaskedBytesDiff_Truncated_bw0k) {
    mers::MaskedBytesDiff Diff { 7 };
    std::vector<uint8_t>
        tmn = { /* size */ 10, /* kmer */ 0b10010011, 0b11001001, 0b00001111, /* edge */ 0b00000000 },
        kmr = {                /* kmer */ 0b00010011, 0b00001001,             /* edge */ 0b00001000 };
    EXPECT_EQ(Diff(tmn.data() + 1, kmr.data()), mers::BW_0_K);
}

TEST_F(GraphFlushingTests, MaskedBytesDiff_Truncated_isk) {
    mers::MaskedBytesDiff Diff { 7 };
    std::vector<uint8_t>
        tmn = { /* size */ 10, /* kmer */ 0b10010011, 0b11001001, 0b00001111, /* edge */ 0b00000000 },
        kmr = {                /* kmer */ 0b10010010, 0b00001001,             /* edge */ 0b00001000 };
    EXPECT_EQ(Diff(tmn.data() + 1, kmr.data()), mers::IS_K);
}

TEST_F(GraphFlushingTests, MaskedBytesDiff_Truncated_is0) {
    mers::MaskedBytesDiff Diff { 7 };
    std::vector<uint8_t>
        tmn = { /* size */ 10, /* kmer */ 0b10010011, 0b11001001, 0b00001111, /* edge */ 0b00000000 },
        kmr = {                /* kmer */ 0b10010011, 0b00001001,             /* edge */ 0b00001000 };
    EXPECT_EQ(Diff(tmn.data() + 1, kmr.data()), mers::IS_0);
}

TEST_F(GraphFlushingTests, LowerBound_KmerQuery_Seen) {
    auto it = std::lower_bound(kmers.begin(), kmers.end(), *terminals.constAt(6), LessThan);
    ASSERT_EQ(it, kmers.at(2));
}

TEST_F(GraphFlushingTests, LowerBound_KmerQuery_Begin) {
    auto it = std::lower_bound(kmers.begin(), kmers.end(), *terminals.constBegin(), LessThan);
    ASSERT_EQ(it, kmers.begin());
}

TEST_F(GraphFlushingTests, LowerBound_KmerQuery_End) {
    auto it = std::lower_bound(kmers.begin(), kmers.end(), *terminals.constAt(8), LessThan);
    ASSERT_EQ(it, kmers.end());
}

TEST_F(GraphFlushingTests, UpperBound_TerminalQuery_First) {
    auto it = std::upper_bound(terminals.constBegin(), terminals.constEnd(), *kmers.begin(), LessThan);
    ASSERT_EQ(it, terminals.constAt(5));
}

TEST_F(GraphFlushingTests, UpperBound_TerminalQuery_Seen) {
    auto it = std::upper_bound(terminals.constBegin(), terminals.constEnd(), *kmers.at(2), LessThan);
    ASSERT_EQ(it, terminals.constAt(7));
}

TEST_F(GraphFlushingTests, UpperBound_TerminalQuery_Begin) {
    std::vector<uint8_t> qry = { 0b00000000, 0b00000000, 0, 0, 0, 0, 0, 0b00000001 };
    auto it = std::upper_bound(terminals.constBegin(), terminals.constEnd(), ndim::Span<uint8_t *>{ qry.data(), 8 }, LessThan);
    ASSERT_EQ(it, terminals.constAt(4));
}

TEST_F(GraphFlushingTests, UpperBound_TerminalQuery_End) {
    std::vector<uint8_t> qry = { 0b00000000, 0b01000000, 0, 0, 0, 0, 0, 0b00000001 };
    auto it = std::upper_bound(terminals.constBegin(), terminals.constEnd(), ndim::Span<uint8_t *>{ qry.data(), 8 }, LessThan);
    ASSERT_EQ(it, terminals.constEnd());
}

TEST_F(GraphFlushingTests, UniqueKmers) {
    auto range = makeRange();
    auto kmerEdges = flush::ranges::counting::uniqueEdges(range.dataBegin(), range.dataEnd(), range.dataB(), 10000);
    ASSERT_EQ(kmerEdges, 5);
}

TEST_F(GraphFlushingTests, UniqueTerminals) {
    auto range = makeRange();
    auto tmlEdges = flush::ranges::counting::uniqueEdges(range.terminalBegin(), range.terminalEnd(), range.terminalB(), 10000);
    ASSERT_EQ(tmlEdges, 9);
}

TEST_F(GraphFlushingTests, InitialValues) {
    auto range = makeRange();

    ASSERT_EQ(range.dataSize(), 8);
    ASSERT_EQ(range.terminalSize(), 9);
    ASSERT_EQ(range.edgePosition(), 0);
}

TEST_F(GraphFlushingTests, StepTo) {
    auto R1 = makeRange();
    [[maybe_unused]] decltype(R1) R2(R1, R1.dataSize(), flush::ranges::detail::pivot_tag);

    ASSERT_EQ(R1.dataSize(), 0);
    ASSERT_EQ(R1.terminalSize(), 0);
    ASSERT_EQ(R1.edgePosition(), 14);
}

TEST_F(GraphFlushingTests, InitialiseRange) {
    auto range = makeRange();

    flush::ranges::detail::SplitRange Split(&range);
    ASSERT_EQ(range.dataSize(), 8);
    ASSERT_EQ(range.terminalSize(), 9);
    ASSERT_EQ(range.edgePosition(), 0);

    flush::ranges::detail::TbbFree{}(Split.nextRange());
    ASSERT_EQ(range.dataSize(), 6);
    ASSERT_EQ(range.terminalSize(), 4);
    ASSERT_EQ(range.edgePosition(), 6);

    flush::ranges::detail::TbbFree{}(Split.nextRange());
    ASSERT_EQ(range.dataSize(), 1);
    ASSERT_EQ(range.terminalSize(), 2);
    ASSERT_EQ(range.edgePosition(), 11);

    flush::ranges::detail::TbbFree{}(Split.nextRange());
    ASSERT_EQ(range.dataSize(), 0);
    ASSERT_EQ(range.terminalSize(), 0);
    ASSERT_EQ(range.edgePosition(), 14);
}

TEST_F(GraphFlushingTests, FlushRange) {
    auto terms = terminals.asRange();
    size_t
        numEdges = 14,
        numColourValues = 17;

    graph::WriteableGraph graph(k, num_colours, num_colours, numEdges);
    auto crs = graph.from(0);

    flush::flush(crs, kmers, terms, 1, (uint8_t)4);

    CheckBufferEdgePositions(graph.colourBuffers);
    auto colours = graph.colourBuffers.combineAllBuffers();

    // check graph structure

    sdsl::int_vector<4> edges = { 0b1001, 0b1011, 0b1100, 0b1001, 0b1001, 0b1010, 0b1001, 0b1100, 0b1100, 0b1001, 0b0001, 0b1010, 0b0010, 0b1100 };
    sdsl::bit_vector wplus =     { 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    ASSERT_THAT(graph.edges_ref(), ContainerEq(edges));
    ASSERT_THAT(graph.wplus_ref(), ContainerEq(wplus));

    auto expColours = graph::ColourBuffer(graph.getColourBitWidth(), numColourValues);
    expColours.colours.assign(   { 0, 0, 0, 0, 0, 0, 1, 0, 0, 3, 0, 7, 1, 5, 0, 0, 0 });
    expColours.boundaries.assign({ 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1 });

    ASSERT_THAT(colours.colours, ContainerEq(expColours.colours));
    ASSERT_THAT(colours.boundaries, ContainerEq(expColours.boundaries));

    // Note that `F` and `C` arrays are not checked here because the terminals have a suffix removed
    // (in order to test for edge cases), but the test structure does not account for this.
    EXPECT_EQ(crs.position(), numEdges);

    // remove temp buffer file
    fs::remove(graph.colourBuffers.filePath());
}

/*

NNN A 1 0 1
NNN C 2 0 1
NNN G 3 0 1
NNN T 4 0 1
NNA G 3 0 1
ACA C 2 1 1
TCA N 0 0 1
NGA C 2 0 1
CGA G 3 2 1
TGA G 3 3 1
NTA C 2 0 1
GTA N 0 0 1
NNC G 3 0 1
CAC T 4 1 1
GAC T 4 3 1
TAC A 1 1 1
TAC T 4 2 1
CTC A 1 3 1
CTC G 3 2 1
GTC N 0 0 1
NNG A 1 0 1
NAG T 4 0 1
GAG T 4 2 1
GAG T 4 3 0
NCG A 1 0 1
TCG N 0 0 1
NTG A 1 0 1
GTG T 4 1 1
NNT A 1 0 1
NNT G 3 0 1
ACT C 2 2 1
ACT C 2 3 0
AGT A 1 2 1
AGT C 2 3 1
AGT G 3 1 1
TGT A 1 1 1

*/
class EgidiGraphTests : public testing::Test {
protected:
    EgidiGraphTests()
        : k(3)
        , genomes(loadGenomes({ STRING(CONTIG_ONE_PATH), STRING(CONTIG_TWO_PATH), STRING(CONTIG_THREE_PATH) }, k))
        , maxColour(genomes.numGenomes())
        , numColourValues(detail::calculateNumColourValuesBound(genomes, k))
        , numEdges(detail::calculateNumEdgesBound(genomes))
    {}

    uint8_t k;
    Dna4GenomeVector genomes;
    size_t 
        maxColour,
        numColourValues,
        numEdges;
};

TEST_F(EgidiGraphTests, BulkFill) {
    graph::WriteableGraph col_graph(k, maxColour, genomes.numGenomes(), numEdges);
    auto p = col_graph.from(0);
    fillGraphBulk(p, genomes, 2);
    col_graph.resize(p);

    // check graph structure

    sdsl::int_vector<4> edges = { 0b1001, 0b1010, 0b1011, 0b1100, 0b1011, 0b1010, 0b1000, 0b1010, 0b1011, 0b0011, 0b1010, 0b1000, 0b1011, 0b1100, 0b0100, 0b1001, 0b0100, 0b1001, 0b1011, 0b1000, 0b1001, 0b1100, 0b0100, 0b1001, 0b1000, 0b1001, 0b1100, 0b1001, 0b1011, 0b1010, 0b1001, 0b1010, 0b1011, 0b0001 };
    sdsl::bit_vector wplus     = { 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1 };

    ASSERT_THAT(col_graph.edges_ref(), ContainerEq(edges));
    ASSERT_THAT(col_graph.wplus_ref(), ContainerEq(wplus));

    CheckBufferEdgePositions(col_graph.colourBuffers);
    auto colours = col_graph.colourBuffers.combineAllBuffers();

    auto expColours = graph::ColourBuffer(col_graph.getColourBitWidth(), 36);
    expColours.colours.assign(   { 0, 0, 0, 0, 0, 1, 0, 0, 2, 3, 0, 0, 0, 1, 3, 1, 2, 3, 2, 0, 0, 0, 2, 3, 0, 0, 0, 1, 0, 0, 2, 3, 2, 3, 1, 1 });
    expColours.boundaries.assign({ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1 });

    ASSERT_THAT(colours.colours, ContainerEq(expColours.colours));
    ASSERT_THAT(colours.boundaries, ContainerEq(expColours.boundaries));

    ASSERT_THAT(p.C, ElementsAreArray({ 1, 8, 6, 7, 4 }));
    ASSERT_THAT(p.F, ElementsAreArray({ 4, 8, 8, 7, 7 }));
    ASSERT_THAT(p.position(), 34);

    // remove temp buffer file
    fs::remove(col_graph.colourBuffers.filePath());
}

TEST_F(EgidiGraphTests, SuffixFill_1) {
    graph::WriteableGraph col_graph(k, maxColour, genomes.numGenomes(), numEdges);
    uint8_t s = 1;
    auto p = col_graph.from(0);
    fillGraphBySuffix(p, genomes, s, 2);
    col_graph.resize(p);
    
    // check graph structure

    ASSERT_THAT(p.C, ElementsAreArray({ 1, 8, 6, 7, 4 }));
    ASSERT_THAT(p.F, ElementsAreArray({ 4, 8, 8, 7, 7 }));
    ASSERT_THAT(p.position(), 34);

    sdsl::int_vector<4> edges = { 0b1001, 0b1010, 0b1011, 0b1100, 0b1011, 0b1010, 0b1000, 0b1010, 0b1011, 0b0011, 0b1010, 0b1000, 0b1011, 0b1100, 0b0100, 0b1001, 0b0100, 0b1001, 0b1011, 0b1000, 0b1001, 0b1100, 0b0100, 0b1001, 0b1000, 0b1001, 0b1100, 0b1001, 0b1011, 0b1010, 0b1001, 0b1010, 0b1011, 0b0001 };
    sdsl::bit_vector wplus     = { 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1 };

    ASSERT_THAT(col_graph.edges_ref(), ContainerEq(edges));
    ASSERT_THAT(col_graph.wplus_ref(), ContainerEq(wplus));

    CheckBufferEdgePositions(col_graph.colourBuffers);
    auto colours = col_graph.colourBuffers.combineAllBuffers();

    auto expColours = graph::ColourBuffer(col_graph.getColourBitWidth(), 36);
    expColours.colours.assign(   { 0, 0, 0, 0, 0, 1, 0, 0, 2, 3, 0, 0, 0, 1, 3, 1, 2, 3, 2, 0, 0, 0, 2, 3, 0, 0, 0, 1, 0, 0, 2, 3, 2, 3, 1, 1 });
    expColours.boundaries.assign({ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1 });

    ASSERT_THAT(colours.colours, ContainerEq(expColours.colours));
    ASSERT_THAT(colours.boundaries, ContainerEq(expColours.boundaries));

    // remove temp buffer file
    fs::remove(col_graph.colourBuffers.filePath());
}

class EgidiGraphGenomeTests : public testing::Test {
protected:
    EgidiGraphGenomeTests()
        : k(3)
        , genome(readEgidiGenome(k))
        , numEdges(detail::calculateNumEdgesBound(genome))
    {
    }

    static ChunkedDna4Genome readEgidiGenome(size_t k) {
        std::string data = ">contig_one\nTACACT\n>contig_two\nTACTCG\n>contig_three\nGACTCA\n";
        std::istringstream datastream(data);

        Dna4Genome genome;
        parseFastaStream(genome, datastream, 0, k);

        return ChunkedDna4Genome(std::move(genome), 10, k);
    }

    uint8_t k;
    ChunkedDna4Genome genome;
    size_t numEdges;
};

TEST_F(EgidiGraphGenomeTests, BulkFill) {
    graph::BaseGraph g(k, 1, numEdges);
    auto p = g.from(0);
    fillGraphBulk(p, genome, 2);
    g.resize(p);

    // check graph structure

    ASSERT_THAT(p.C, ElementsAreArray({ 1, 8, 6, 7, 4 }));
    ASSERT_THAT(p.F, ElementsAreArray({ 4, 8, 8, 7, 7 }));
    ASSERT_THAT(p.position(), 34);

    sdsl::int_vector<4> edges = { 0b1001, 0b1010, 0b1011, 0b1100, 0b1011, 0b1010, 0b1000, 0b1010, 0b1011, 0b0011, 0b1010, 0b1000, 0b1011, 0b1100, 0b0100, 0b1001, 0b0100, 0b1001, 0b1011, 0b1000, 0b1001, 0b1100, 0b0100, 0b1001, 0b1000, 0b1001, 0b1100, 0b1001, 0b1011, 0b1010, 0b1001, 0b1010, 0b1011, 0b0001 };
    sdsl::bit_vector wplus     = { 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1 };

    ASSERT_THAT(g.edges_ref(), ContainerEq(edges));
    ASSERT_THAT(g.wplus_ref(), ContainerEq(wplus));
}

TEST_F(EgidiGraphGenomeTests, SuffixFill_1) {
    graph::BaseGraph g(k, 1, numEdges);
    uint8_t s = 1;
    auto p = g.from(0);
    fillGraphBySuffix(p, genome, s, 2);
    g.resize(p);
    
    // check graph structure

    ASSERT_THAT(p.C, ElementsAreArray({ 1, 8, 6, 7, 4 }));
    ASSERT_THAT(p.F, ElementsAreArray({ 4, 8, 8, 7, 7 }));
    ASSERT_THAT(p.position(), 34);

    sdsl::int_vector<4> edges = { 0b1001, 0b1010, 0b1011, 0b1100, 0b1011, 0b1010, 0b1000, 0b1010, 0b1011, 0b0011, 0b1010, 0b1000, 0b1011, 0b1100, 0b0100, 0b1001, 0b0100, 0b1001, 0b1011, 0b1000, 0b1001, 0b1100, 0b0100, 0b1001, 0b1000, 0b1001, 0b1100, 0b1001, 0b1011, 0b1010, 0b1001, 0b1010, 0b1011, 0b0001 };
    sdsl::bit_vector wplus     = { 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1 };
}

class EgidiGraphReadTests : public testing::Test {
protected:
    EgidiGraphReadTests()
        : k(3)
        , data(3)
        , numEdges()
    {
        data[0]
            .emplace_back()
            .emplace_back("contig_one")
            .push("TACACT"_dna5, k);
        data[1]
            .emplace_back()
            .emplace_back("contig_two")
            .push("TACTCG"_dna5, k);
        data[2]
            .emplace_back()
            .emplace_back("contig_three")
            .push("GACTCA"_dna5, k);
        numEdges = detail::calculateNumEdgesBound(data);
    }

    uint8_t k;
    std::vector<reads::ReadChunks> data;
    size_t numEdges;
};

TEST_F(EgidiGraphReadTests, BulkFill) {
    graph::BaseGraph g(k, 1, numEdges);
    auto p = g.from(0);
    fillGraphBulk(p, data, 2);
    g.resize(p);

    // check graph structure
    
    ASSERT_THAT(p.position(), 34);
    ASSERT_THAT(p.C, ElementsAreArray({ 1, 8, 6, 7, 4 }));
    ASSERT_THAT(p.F, ElementsAreArray({ 4, 8, 8, 7, 7 }));

    sdsl::int_vector<4> edges = { 0b1001, 0b1010, 0b1011, 0b1100, 0b1011, 0b1010, 0b1000, 0b1010, 0b1011, 0b0011, 0b1010, 0b1000, 0b1011, 0b1100, 0b0100, 0b1001, 0b0100, 0b1001, 0b1011, 0b1000, 0b1001, 0b1100, 0b0100, 0b1001, 0b1000, 0b1001, 0b1100, 0b1001, 0b1011, 0b1010, 0b1001, 0b1010, 0b1011, 0b0001 };
    sdsl::bit_vector wplus    = { 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1 };

    ASSERT_THAT(g.edges_ref(), ContainerEq(edges));
    ASSERT_THAT(g.wplus_ref(), ContainerEq(wplus));
}

TEST_F(EgidiGraphReadTests, SuffixFill_1) {
    graph::BaseGraph g(k, 1, numEdges);
    uint8_t s = 1;
    auto p = g.from(0);
    fillGraphBySuffix(p, data, s, 2);
    g.resize(p);
    
    // check graph structure

    ASSERT_THAT(p.C, ElementsAreArray({ 1, 8, 6, 7, 4 }));
    ASSERT_THAT(p.F, ElementsAreArray({ 4, 8, 8, 7, 7 }));
    ASSERT_THAT(p.position(), 34);

    sdsl::int_vector<4> edges = { 0b1001, 0b1010, 0b1011, 0b1100, 0b1011, 0b1010, 0b1000, 0b1010, 0b1011, 0b0011, 0b1010, 0b1000, 0b1011, 0b1100, 0b0100, 0b1001, 0b0100, 0b1001, 0b1011, 0b1000, 0b1001, 0b1100, 0b0100, 0b1001, 0b1000, 0b1001, 0b1100, 0b1001, 0b1011, 0b1010, 0b1001, 0b1010, 0b1011, 0b0001 };
    sdsl::bit_vector wplus     = { 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1 };
}

// larger k-mer and suffix size

class SmallFastaTests : public testing::Test {
protected:
    SmallFastaTests() : 
        k(9),
        genomes(loadGenomes({ STRING(SMALL_SEQ) }, k)),
        num_colours(genomes.numGenomes()),
        total_bp(detail::calculateNumEdgesBound(genomes))
    {}

    uint8_t k;
    Dna4GenomeVector genomes;
    size_t num_colours, total_bp;
};

TEST_F(SmallFastaTests, SuffixFill_3) {
    auto col_graph = graph::WriteableGraph(k, num_colours, genomes.numGenomes(), total_bp);
    uint8_t s = 3;
    auto p = col_graph.from(0);
    fillGraphBySuffix(p, genomes, s, 2);
    col_graph.resize(p);

    // check graph structure

    sdsl::int_vector<4> edges = { 0b1010, 0b1011, 0b1011, 0b1011, 0b1010, 0b1100, 0b1010, 0b1011, 0b1010, 0b1100, 0b1010, 0b1100, 0b1100, 0b1011, 0b1001, 0b1100, 0b1001, 0b1011, 0b1010, 0b1010, 0b1100, 0b1011, 0b1001, 0b1011, 0b1011, 0b1001, 0b1001, 0b1011, 0b1010, 0b1001, 0b1011, 0b1100, 0b1000, 0b1010, 0b1001, 0b1001, 0b1100, 0b1011, 0b1100, 0b1001, 0b1011, 0b1001, 0b1001, 0b1001, 0b1100, 0b1011, 0b1010, 0b1100, 0b1100, 0b1100, 0b1001, 0b1100, 0b1010, 0b1011, 0b1010, 0b1001, 0b1011, 0b1010, 0b1011, 0b1011, 0b1100, 0b1010, 0b1100, 0b1001, 0b1100, 0b1010, 0b1010, 0b1010, 0b1011, 0b1010, 0b1010, 0b1100, 0b1001, 0b1010, 0b1001, 0b1011, 0b1001, 0b1100, 0b1001, 0b1011, 0b1010, 0b1100, 0b1001, 0b1001, 0b1001, 0b1010, 0b1001, 0b1011, 0b1100, 0b1011, 0b1010, 0b1000, 0b1100, 0b1011, 0b1010, 0b1011, 0b1011, 0b1010, 0b1011, 0b1100, 0b1011, 0b1100, 0b1001, 0b1010, 0b1010, 0b1011, 0b1010, 0b1010, 0b1011, 0b1100, 0b1011, 0b1100, 0b1001, 0b1100, 0b1011, 0b1010, 0b1010, 0b1001, 0b1010, 0b1011, 0b1001, 0b1010 };
    sdsl::bit_vector wplus =     { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    
    ASSERT_THAT(col_graph.edges_ref(), ContainerEq(edges));
    ASSERT_THAT(col_graph.wplus_ref(), ContainerEq(wplus));

    CheckBufferEdgePositions(col_graph.colourBuffers);
    auto colours = col_graph.colourBuffers.combineAllBuffers();

    graph::ColourBuffer expColours(col_graph.getColourBitWidth(), p.position());
    expColours.colours.assign({ 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1 });
    sdsl::util::set_to_value(expColours.boundaries, 1);

    ASSERT_THAT(colours.colours, ContainerEq(expColours.colours));
    ASSERT_THAT(colours.boundaries, ContainerEq(expColours.boundaries));

    ASSERT_THAT(p.C, ElementsAreArray({ 1, 27, 33, 33, 27 }));
    ASSERT_THAT(p.F, ElementsAreArray({ 2, 27, 33, 33, 27 }));
    ASSERT_THAT(p.position(), 122);

    // remove temp buffer file
    fs::remove(col_graph.colourBuffers.filePath());
}

class GraphIOTests : public testing::Test {
protected:
    GraphIOTests()
        :
            k(9),
            maxColour(7),
            cwidth(3),
            edges({0b1000, 0b1001, 0b0010, 0b1011, 0b1100, 0b0000, 0b1001, 0b1010, 0b1011, 0b0100}),
            wplus({0, 0, 1, 0, 0, 1, 0, 0, 1, 0}),
            F({0, 1, 2, 3, 4}),
            C({4, 3, 2, 1, 0}),
            testGraph(k, maxColour, maxColour, sdsl::int_vector<4U>(edges), sdsl::bit_vector(wplus), graph::ColourBufferRegistry(), F, C),
            tempTestDir(tempio::create_temporary_directory())
    {}

    ~GraphIOTests()
    {
        LOG(INFO) << "removing temporary directory " << tempTestDir;
        fs::remove_all(tempTestDir);
    }

    uint8_t
        k,
        maxColour,
        cwidth;
    sdsl::int_vector<4> edges;
    sdsl::bit_vector wplus;
    std::array<size_t, 5> F,
                          C;
    graph::WriteableGraph testGraph;
    fs::path tempTestDir;
};

TEST_F(GraphIOTests, LoadFromDisk) {
    graph::saveToDisk(testGraph, tempTestDir);

    graph::WriteableGraph loadedGraph;
    graph::loadFromDisk(tempTestDir, loadedGraph);

    EXPECT_EQ(loadedGraph.kmer_size(), k);
    EXPECT_EQ(loadedGraph.getColourBitWidth(), cwidth);

    ASSERT_THAT(loadedGraph.edges_ref(), ContainerEq(edges));
    ASSERT_THAT(loadedGraph.wplus_ref(), ContainerEq(wplus));

    EXPECT_EQ(loadedGraph.C, C);
    EXPECT_EQ(loadedGraph.F, F);
}

// main

int main(int argc, char **argv) {
    FLAGS_logtostderr = USE_LOGS;
    google::InitGoogleLogging("test_graph");
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
