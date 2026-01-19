#include <filesystem>
#include <vector>
#include <glog/logging.h>
#include <gtest/gtest.h>
#include <maki/fasta.hpp>
#include <maki/graph.hpp>
#include <maki/maki.h>
#include <maki/utils.hpp>

namespace fs = std::filesystem;

class MainTests : public testing::Test {
protected:
    MainTests() : bufferFile((fs::temp_directory_path() / "__temp_ColourBufferRegistry.h5").string()) {}

    ~MainTests() { fs::remove(bufferFile); }

    uint8_t
        k = 15,
        s = 1;
    uint32_t threads = 10;
    size_t rr = 2;
    // genomes
    std::vector<fs::path> genome_files = {
        STRING(GCF_000005845_FNA),
        STRING(GCF_000006765_FNA),
        STRING(GCF_000006925_FNA),
        STRING(GCF_000006945_FNA),
        STRING(GCF_000007765_FNA),
        STRING(GCF_000008725_FNA)
    };
    Dna4GenomeVector genomes = loadGenomes(genome_files, k);
    std::string bufferFile;
};

TEST_F(MainTests, Bulk)
{
    auto garr = graph::initialiseEmptyGraph(
        genomes,
        k,
        detail::calculateNumEdgesBound(genomes),
        bufferFile);
    graph::makeGraphBulk(garr, genomes, threads);
}

TEST_F(MainTests, Suffix)
{
    auto garr = graph::initialiseEmptyGraph(
        genomes,
        k,
        detail::calculateNumEdgesBound(genomes),
        bufferFile);
    graph::makeGraphBySuffix(garr, genomes, k, s, threads);
}

TEST_F(MainTests, Detect)
{
    auto garr = graph::initialiseEmptyGraph(
        genomes,
        k,
        detail::calculateNumEdgesBound(genomes),
        bufferFile);
    graph::makeGraphDetect(
        garr,
        genomes,
        k,
        threads,
        rr * 1000000000ull,
        _mkimem_details::record_size(
            _mkimem_details::key_size(k),
            _mkimem_details::value_size(genomes.getFeatureWidth())));
}
