#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <maki/maki.h>
#include <maki/fasta.hpp>
#include <maki/graph.hpp>
#include <maki/dist.hpp>
#include <maki/classify.hpp>

using namespace seqan3::literals;

using ::testing::ContainerEq;
using ::testing::ElementsAreArray;

// ==== Utilities to support testing =================================================

/**
 * Reverse complement of sequence.
 */
std::string ReverseComplement(const std::string& seq) {
    std::string rcomp = seq;
    // Reverse the sequence
    std::reverse(rcomp.begin(), rcomp.end());
    // Complement the bases
    for (char& nt : rcomp) {
        switch (nt) {
            case 'A': nt = 'T'; break;
            case 'T': nt = 'A'; break;
            case 'C': nt = 'G'; break;
            case 'G': nt = 'C'; break;
        }
    }
    return rcomp;
}

/**
 * Parse FASTA text into a vector of sequences. Lines starting with '>' are headers.
 * Sequences can span multiple lines; boundaries DO NOT cross (each sequence is independent).
 */
static std::vector<std::string> ParseFastaSequences(const std::string& fasta) {
    std::vector<std::string> sequences;
    std::string current;
    std::string line;
    for (size_t i = 0, n = fasta.size(); i <= n; ++i) {
        char c = (i < n ? fasta[i] : '\n');
        if (c == '\r') continue; // normalize CRLF
        if (c != '\n') {
            line.push_back(c);
        } else {
            if (!line.empty()) {
                if (!line.empty() && line[0] == '>') {
                    if (!current.empty()) {
                        sequences.push_back(current);
                        current.clear();
                    }
                    // header line; ignore
                } else {
                    // append sequence line (strip spaces, to upper)
                    for (char ch : line) {
                        if (!std::isspace(static_cast<unsigned char>(ch))) {
                            current.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
                        }
                    }
                }
                line.clear();
            } else {
                // blank line: ignore
            }
        }
    }
    if (!current.empty()) {
        sequences.push_back(current);
    }
    return sequences;
}

/**
 * Compute all UNIQUE k-mers for a single sequence (does not cross boundaries).
 */
static std::unordered_set<std::string> UniqueKMersForSequence(const std::string& seq, size_t k) {
    std::unordered_set<std::string> s;
    if (k == 0 || seq.size() < k) return s;
    for (size_t i = 0; i + k <= seq.size(); ++i) {
        s.insert(seq.substr(i, k));
    }
    return s;
}

/**
 * Compute all UNIQUE k-mers from a FASTA (union across sequences; no cross-boundary k-mers).
 */
static std::unordered_set<std::string> UniqueKMersFromFasta(const std::string& fasta, size_t k) {
    auto seqs = ParseFastaSequences(fasta);
    std::unordered_set<std::string> all;
    for (const auto& s : seqs) {
        // Forward Sequences
        {
            auto ks = UniqueKMersForSequence(s, k);
            all.insert(ks.begin(), ks.end());
        }
        // Reverse Complement
        {
            auto rks = UniqueKMersForSequence(ReverseComplement(s), k);
            all.insert(rks.begin(), rks.end());
        }
    }
    return all;
}

/**
 * Helper to compute edge count.
 */
static size_t edge_count(const graph::TraversableGraph& g) {
    return g.edge_count();
}

/**
 * Build a FASTA string from sequences with trivial generated headers.
 */
static std::string BuildFasta(const std::vector<std::string>& seqs) {
    std::string fasta;
    for (size_t i = 0; i < seqs.size(); ++i) {
        fasta += ">seq" + std::to_string(i) + "\n";
        fasta += seqs[i] + "\n";
    }
    return fasta;
}

/**
 * Random DNA generator (A/C/G/T) with reproducible seed.
 */
static std::string RandomDNA(size_t length, std::mt19937& rng) {
    static const char bases[4] = {'A','C','G','T'};
    std::uniform_int_distribution<int> dist(0, 3);
    std::string s;
    s.reserve(length);
    for (size_t i = 0; i < length; ++i) s.push_back(bases[dist(rng)]);
    return s;
}

Dna4GenomeVector readGenome(std::string &data, size_t k) {
    Dna4GenomeVector result(1);
    std::istringstream datastream(data);
    Dna4Genome &genome = result[0];
    parseFastaStream(genome, datastream, 1, k);
    result.setMaxFeatureId(getMaxFeatureId(result));
    result.setBaseFeatureWidth(1);
    return result;
}

graph::TraversableGraph makeGraph(const Dna4GenomeVector &fna, size_t k, fs::path &bufferPath) {
    graph::WriteableGraph g
        = graph::initialiseEmptyGraph(fna, k-1, detail::calculateNumEdgesBound(fna), bufferPath);
    graph::makeGraphBulk(g, fna, 1);
    return graph::TraversableGraph(std::move(g));
}

graph::TraversableGraph makeGraph(std::string data, size_t k, fs::path &bufferPath) {
    Dna4GenomeVector fna = readGenome(data, k);
    return makeGraph(fna, k, bufferPath);
}

// ==== Core checker that targets both correctness and bwd sanity ====================

/**
 * Build graph from FASTA & k, then:
 *  - Check num_edges matches expected unique k-mers.
 *  - For each edge i:
 *      - kmer(i) length == k
 *      - kmer(i) is one of expected k-mers
 *      - bwd chain never yields an ID >= num_edges
 *  - Check the set of returned kmer(i) across all edges equals the expected set.
 */
static void BuildAndCheckGraph(const std::string& fasta, size_t k, fs::path &bufferPath) {
    graph::TraversableGraph g = makeGraph(fasta, k, bufferPath);

    const size_t E = edge_count(g);
    const auto expected_set = UniqueKMersFromFasta(fasta, k);

    // If your BOSS includes sentinel ($) edges or duplicates, this equality may need relaxing.
    // Otherwise, assert exact match between unique k-mers and edges stored.
    ASSERT_GE(E, expected_set.size()) << "Edge count != unique k-mers in FASTA; adjust if graph retains multiplicity/sentinels.";

    std::unordered_set<std::string> observed;
    observed.reserve(E);

    for (size_t i = 0; i < E; ++i) {
        std::string km = toString(g.kmer(i));
        if (km.size() < k) continue; // Skip terminals

        ASSERT_EQ(km.size(), k) << "kmer(" << i << ") has length " << km.size() << ", expected " << k;
        ASSERT_TRUE(expected_set.find(km) != expected_set.end())
            << "kmer(" << i << ")='" << km << "' not found among expected k-mers.";

        // Targeted bwd sanity check: for up to k-1 steps, ensure any returned ID is < E.
        size_t curr = i;
        for (size_t step = 0; step + 1 < k; ++step) {  // up to k-1 backwards steps
            std::optional<size_t> prev = g.bwd(curr);
            if (!prev.has_value()) break;  // reached a boundary; valid
            ASSERT_LT(prev.value(), E) << "bwd(" << curr << ") returned out-of-range edge id " << prev.value()
                                       << " (E=" << E << ").";
            curr = prev.value();
        }

        observed.insert(km);
    }

    // Ensure bijection: all expected k-mers appear as edge labels
    ASSERT_EQ(observed.size(), expected_set.size());
    for (const auto& km : expected_set) {
        ASSERT_TRUE(observed.find(km) != observed.end())
            << "Expected k-mer '" << km << "' not found among graph edges.";
    }
}

/**
 * Check buffer structure.
 */
static void CheckBufferEdgePositions(const graph::ColourBufferRegistry &buffers) {
    std::cout << "checking " << buffers.numBuffers() << " colour buffers" << std::endl;
    ZstdDecompressor zstd;
    // trace edges
    size_t edgeCount = 0;
    for (const auto &[key, bufferId] : buffers) {
        graph::ColourBuffer src = buffers.readColours(bufferId, zstd);
        for (auto b : src.boundaries) edgeCount += b;
        ASSERT_EQ(edgeCount, key);
    }
}

// ==== Tests ========================================================================

class TraversableGraph_Kmer : public testing::Test {
protected:
    TraversableGraph_Kmer()
        : bufferPath(tempio::create_temporary_directory())
        , graphBuffer(bufferPath / "ref.h5")
    {}

    ~TraversableGraph_Kmer() { fs::remove_all(bufferPath); }

    fs::path bufferPath, graphBuffer;
};

TEST_F(TraversableGraph_Kmer, SimplePath_NoRepeats) {
    std::string fasta = BuildFasta({"ACGTACGT"});
    size_t k = 3;
    BuildAndCheckGraph(fasta, k, graphBuffer);
}

TEST_F(TraversableGraph_Kmer, Homopolymer_Repeats) {
    std::string fasta = BuildFasta({"AAAAAAAAAA"});  // 10 As
    for (size_t k = 2; k <= 9; ++k) {
        BuildAndCheckGraph(fasta, k, graphBuffer);
    }
}

TEST_F(TraversableGraph_Kmer, Branching) {
    std::string fasta = BuildFasta({"ACAAA", "ACGAA"});
    size_t k = 3;
    BuildAndCheckGraph(fasta, k, graphBuffer);
}

TEST_F(TraversableGraph_Kmer, Palindromic) {
    std::string fasta = BuildFasta({"ATGCAT"});
    size_t k = 3;
    BuildAndCheckGraph(fasta, k, graphBuffer);
}

TEST_F(TraversableGraph_Kmer, MultipleSequences_BoundariesRespected) {
    std::string fasta = BuildFasta({"AAA", "TTT"});
    size_t k = 2;
    BuildAndCheckGraph(fasta, k, graphBuffer);
}

TEST_F(TraversableGraph_Kmer, KEqualsSequenceLength) {
    std::string fasta = BuildFasta({"ACGT"});
    size_t k = 4;
    BuildAndCheckGraph(fasta, k, graphBuffer);
}

TEST_F(TraversableGraph_Kmer, KGreaterThanSequenceLength_ShouldYieldNoEdges) {
    std::string fasta = BuildFasta({"ACGT"});
    size_t k = 5;
    graph::TraversableGraph g = makeGraph(fasta, k, graphBuffer);
    ASSERT_EQ(edge_count(g), 0u) << "Graph should have 0 edges when k > sequence length.";
}

TEST_F(TraversableGraph_Kmer, RandomizedFuzz_Small) {
    std::mt19937 rng(123456); // reproducible
    std::uniform_int_distribution<int> numSeqDist(1, 3);
    std::uniform_int_distribution<int> lenDist(10, 100);

    for (int iter = 0; iter < 50; ++iter) {
        int numSeqs = numSeqDist(rng);
        std::vector<std::string> seqs;
        seqs.reserve(numSeqs);
        size_t minLen = SIZE_MAX;
        for (int s = 0; s < numSeqs; ++s) {
            size_t len = static_cast<size_t>(lenDist(rng));
            seqs.push_back(RandomDNA(len, rng));
            minLen = std::min(minLen, len);
        }
        std::string fasta = BuildFasta(seqs);

        std::uniform_int_distribution<size_t> kDist(2, std::min<size_t>(10, minLen));
        size_t k = kDist(rng);

        BuildAndCheckGraph(fasta, k, graphBuffer);
    }
}

TEST_F(TraversableGraph_Kmer, Stress_LongSequence_ModerateK) {
    std::mt19937 rng(424242);
    std::string seq = RandomDNA(2000, rng);
    std::string fasta = BuildFasta({seq});
    size_t k = 21;

    BuildAndCheckGraph(fasta, k, graphBuffer);
}

/*

   Node    Edge    W(-)    W(+)
 0 NNN     G       1       0
 1 NNN     T       1       1
 2 TCA     N       1       1
 3 NGA     C       1       1
 4 TGA     G       1       1
 5 GAC     T       1       1
 6 CTC     A       1       1
 7 GTC     N       1       1
 8 NNG     A       1       1
 9 GAG     T       1       1
10 NTG     A       1       1
11 NNT     G       1       1
12 ACT     C       1       1
13 AGT     C       1       1

*/
class GraphTraverseTests : public testing::Test {
protected:
    GraphTraverseTests()
        : k(4)
        , bufferPath(tempio::create_temporary_directory())
        , graphBuffer(bufferPath / "ref.h5")
        , g(makeGraph(">contig_one\nGACTCA\n>contig_two\nTGAGTC\n", k, graphBuffer))
    {}

    ~GraphTraverseTests() { fs::remove_all(bufferPath); }

    uint8_t k;
    fs::path bufferPath, graphBuffer;
    graph::TraversableGraph g;
};

TEST_F(GraphTraverseTests, Fwd0) {
    auto i = g.fwd(0);
    ASSERT_EQ(*i, 8);
}

TEST_F(GraphTraverseTests, Fwd1) {
    auto i = g.fwd(1);
    ASSERT_EQ(*i, 11);
}

TEST_F(GraphTraverseTests, Fwd4) {
    auto i = g.fwd(4);
    ASSERT_EQ(*i, 9);
}

TEST_F(GraphTraverseTests, NoFwd) {
    auto i = g.fwd(2);
    ASSERT_FALSE(i.has_value());
}

TEST_F(GraphTraverseTests, Bwd9) {
    auto i = g.bwd(9);
    ASSERT_EQ(*i, 4);
}

TEST_F(GraphTraverseTests, BwdRoot) {
    auto i = g.bwd(8);
    ASSERT_EQ(*i, 0);
}

TEST_F(GraphTraverseTests, BwdNearRoot) {
    auto i = g.bwd(11);
    ASSERT_EQ(*i, 1);
}

TEST_F(GraphTraverseTests, Bwd0) {
    auto i = g.bwd(0);
    ASSERT_FALSE(i.has_value());
}

TEST_F(GraphTraverseTests, Kmer6) {
    auto s = g.kmer(6);
    ASSERT_EQ(s, Dna4Sequence("CTCA"_dna4));
}

TEST_F(GraphTraverseTests, Kmer7) {
    auto s = g.kmer(7);
    ASSERT_EQ(s, Dna4Sequence("GTC"_dna4));
}

TEST_F(GraphTraverseTests, Kmer8) {
    auto s = g.kmer(8);
    ASSERT_EQ(s, Dna4Sequence("GA"_dna4));
}

TEST_F(GraphTraverseTests, Kmer0) {
    auto s = g.kmer(0);
    ASSERT_EQ(s, Dna4Sequence("G"_dna4));
}

/*

   Node    Edge    W(-)    W(+)
 0 NNN     A       1       0
 1 NNN     C       1       1
 2 NNA     C       1       1
 3 CGA     C       1       1
 4 NNC     C       1       1
 5 NAC     G       1       1
 6 GAC     C       1       1
 7 NCC     G       1       1
 8 ACC     G       0       1
 9 GTC     G       1       1
10 ACG     G       1       1
11 CCG     A       1       0
12 CCG     T       1       1
13 TCG     G       0       1
14 CGG     T       1       1
15 GGT     C       1       1

*/
class GraphTraverseRepeatTests : public testing::Test {
protected:
    GraphTraverseRepeatTests()
        : k(4)
        , bufferPath(tempio::create_temporary_directory())
        , graphBuffer(bufferPath / "other.h5")
        , g(makeGraph(">contig_one\nACGGTCGG\n", k, graphBuffer))
    {}

    ~GraphTraverseRepeatTests() { fs::remove_all(bufferPath); }

    uint8_t k;
    fs::path bufferPath, graphBuffer;
    graph::TraversableGraph g;
};

TEST_F(GraphTraverseRepeatTests, Fwd7) {
    auto i = g.fwd(7);
    ASSERT_EQ(*i, 12);
}

TEST_F(GraphTraverseRepeatTests, Fwd8) {
    auto i = g.fwd(8);
    ASSERT_EQ(*i, 12);
}

TEST_F(GraphTraverseRepeatTests, Bwd12) {
    auto i = g.bwd(12);
    ASSERT_EQ(*i, 7);
}

TEST_F(GraphTraverseRepeatTests, Bwd11) {
    auto i = g.bwd(11);
    ASSERT_EQ(*i, 7);
}

TEST_F(GraphTraverseRepeatTests, Kmer8) {
    auto s = g.kmer(8);
    ASSERT_EQ(s, Dna4Sequence("ACCG"_dna4));
}

TEST_F(GraphTraverseRepeatTests, Kmer8Str) {
    auto s = toString(g.kmer(8));
    ASSERT_EQ(s, std::string("ACCG"));
}

class TraversableGraph_Colour : public testing::Test {
protected:
    TraversableGraph_Colour()
        : bufferPath(tempio::create_temporary_directory())
        , graphBuffer(bufferPath / "ref.h5")
    {}

    ~TraversableGraph_Colour() { fs::remove_all(bufferPath); }

    fs::path bufferPath, graphBuffer;
};

/**
 * Collect requested colours from graph
 */
TEST_F(TraversableGraph_Colour, RandomSequences) {
    // generate random genomes
    std::mt19937 rng(123456); // reproducible
    std::uniform_int_distribution<int> numSeqDist(5, 20);
    int numSeqs = numSeqDist(rng);
    auto fna = fasta_utils::randomGenomes(numSeqs, 42u, 100u, 10000u);

    // build graph
    size_t k = 3;
    auto g = makeGraph(fna, k, graphBuffer);
    size_t edges = edge_count(g);
    CheckBufferEdgePositions(g.colourBuffers);

    // request all colours
    classify::ResultVector req(edges);
    std::iota(req.begin(), req.end(), 0u);
    auto result = classify::retrieveColours(std::move(req), g);

    // count all colours manually
    auto manual = dist::cooc::countRaw(g, 1u);

    // compare
    ASSERT_EQ(result, manual);
}

