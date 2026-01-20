#pragma once
#include <array>
#include <string>
#include <vector>
#include <unordered_map>
#include <seqan3/argument_parser/all.hpp>
#include <glog/logging.h>
#include <gtl/phmap.hpp>
#include <oneapi/tbb.h>
#include <filesystem>
#include <maki/maki.h>
#include <maki/graph.hpp>
#include <maki/classify.hpp>
#include <maki/dist.hpp>
#include <maki/fasta.hpp>
#include <maki/tree.hpp>
#include <maki/timing.hpp>
#include <maki/utils.hpp>

// Base command interface
class BaseCommand {
public:
    virtual ~BaseCommand() =default;
    virtual void registerOptions(seqan3::argument_parser &parser) = 0;
    virtual int execute() = 0;
    int run(const char *execName, seqan3::argument_parser &parser);
protected:
    std::string logFile;
    size_t threads = 1;
};

// Command declarations
class IndexCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string queryFile, outputFolder;
    uint8_t kmerSize = 0;
    bool overrideDatabase = false;
    long double resourceLimit = 0.0;
    bool useBulk = false;
    uint8_t suffixSize = 0;
};

class MakeFilterCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string outputFolder;
    uint8_t kmerSize = 0;
    bool overrideDatabase = false;
    long double resourceLimit = 0.0;
    bool useBulk = false;
    uint8_t suffixSize = 0;
};

class DbStatsCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string databasePath;
    bool loadBig = false;
};

class ClassifyCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string databasePath, queryPath, outputFile;
    size_t minReadLength = 50;
    long double resourceLimit = 0.0;
    bool useBulk = false;
    uint8_t suffixSize = 0;
};

class DistanceCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string databasePath, outputFile, tempDir;
};

class CountCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string databasePath, outputFile;
    bool useWhole = false;
};

class CoverageCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string inputFile, countFile, outputFile;
    bool skipHeader;
    size_t valueColumn = 1u;
};

class SubcountCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string databasePath, countString, outputFile;
    size_t reps = 1;
    std::vector<size_t> parseSubsamples() const;
    std::vector<size_t> extractSubsamples(uint64_t) const;
};

class CompressCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string databasePath, outputFile;
};

class JaccardCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string databasePath, distFile, countFile, nwkFile, outputFile;
    size_t take = 0;
    double alpha = 0.05;
};

class LcaCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string databasePath, clusterFile, countFile, nwkFile, outputFile;
    bool printColours = false;
    float minFamDiv = 0.0, maxFamDiv = 1.0, minPhyCons = 0.0, maxPhyCons = 1.0;
};

class GlcaCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string databasePath, clusterFile, nwkFile, outputFile;
};

class SamdistCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string databasePath, countFile, nwkFile, outputFile;
};

class QueryCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string queryFile, databasePath, countFile, nwkFile, outputFile;
};

class KmersCommand : public BaseCommand {
public:
    void registerOptions(seqan3::argument_parser &parser) override;
    int execute() override;
private:
    std::string coloursFile, databasePath, outputFile;
    bool skipHeader;
};

// Utilities
namespace LoggingUtils { void setupLogging(const std::string &logFile, const char *execName); }
namespace ThreadUtils { void configureThreads(size_t threads); }
namespace FilesystemUtils {
    bool checkPathExists(const std::string &path, const std::string &description);
    bool ensureDirectoryExists(const std::string &path, bool allowOverwrite);
}
namespace ValidationUtils { bool validateIndexingFlags(bool useBulk, bool suffixSize, long double resourceLimit); }
namespace Metadata {
    struct GenomeManifest { std::unordered_map<size_t,GenomeToken> ids; };
    struct AnnotationManifest { std::unordered_map<uint64_t,Gff3Record> ids; };
    GenomeManifest loadGenomeManifest(const std::string &file);
    AnnotationManifest loadAnnotationManifest(const std::string &file);
    using GffRequestSet = std::unordered_map<Gff3Record,uint64_t>;
    using AccessionManifest = std::unordered_map<std::string,uint64_t>;
    AccessionManifest collectAccessions(const std::string &file, uint8_t);
    std::unordered_map<std::string,GffRequestSet> retrieveAllSources(uint8_t baseIdWidth, const GenomeManifest &genomeIds, const AnnotationManifest &annotIds);
}
namespace PhylogenyUtils {
    struct Phylogeny {
        maki_tree::Tree tree;
        maki_tree::lca::LcaTable lcaIndex;
        Phylogeny(const std::string&);
        Phylogeny(maki_tree::Tree&&);
        template <typename ...Args> size_t lca(Args &&...args) const { return lcaIndex.lca(std::forward<Args>(args)...); }
    };
    Phylogeny loadAndReduce(const std::string &nwkFile, const Metadata::GenomeManifest &genomes);
    std::string getNearestNamedLca(const std::string &n1, const std::string &n2, const Phylogeny &phy);
    std::string getNearestNamedAncestor(size_t nodeId, const Phylogeny &phy);
    void colourToNodeId(std::vector<size_t>&, const dist::bitmask&, const Metadata::GenomeManifest&, const Phylogeny&);
}
namespace ClusterUtils {
    class Clustering {
        struct StrHash { // based on https://en.cppreference.com/w/cpp/container/unordered_set/find.html#Example
            using hash_type = std::hash<std::string_view>;
            using is_transparent = void;
            size_t operator()(const char* s) const;
            size_t operator()(std::string_view s) const;
            size_t operator()(const std::string& s) const;
        };
        uint8_t baseIdWidth;
        std::unordered_set<std::string,StrHash,std::equal_to<>> seeds; // collection of seed labels
        std::unordered_map<uint64_t,const std::string*> lookup; // lookup seed label for query
        std::unordered_map<const std::string*,const std::string*> reps; // lookup representative seed
    public:
        Clustering(uint8_t baseIdWidth);
        Clustering(Clustering&&) =delete;
        Clustering(const Clustering&) =delete;
        Clustering& operator=(Clustering&&) =delete;
        Clustering& operator=(const Clustering&) =delete;
        const std::string* seed(uint64_t qry) const; // seed that a query belongs to
        const std::string* rep(uint64_t qry) const; // cluster that a query belongs to
        void insert(uint64_t qry, std::string &&seed); // map a query to a seed
        void assignIfExist(std::string_view rep, std::string_view seed); // assign a seed to a cluster
        void insertGenomeSeeds(const Metadata::GenomeManifest&);
        void loadAllSeeds(const Metadata::GenomeManifest&, const Metadata::AnnotationManifest&);
        void loadSeeds(const std::unordered_set<size_t>&, const Metadata::GenomeManifest&, const Metadata::AnnotationManifest&);
        void retrieveRequests(const std::unordered_map<std::string, Metadata::GffRequestSet> &requests);
        void loadClusters(const std::string &);
    };
    std::string commonCluster(const Clustering&, uint64_t, uint64_t);
    using GenomeClusters = std::unordered_map<std::string,std::vector<uint64_t>>;
    GenomeClusters groupByGenomeId(const std::string&, const Metadata::AccessionManifest&, const Metadata::GenomeManifest&, const PhylogenyUtils::Phylogeny&);
}
namespace ParseUtils {
    template <typename Stream, typename Fn> void parseFromStream(Stream &fh, Fn f);
    template <typename Stream, typename Fn> void parseFromStream(Stream &fh, size_t skip, Fn f);
    template <class ...Args> void parseGz(const std::string &fpath, Args&&...args);
    template <class ...Args> void parseRaw(const std::string &fpath, Args&&...args);
    template <class ...Args> void parse(const std::string &fpath, Args&&...args);
    size_t tab(const std::string &s, size_t t);
    size_t tabFind(const std::string &s, const char *qry);
    std::pair<size_t,size_t> tabSelect(const std::string &s, size_t t);
    size_t tabSelectAsSize_t(const std::string &s, size_t t);
    std::string tabExtractString(const std::string &s, size_t t);
    std::unordered_map<size_t,size_t> parseCounts(const std::string &f);
    template <typename T> void unique(std::vector<T> &data);
    template <typename C> void parseNumbersToList(C&, const std::string &);
    namespace Lca1 {
        struct CompressItem {
            size_t itemValue, numCoding, numColours, cladeSize;
            std::string lcaName;
            std::vector<uint64_t> values;
            std::vector<const std::string*> clusters;
            size_t numClusters() const noexcept;
            size_t numGenomes() const noexcept;
            void parseColours(const std::string &line, const ClusterUtils::Clustering &lookup);
            void collectPhylogenyNodes(const dist::bitmask &msk, const Metadata::GenomeManifest &genomeIds, const PhylogenyUtils::Phylogeny &phy);
            std::ostream& print(std::ostream &os) const;
        };
    }
    dist::cooc::RankMap parseColourProfiles(const std::string &colourFile, size_t bitWidth, bool skipHeader);
}
std::ostream& operator<<(std::ostream& os, const ParseUtils::Lca1::CompressItem& obj);
namespace GraphBuilder {
    struct ArgumentInvalidationException : public std::exception {};
    struct IndexParams {
        std::string queryFile, outputFolder;
        uint8_t kmerSize;
        size_t threads;
        long double resourceLimit;
        bool useBulk;
        uint8_t suffixSize;
        size_t edgeBound;
        size_t parsedBytes;
        size_t graphBytes;
        size_t roundedResourceLimit;
        IndexParams(const std::string&, const std::string&, uint8_t, size_t, long double, bool, uint8_t);
        bool validateBytesLimit(); // calculates `roundedResourceLimit`, assuming `edgeBound`, `parsedBytes`, `graphBytes`
    };
    template <class Container> graph::WriteableGraph buildIndexInMemory(IndexParams &params, const Container &genomes, const graph::GraphDiskFileConfig &config);
    template <class Container> graph::BaseGraph buildFilterInMemory(IndexParams &params, const Container &xgenome);
    bool buildIndex(IndexParams &params);
    bool buildFilter(IndexParams &params);
}
namespace StatCalculator {
    struct SummariseColourPairData {
        static constexpr size_t capacity = 10000;
        std::array<stats::SamdistResult,capacity> records;
        size_t size, baseIdWidth, k;
        const std::unordered_map<size_t, size_t> &uniqKmers;
        const Metadata::GenomeManifest &genomeIds;
        const PhylogenyUtils::Phylogeny &phy;
        const ClusterUtils::Clustering &clusters;
        std::ostream &os;
        SummariseColourPairData(size_t baseIdWidth_, size_t k_, const std::unordered_map<size_t,size_t> &uniqKmers_, const Metadata::GenomeManifest &genomeIds_, const PhylogenyUtils::Phylogeny &phy_, const ClusterUtils::Clustering &clusters_, std::ostream &os_);
        SummariseColourPairData(const SummariseColourPairData&);
        ~SummariseColourPairData();
        void operator()(size_t gid1, size_t gid2, size_t fid1, size_t fid2, size_t shared);
        void operator()(size_t gid1, size_t gid2, const dist::pairwise::ColourSampler::KeyValuePair &kv);
        bool full() const;
        void flush();
    };
};
namespace ReadsUtils {
    std::vector<reads::DataFilePair> pairReads(const std::string &readsFile);
    std::vector<reads::ReadChunks> chunkReads(const reads::DataFilePair &fp, size_t minReadLength, size_t minFragLength, size_t threads);
}
namespace QueryUtils {
    dist::query::QuerySet parseQueries(const std::string &queryFile, const Metadata::GenomeManifest &genomeIds);
}

// template definitions
namespace ParseUtils {
    template <typename Stream, typename Fn>
    void parseFromStream(Stream &fh, Fn f) {
        std::string line;
        while (std::getline(fh, line)) f(line);
        fh.close();
    }

    template <typename Stream, typename Fn>
    void parseFromStream(Stream &fh, size_t skip, Fn f) {
        size_t skipped = 0;
        // skip lines
        std::string line;
        while ((skipped < skip) && std::getline(fh, line)) ++skipped;
        // read rest of file
        while (std::getline(fh, line)) f(line);
        fh.close();
    }

    template <class ...Args>
    void parseGz(const std::string &fpath, Args&&...args) {
        zstr::ifstream fh(fpath);
        if (fh.is_open()) parseFromStream(fh, std::forward<Args>(args)...);
        else LOG(ERROR) << "Error opening file " << fpath;
    }

    template <class ...Args>
    void parseRaw(const std::string &fpath, Args&&...args) {
        std::ifstream fh(fpath);
        if (fh.is_open()) parseFromStream(fh, std::forward<Args>(args)...);
        else LOG(ERROR) << "Error opening file " << fpath;
    }

    template <class ...Args>
    void parse(const std::string &fpath, Args&&...args) {
        if (fpath.ends_with(".gz")) parseGz(fpath, std::forward<Args>(args)...);
        else parseRaw(fpath, std::forward<Args>(args)...);
    }

    template <typename T>
    void unique(std::vector<T> &data) {
        std::sort(data.begin(), data.end());
        data.erase(std::unique(data.begin(), data.end()), data.end());
    }

    template <typename C>
    void parseNumbersToList(C &out, const std::string &s) {
        size_t l = 0, k, r = s.size(), x;
        while (l < r) {
            k = std::min(r, s.find(',', l));
            try { x = std::stoull(s.substr(l, k - l)); }
            catch (const std::invalid_argument& e) { throw std::runtime_error("Could not parse number " + s.substr(l, k - l)); }
            out.push_back(x);
            l = k + 1;
        }
    }
}

namespace GraphBuilder {
    template <class Container>
    graph::WriteableGraph buildIndexInMemory(IndexParams &params, const Container &data, const graph::GraphDiskFileConfig &config) {
        params.edgeBound = detail::calculateNumEdgesBound(data);
        params.parsedBytes = detail::bytesUsed(data);
        params.graphBytes = graph::WriteableGraph::maxRSS(params.edgeBound, params.threads, detail::calculateNumColourValuesBound(data, params.kmerSize), data.getFeatureWidth());
        if (!params.validateBytesLimit()) throw ArgumentInvalidationException{};
        graph::WriteableGraph garr = graph::initialiseEmptyGraph(data, params.kmerSize, params.edgeBound, config.buffers);
        if (params.useBulk) graph::makeGraphBulk(garr, data, params.threads);
        else if (params.suffixSize) graph::makeGraphBySuffix(garr, data, params.kmerSize, params.suffixSize, params.threads);
        else if (params.roundedResourceLimit) graph::makeGraphDetect(garr, data, params.kmerSize, params.threads, params.roundedResourceLimit, _mkimem_details::record_size(_mkimem_details::key_size(params.kmerSize), data.getFeatureWidth()));
        else graph::makeGraphBulk(garr, data, params.threads);
        return garr;
    }

    template <class Container>
    graph::BaseGraph buildFilterInMemory(IndexParams &params, const Container &data) {
        params.edgeBound = detail::calculateNumEdgesBound(data);
        params.parsedBytes = detail::bytesUsed(data);
        params.graphBytes = graph::BaseGraph::maxRSS(params.edgeBound);
        if (!params.validateBytesLimit()) throw ArgumentInvalidationException{};
        graph::BaseGraph garr = graph::initialiseEmptyGraph(params.kmerSize, params.edgeBound);
        if (params.useBulk) graph::makeGraphBulk(garr, data, params.threads);
        else if (params.suffixSize) graph::makeGraphBySuffix(garr, data, params.kmerSize, params.suffixSize, params.threads);
        else if (params.roundedResourceLimit) graph::makeGraphDetect(garr, data, params.kmerSize, params.threads, params.roundedResourceLimit, _mkimem_details::record_size(_mkimem_details::key_size(params.kmerSize), 1));
        else graph::makeGraphBulk(garr, data, params.threads);
        return garr;
    }
}

namespace Metadata {
    template<typename Container>
    std::unordered_map<std::string,GffRequestSet> retrieveSources(const Container &ids, uint8_t baseIdWidth, const GenomeManifest &genomeIds, const AnnotationManifest &annotIds) {
        std::unordered_map<std::string,GffRequestSet> requests;
        dist::bitmask msk(baseIdWidth);
        for (const uint64_t &id : ids) {
            if (id == msk(id)) continue; // only searching for annotation IDs, not genome IDs
            auto grec = genomeIds.ids.find(msk(id));
            if (grec == genomeIds.ids.cend()) throw std::runtime_error("Unrecognised genome ID " + std::to_string(msk(id)));
            auto arec = annotIds.ids.find(id);
            if (arec == annotIds.ids.cend()) throw std::runtime_error("Unrecognised annot ID " + std::to_string(id));
            requests[grec->second.source][arec->second] = id;
        }
        return requests;
    }
}
