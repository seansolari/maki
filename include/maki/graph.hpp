#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <stack>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_pipeline.h>
#include <sdsl/vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <highfive/highfive.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>
#include <zstd.h> // presumes zstd library is installed
#include <maki/maki.h>
#include <maki/bitvectors.hpp>
#include <maki/fasta.hpp>
#include <maki/reads.hpp>
#include <maki/matrix.hpp>
#include <maki/suffix.hpp>
#include <maki/mers.hpp>
#include <maki/radixSort.hpp>
#include <maki/terminals.hpp>
#include <maki/timing.hpp>
#include <maki/utils.hpp>
#include <maki/vectorbuffer.hpp>
#include <maki/flush.hpp>

namespace fs = std::filesystem;
namespace pv = bv::threaded;

using oneapi::tbb::filter_mode::parallel;
using oneapi::tbb::filter_mode::serial_in_order;

namespace detail {
    size_t calculateNumEdgesBound(const Dna4GenomeVector &genomes);
    size_t calculateNumEdgesBound(const ChunkedDna4Genome &genome);
    size_t calculateNumEdgesBound(const std::vector<reads::ReadChunks> &reads);
    size_t calculateNumColourValuesBound(const Dna4GenomeVector &genomes, uint8_t k);
    size_t bytesUsed(const Dna4GenomeVector &vec);
    size_t bytesUsed(const ChunkedDna4Genome &genome);
    size_t bytesUsed(const std::vector<reads::ReadChunks> &reads);

} // namespace detail

namespace graph {
    // forward class defs
    class ColourBuffer;
    class ColourBufferRegistry;
    class BaseGraph;
    class WriteableGraph;
    struct TraversableGraph;

    struct GraphDiskFileConfig {
        fs::path base, edges, wplus, buffers, meta, genomesids, annotids;
        GraphDiskFileConfig(std::string _base);
        bool validate() const;
    };

    struct StaticGraphDiskFileConfig {
        fs::path base, last, lastRank, lastSelect, W, buffers, meta, genomesids, annotids, tempW, tempDir;
        StaticGraphDiskFileConfig(std::string _base);
        bool validate() const;
    };

    struct ColourBufferIterator {
        using difference_type = typename sdsl::bit_vector::difference_type;

        typename sdsl::int_vector<0>::iterator colour;
        typename sdsl::bit_vector::iterator    boundary;

        bool operator==(const ColourBufferIterator &rhs) const;
        bool operator!=(const ColourBufferIterator &rhs) const;

        difference_type operator-(const ColourBufferIterator &rhs) const;

        ColourBufferIterator &operator+=(size_t dx);
        ColourBufferIterator &operator++();
    };

    struct ColourBufferConstIterator {
        using difference_type = typename sdsl::bit_vector::difference_type;

        typename sdsl::int_vector<0>::const_iterator colour;
        typename sdsl::bit_vector::const_iterator    boundary;

        bool operator==(const ColourBufferConstIterator &rhs) const;
        bool operator!=(const ColourBufferConstIterator &rhs) const;

        difference_type operator-(const ColourBufferConstIterator &rhs) const;

        ColourBufferConstIterator& operator+=(size_t dx);
        ColourBufferConstIterator& operator++();
        ColourBufferConstIterator operator+(size_t dx) const;
    };

    struct ColourBuffer {
        sdsl::int_vector<0> colours;
        sdsl::bit_vector    boundaries;

        void reserve(size_t size);
        void resize(size_t newSize);
        void clear();
        void bit_compress();
        void set_width(size_t width);
        size_t size() const noexcept;
        bool empty() const;

        ColourBuffer() =default;
        ColourBuffer(uint8_t colourWidth, size_t numEntries);
        ColourBuffer(uint8_t colourWidth);

        ColourBufferIterator begin();
        ColourBufferIterator end();
        ColourBufferConstIterator cbegin() const;
        ColourBufferConstIterator cend() const;
    };

    // iteration

    struct _BaseCursor {
        _BaseCursor(size_t i_);
        _BaseCursor() : _BaseCursor(0) {}

        void resetBaseCursor(size_t i_);

        // tracking members
        int64_t                currentBlock;
        std::array<int64_t, 5> lastBlock;
        std::array<size_t,  5> F,
                               C;
    protected:
        size_t                   _i;
        std::vector<BufferValue> _tmpValues;

    public:
        void operator++();
        void operator+=(size_t di_);
        void operator=(size_t rhs_);
        size_t operator-(const _BaseCursor &other) const noexcept;
        size_t operator-(size_t other) const noexcept;
        size_t position() const noexcept;
        operator size_t() const noexcept;
        operator int64_t() const;
        void setToCurrent() noexcept;
        void setToCurrent(uint8_t e) noexcept;
        bool notCurrentBlock(uint8_t e) const noexcept;
    };

    struct BaseWriteCursor : public _BaseCursor
    {
        friend class BaseGraph;
        
    private:
        pv::SpinRegionManager::Accessor locks;
        sdsl::int_vector<4>             &edges;
        sdsl::bit_vector                &wplus;
        uint8_t                         __k;

    protected:
        auto sortTempValues();

    public:
        BaseWriteCursor(BaseGraph &graph, size_t i_);
        BaseWriteCursor(BaseGraph &graph);
        BaseWriteCursor(BaseWriteCursor const &rhs, size_t i);

        BaseWriteCursor moveTo(size_t i_) const;
        void emplace(BufferValue &&v);
        void writeNode(uint8_t msb);
        uint8_t k() const noexcept;
    };

    struct GraphWriteCursor : public _BaseCursor
    {
        friend class WriteableGraph;

        GraphWriteCursor(WriteableGraph &graph, size_t i_, uint8_t colWidth_);
        GraphWriteCursor(WriteableGraph &graph, size_t i_, uint8_t colWidth_, size_t, ZstdCompressor*);
        GraphWriteCursor(WriteableGraph &graph);
        GraphWriteCursor(GraphWriteCursor const&, size_t, size_t, ZstdCompressor*);
        GraphWriteCursor(GraphWriteCursor const&, size_t);
        ~GraphWriteCursor();

    private:
        pv::SpinRegionManager::Accessor locks;
        sdsl::int_vector<4>             &edges;
        sdsl::bit_vector                &wplus;
        ColourBufferRegistry            &colourRegistry;
        uint8_t                         __k;
        ColourBuffer                    colours;
        ZstdCompressor                  *zcx;
    protected:
        uint8_t                         colourBitWidth;

    private:
        static constexpr size_t MAX_BUFFER_SIZE = GRAPH_FLUSH_MAX_TEMP_BYTES / 8;

    protected:
        auto sortTempValues();

    public:
        GraphWriteCursor moveTo(size_t i_, size_t reserve, ZstdCompressor *cx);
        GraphWriteCursor moveTo(size_t i_);
        void setCompressor(ZstdCompressor *zcx_);
        void emplace(BufferValue &&v);
        void compressColourBuffer();
        void resetColourBuffer();
        void writeNode(uint8_t msb);
        void flushColours();
        size_t getColourPosition() const;
        uint8_t k() const noexcept;
        sdsl::int_vector<0> const& getColours() const noexcept;
        sdsl::bit_vector const& getBoundaries() const noexcept;
    };

    template <typename Graph_>
    struct graph_traits
    {
        typename Graph_::cursor_type cursor_type;
    };

    template <typename Graph_>
    using graph_cursor_t = graph_traits<Graph_>::cursor_type;

    // bulk fill

    void fillGraphBulk(GraphWriteCursor &crs, Dna4GenomeVector const               &genomes_, uint32_t threads_);
    void fillGraphBulk(BaseWriteCursor  &crs, ChunkedDna4Genome const              &genome_,  uint32_t threads_);
    void fillGraphBulk(BaseWriteCursor  &crs, std::vector<reads::ReadChunks> const &reads_,   uint32_t threads_);

    // suffix-wise fill

    void fillGraphBySuffix(GraphWriteCursor &crs, Dna4GenomeVector const               &genomes_, uint8_t s_, uint32_t threads_);
    void fillGraphBySuffix(BaseWriteCursor  &crs, ChunkedDna4Genome const              &genome_,  uint8_t s_, uint32_t threads_);
    void fillGraphBySuffix(BaseWriteCursor  &crs, std::vector<reads::ReadChunks> const &reads_,   uint8_t s_, uint32_t threads_);

    // construct API

    WriteableGraph  initialiseEmptyGraph(const Dna4GenomeVector &genomes, uint8_t k_, size_t edges_, fs::path bufferPath_);    
    BaseGraph       initialiseEmptyGraph(uint8_t k_, size_t edges_);

    template <typename Graph, typename Cursor>
    void finaliseGraphShape(Graph &graph_, Cursor &p) {
        graph_.resize(p);
        graph_.F[0] = 0u;
        graph_.F[1] = p.F[0];
        graph_.F[2] = p.F[0] + p.F[1];
        graph_.F[3] = p.F[0] + p.F[1] + p.F[2];
        graph_.F[4] = p.F[0] + p.F[1] + p.F[2] + p.F[3];
        graph_.C[0] = 0u;
        graph_.C[1] = p.C[0];
        graph_.C[2] = p.C[0] + p.C[1];
        graph_.C[3] = p.C[0] + p.C[1] + p.C[2];
        graph_.C[4] = p.C[0] + p.C[1] + p.C[2] + p.C[3];
    }

    template <typename Graph, typename Container>
    void makeGraphBulk(Graph &graph_, Container const &data_, uint32_t threads_) {
        LOG(INFO) << "using bulk k-mer sorting strategy";
        auto crs = graph_.from(0);
        fillGraphBulk(crs, data_, threads_);
        finaliseGraphShape(graph_, crs);
    }

    template <typename Graph, typename Container>
    void makeGraphBySuffix(Graph &graph_, Container const &data_, uint8_t k_, uint8_t s_, uint32_t threads_) {
        if ((int)k_ - (int)s_ <= 1) {
            LOG(ERROR)
                << "required suffix size s=" << (int)s_
                << " is within 1 of k-mer size k=" << (int)k_
                << ",\n\tthis case is not implemented, as it suggests you do not have enough compute";
            exit(1);
        }

        LOG(INFO) << "using suffix-wise k-mer sorting strategy, with suffix size s=" << (int)s_;
        auto crs = graph_.from(0);
        fillGraphBySuffix(crs, data_, s_, threads_);
        finaliseGraphShape(graph_, crs);
    }

    template <typename Graph, typename Container>
    void makeGraphDetect(Graph &graph_, Container const &data_, uint8_t k_, uint32_t threads_, size_t resourceLimitBytes, size_t recordSize) {
        size_t indexingBytes = 2 * graph_.edge_count() * recordSize;

        if (indexingBytes < resourceLimitBytes)
            return makeGraphBulk(graph_, data_, threads_);
        else
        {
            uint8_t s = static_cast<uint8_t>(std::ceil(
                (
                    1.0l
                        + std::log2(graph_.edge_count())
                        + std::log2(recordSize)
                        - std::log2(resourceLimitBytes)
                ) / 2.0l
            ));

            if ((int)k_ - (int)s <= 1) {
                LOG(FATAL)
                    << "required suffix size s=" << (int)s
                    << " is within 1 of k-mer size k=" << (int)k_
                    << ",\n\tthis case is not implemented, as it suggests you do not have enough compute";
            }
            return makeGraphBySuffix(graph_, data_, k_, s, threads_);
        }
    }

    // graph IO

    void saveToDisk(const WriteableGraph &, std::string);
    void saveToDisk(const TraversableGraph &, std::string);
    void loadFromDisk(std::string, WriteableGraph &);
    void loadFromDisk(GraphDiskFileConfig const &, WriteableGraph &);
    void loadFromDisk(std::string, TraversableGraph &);
    void loadFromDisk(StaticGraphDiskFileConfig const &, TraversableGraph &);

    template <typename Graph>
    void saveSmallBuffers(fs::path const &file, Graph &g) {
        std::ofstream os(file, std::ios::binary);
        cereal::BinaryOutputArchive oarchive(os);
        oarchive( g );
    }

    template <typename Graph>
    void loadSmallBuffers(fs::path const &file, Graph &out) {
        std::ifstream is(file, std::ios::binary);
        cereal::BinaryInputArchive iarchive(is);
        iarchive( out );
    }

    // Graph Definition

    struct BufferIdPair {
        size_t  edgePosition,
                colourBufferId,
                boundaryBufferId,
                bufferLength;
        uint8_t colourWidth;

        BufferIdPair(size_t ePos_, size_t colId_, size_t bndId_, size_t len_, uint8_t w_);
        BufferIdPair() =default;
        BufferIdPair(const BufferIdPair &) =default;
        BufferIdPair(BufferIdPair &&) =default;
        BufferIdPair &operator=(const BufferIdPair &) =default;
        BufferIdPair &operator=(BufferIdPair &&) =default;
        bool operator==(const BufferIdPair &) const =default;

        std::string getColourBufferId() const;
        std::string getBoundaryBufferId() const;

        template <class Archive>
        void serialize(Archive &ar)
        { ar(edgePosition, colourBufferId, boundaryBufferId, bufferLength, colourWidth); }

        friend std::ostream &operator<<(std::ostream &, const BufferIdPair &obj);
        friend void PrintTo(const BufferIdPair &obj, std::ostream *os);
    };

    extern std::mutex registerLock;

    class ColourBufferRegistry {
    protected:
        std::string                     h5file;
        std::unique_ptr<HighFive::File> h5p;
        std::map<size_t, BufferIdPair>  buffers;

    public:
        static constexpr struct temp_flag_t{} temp_flag{};

        // empty registry for delayed initialisation
        ColourBufferRegistry() {}

        // create a new buffer
        ColourBufferRegistry(std::string_view h5file_);

        // temp-backed buffer
        ColourBufferRegistry(temp_flag_t)
            : ColourBufferRegistry((std::filesystem::temp_directory_path() / "__temp_ColourBufferRegistry.h5").string())
        {}

        template <class Archive>
        void save(Archive &ar) const
        {
            ar(h5file, buffers);
        }

        template <class Archive>
        void load(Archive &ar)
        {
            ar(h5file, buffers);

            if (!h5file.empty())
            {
                try {
                    h5p = std::make_unique<HighFive::File>(h5file, HighFive::File::ReadOnly);
                } catch(const HighFive::FileException &e) {
                    LOG(ERROR) << e.what() << '\n';
                }
            }
        }

    private:
        void saveBuffer(const std::string &bufferId, const CompressedIntVector &data);
        void readBuffer(const std::string &bufferId, CompressedIntVector &vec) const;
        BufferIdPair registerBuffer(size_t edgeIndex, size_t numElements, uint8_t elementWidth);

    public:
        void saveColours(size_t edgeIndex, ColourBuffer const &data, ZstdCompressor &zstd);
        void readColours(ColourBuffer &out, const BufferIdPair &ids, ZstdDecompressor &zstd) const;
        ColourBuffer readColours(const BufferIdPair &ids, ZstdDecompressor &zstd) const;
        void readBoundaries(sdsl::bit_vector &out, const BufferIdPair &ids, ZstdDecompressor &zstd) const;
        friend std::ostream &operator<<(std::ostream &, const ColourBufferRegistry &obj);

        friend void PrintTo(const ColourBufferRegistry &obj, std::ostream *os)
        {
            (*os) << obj;
            return;
        }

        std::map<size_t, BufferIdPair>::iterator begin();
        std::map<size_t, BufferIdPair>::iterator end();
        std::map<size_t, BufferIdPair>::const_iterator begin() const;
        std::map<size_t, BufferIdPair>::const_iterator end() const;
        std::map<size_t, BufferIdPair>::const_iterator cbegin() const;
        std::map<size_t, BufferIdPair>::const_iterator cend() const;
        std::reverse_iterator<std::map<size_t, graph::BufferIdPair>::const_iterator> rbegin() const;
        std::reverse_iterator<std::map<size_t, graph::BufferIdPair>::const_iterator> rend() const;
        std::map<size_t, BufferIdPair>::const_iterator getKey(size_t id) const;

        uint8_t getMaxColourWidth() const;
        size_t getTotalSize() const;
        ColourBuffer combineAllBuffers() const;
        size_t numBuffers() const noexcept;

        std::string_view filePath() const noexcept;

        bool operator==(const ColourBufferRegistry &) const =default;
    };

    class BaseGraph {
        friend struct BaseWriteCursor;

    public:
        using cursor_type = BaseWriteCursor;

        // bound on space usage of graph, in bytes
        inline static constexpr size_t maxRSS(size_t edges)  {
            return ((5ul /*bits-per-edge in graph*/ * edges) + 7ul) / 8;
        }

        static constexpr size_t LOCK_SAMPLE_RATE = 8 * 1024;

    public:
        BaseGraph() = default;

        BaseGraph(size_t _k,
                  uint64_t _numGenomes,
                  size_t _reserve);

        BaseGraph(uint8_t k_,
                  uint64_t _numGenomes,
                  sdsl::int_vector<4> &&edges_,
                  sdsl::bit_vector &&wplus_,
                  std::array<size_t,5> F_,
                  std::array<size_t,5> C_);

    protected:
        uint8_t               k;
        uint64_t              numGenomes;
        pv::SpinRegionManager locks;
        sdsl::int_vector<4>   edges;
        sdsl::bit_vector      wplus;

    public:
        std::array<size_t, 5> F = {0, 0, 0, 0, 0},
                              C = {0, 0, 0, 0, 0};

        template <class Archive>
        void serialize(Archive &ar)
        { ar(k, numGenomes, F, C); }

    public:
        inline uint8_t kmer_size() const noexcept { return k; }

        inline size_t edges_size() const noexcept { return edges.size(); }

        inline size_t wplus_size() const noexcept { return wplus.size(); }

        inline size_t edge_count() const { return edges.size(); }

        inline sdsl::int_vector<4> const& edges_ref() const noexcept
        { return edges; }

        inline sdsl::bit_vector const& wplus_ref() const noexcept
        { return wplus; }

        inline void resize(size_t _new_size)
        {
            edges.resize(_new_size);
            wplus.resize(_new_size);
        }

        inline cursor_type from(size_t pos) { return BaseWriteCursor(*this, pos); }

        friend void makeStaticGraphOnDisk(BaseGraph &g, fs::path const &databaseDir);

        friend void makeStaticGraphOnDisk(WriteableGraph &g, fs::path const &databaseDir);

        friend struct TraversableGraph;
    };

    // represent (ordered) alphabet [ $, A, C, G, T ] in BITS_PER_EDGE bits
    class WriteableGraph : public BaseGraph
    {
        friend struct GraphWriteCursor;

    public:
        using cursor_type = GraphWriteCursor;

        // bound on space usage of graph, in bytes
        static constexpr size_t maxRSS(size_t edges, size_t threads, size_t cedges, size_t cwidth)  {
            size_t
                strBytes = (((5ul /*bits-per-edge in graph*/ * edges) + 7ul) / 8),
                colBytes = std::min(
                    threads * (size_t)GRAPH_FLUSH_MAX_TEMP_BYTES,
                    ((cwidth * cedges + 7ul) / 8));
            return strBytes + colBytes;
        }

    public:
        WriteableGraph() = default;

        WriteableGraph(size_t _k,
                       uint64_t _maxColour,
                       uint64_t _numGenomes,
                       uint8_t baseColourWidth_,
                       size_t _reserve);

        WriteableGraph(size_t _k,
                       uint64_t _maxColour,
                       uint64_t _numGenomes,
                       uint8_t baseColourWidth_,
                       size_t _reserve,
                       std::string _bufferPath);

        WriteableGraph(size_t _k,
                       uint64_t _maxColour,
                       uint64_t _numGenomes,
                       size_t _reserve)
            : WriteableGraph(_k, _maxColour, _numGenomes, required_bits(_maxColour), _reserve)
        {
        }

        WriteableGraph(size_t _k,
                       uint64_t _maxColour,
                       uint64_t _numGenomes,
                       size_t _reserve,
                       std::string _bufferPath)
            : WriteableGraph(_k, _maxColour, _numGenomes, required_bits(_maxColour), _reserve, _bufferPath)
        {
        }

        WriteableGraph(uint8_t k_,
                       uint64_t maxColour_,
                       uint64_t _numGenomes,
                       uint8_t baseColourWidth_, // number of LSB in colour IDs corresponding to sequence ID
                       sdsl::int_vector<4> &&edges_,
                       sdsl::bit_vector &&wplus_,
                       ColourBufferRegistry &&colourReg_,
                       std::array<size_t,5> F_,
                       std::array<size_t,5> C_)
            : BaseGraph(k_, _numGenomes, std::move(edges_), std::move(wplus_), F_, C_)
            , maxColour(maxColour_)
            , colourBitWidth(required_bits(maxColour_))
            , baseIdBitWidth(baseColourWidth_)
            , colourBuffers(std::move(colourReg_))
        {
        }

        WriteableGraph(uint8_t k_,
                       uint64_t colours_,
                       uint64_t _numGenomes,
                       sdsl::int_vector<4> &&edges_,
                       sdsl::bit_vector &&wplus_,
                       ColourBufferRegistry &&colourReg_,
                       std::array<size_t,5> F_,
                       std::array<size_t,5> C_)
            : BaseGraph(k_, _numGenomes, std::move(edges_), std::move(wplus_), F_, C_)
            , maxColour(colours_)
            , colourBitWidth(required_bits(colours_))
            , baseIdBitWidth(colourBitWidth)
            , colourBuffers(std::move(colourReg_))
        {
        }

    protected:
        uint64_t maxColour;
        uint8_t  colourBitWidth,
                 baseIdBitWidth;

    public:
        ColourBufferRegistry colourBuffers;

    public:
        template <class Archive>
        void serialize(Archive &ar)
        {
            BaseGraph::serialize(ar);
            ar(maxColour, colourBitWidth, baseIdBitWidth, colourBuffers);
        }

    public:
        inline uint64_t getMaxColour() const noexcept { return maxColour; }
        inline uint8_t getColourBitWidth() const noexcept { return colourBitWidth; }
        inline uint8_t getBaseColourBitWidth() const noexcept { return baseIdBitWidth; }
        inline cursor_type from(size_t pos) { return GraphWriteCursor(*this, pos, colourBitWidth); }
        inline cursor_type from(size_t pos, size_t reserve, ZstdCompressor *cx) { return GraphWriteCursor(*this, pos, colourBitWidth, reserve, cx); }

        friend void makeStaticGraphOnDisk(BaseGraph &g, fs::path const &databaseDir);
        friend void makeStaticGraphOnDisk(WriteableGraph &g, fs::path const &databaseDir);
        friend void saveToDisk(WriteableGraph const&, std::string);
        friend void loadFromDisk(GraphDiskFileConfig const &, WriteableGraph &);
        friend struct TraversableGraph;
    };

    namespace detail
    {
        
        template <typename Cursor_, typename Unwind_, typename Kmers_>
        void suffixwiseFill(Cursor_ &crs_,
                            std::vector<range_t<Unwind_>> const &data,
                            Unwind_ Apply,
                            std::vector<suffix::SuffixTable> const &blocks,
                            Kmers_ &kmers,
                            Kmers_ &temp,
                            terminals::TerminalRange &tmls,
                            suffix::SmallRollingNuclSeq suffix_,
                            uint8_t s_,
                            uint32_t threads);

        template <typename Unwind_, typename Kmers_>
        void extractAndSortSuffix(std::vector<range_t<Unwind_>> const &data,
                                  Unwind_ Apply,
                                  std::vector<suffix::SuffixTable> const &counts,
                                  suffix::SmallRollingNuclSeq suffix_,
                                  Kmers_ &kmers,
                                  Kmers_ &temp,
                                  uint32_t threads);

        template <typename Cursor_>
        void insertPartialSuffix(Cursor_ &crs_,
                                 terminals::TerminalRange &tmls,
                                 suffix::SmallRollingNuclSeq suffix_);

    } // namespace detail
    
    extern std::mutex lockForBV;

    struct TraversableGraph
    {
        friend void makeStaticGraphOnDisk(BaseGraph &g, fs::path const &databaseDir);
        friend void makeStaticGraphOnDisk(WriteableGraph &g, fs::path const &databaseDir);

        TraversableGraph() = default;
        explicit TraversableGraph(BaseGraph &&g);
        explicit TraversableGraph(BaseGraph &&g, fs::path const &tempFile);
        explicit TraversableGraph(WriteableGraph &&g);
        explicit TraversableGraph(WriteableGraph &&g, fs::path const &tempFile);

        TraversableGraph(uint8_t k,
                         uint64_t maxColour,
                         uint64_t numGenomes,
                         uint8_t colourBitWidth,
                         uint8_t baseIdBitWidth,
                         ColourBufferRegistry &&colourBuffers,
                         std::array<size_t,5> const &F,
                         std::array<size_t,5> const &C);

        // member functions

        template <class Archive>
        void serialize(Archive &ar)
        {
            ar(k, maxColour, numGenomes, colourBitWidth, baseIdBitWidth, colourBuffers, F, C);
        }

        // member access

        inline uint64_t getMaxColour() const noexcept { return maxColour; }
        inline uint64_t getNumGenomes() const noexcept { return numGenomes; }
        inline uint8_t getColourBitWidth() const noexcept { return colourBitWidth; }
        inline uint8_t getBaseColourBitWidth() const noexcept { return baseIdBitWidth; }

    public:
        struct IndexRange
        {
            IndexRange(size_t begin_, size_t end_) : begin(begin_), end(end_) {}

            size_t begin, end;

            inline size_t size() const noexcept { return end - begin; }
            inline bool empty() const noexcept { return begin == end; }
        };

        // traversal
    protected:
        uint8_t block(size_t i) const;

        // Return index of first occurrence of `1` at or before `i`.
        inline size_t lastPred(size_t i) const { return lastSelectSupport(lastRankSupport(i + 1)); }

        // Return index of first occurrence of `1` at or after `i`.
        inline size_t lastSucc(size_t i) const { return lastSelectSupport(lastRankSupport(i) + 1); }

        // Return index of first occurrence of `c` at or before `i`.
        inline size_t WPred(size_t i, uint8_t c) const { return W.select(W.rank(i + 1, c), c); }

        // Return index of first occurrence of `c` at or after `i`.
        inline size_t WSucc(size_t i, uint8_t c) const { return W.select(W.rank(i, c) + 1, c); }

    public:
        // step forward from `i`
        std::optional<size_t> fwd(size_t i) const;

        // step backward from `i`
        std::optional<size_t> bwd(size_t i) const;

        // returns edge value at position
        uint8_t edge(size_t i) const;

        // return k-mer represented by edge `i`
        Dna4Sequence kmer(size_t i) const;

        inline IndexRange getRoot() const { return IndexRange(0u, lastSucc(0u) + 1u); }

        // get Node/k-mer containing the edge at `i` in the form `[begin, end)`
        inline IndexRange getNode(size_t i) const
        {
            assert(i > 0u);
            return IndexRange(lastPred(i - 1) + 1u, lastSucc(i) + 1u);
        }

        // returns index of node by taking edge `e` from node at `i` (if the edge exists)
        std::optional<size_t> outgoing(size_t i, uint8_t e) const;

    protected:
        struct _IndexRangeDepth : public IndexRange {
            _IndexRangeDepth(size_t begin_, size_t end_, size_t depth = 0)
                : IndexRange(begin_, end_),
                  depth(depth)
            {
            }

            size_t depth;
        };

        using _IndexRangeDepth_Stack = std::stack<_IndexRangeDepth,
                                                  std::vector<_IndexRangeDepth>>;

        struct FillTerminalBitVector {
        public:
            FillTerminalBitVector(TraversableGraph const *g_, sdsl::bit_vector *bv_)
                : g(g_),
                  bv(bv_)
            {
            }

        public:
            template <typename Range>
            inline void operator()(Range &r) const
            {
                for (auto &&item : r)
                    operator()(std::move(item));
            }

            void operator()(_IndexRangeDepth &&b_) const;

        protected:
            void flushCache(std::vector<size_t> &cache) const;

        protected:
            TraversableGraph const *g;
            sdsl::bit_vector *bv;

            static const size_t maxCacheSize = 1000u;
        };

    public:
        // return bit-vector whose length is the number of nodes in the graph `m`,
        // and where an item is `1` if the corresponding node is a terminal
        sdsl::bit_vector getTerminals() const;

    public:
        class SuccinctTerminalIndex
        {
        public:
            SuccinctTerminalIndex(TraversableGraph const &g_);

        public:
            inline size_t size() const { return sdv.size(); }

            template <class... Args>
            inline auto rank(Args &&...args) const { return sdvRank.rank(std::forward<Args>(args)...); }

            template <class... Args>
            inline auto select(Args &&...args) const { return sdvSelect.select(std::forward<Args>(args)...); }

        protected:
            sdsl::sd_vector<> sdv;
            sdsl::sd_vector<>::rank_1_type sdvRank;
            sdsl::sd_vector<>::select_1_type sdvSelect;
        };

        class KmerRangeGenerator
        {
        public:
            KmerRangeGenerator(TraversableGraph const &g_);

            KmerRangeGenerator(KmerRangeGenerator const &g_, size_t begin, size_t end);

            std::optional<IndexRange> next();

            // `begin` and `end` specify range in terms of edges, not nodes
            inline KmerRangeGenerator subrange(size_t begin, size_t end) const { return KmerRangeGenerator(*this, begin, end); }

        protected:
            std::shared_ptr<SuccinctTerminalIndex> sdv;
            sdsl::rank_support_v5<1,1>    const    &lastRank;
            sdsl::select_support_mcl<1,1> const    &lastSelect;

            size_t i,           // current edge index
                   endI,        // one-past-end index of last colour
                   tI,          // current terminal index
                   maxTerminal; // terminal at which to stop
        };

        inline KmerRangeGenerator nonTerminalNodes() const { return KmerRangeGenerator(*this); }

    public:
        struct ColourBufferRange {
            ColourBufferRange(const ColourBuffer &buffer, size_t grainsize, size_t edgeEnd);
            ColourBufferRange(const ColourBuffer &buffer, size_t grainsize);
            ColourBufferRange(const ColourBufferRange &rhs) = default;
            ColourBufferRange(ColourBufferRange &lhs, oneapi::tbb::proportional_split p);
            ColourBufferRange(ColourBufferRange &lhs, oneapi::tbb::split);
        protected:
            ColourBufferRange(const ColourBufferRange &rhs, size_t begin, size_t end);
        public:
            ColourBufferRange subrange(size_t begin, size_t end) const;
            ColourBufferConstIterator begin() const;
            ColourBufferConstIterator end() const;
            bool empty() const;
            size_t size() const;
            size_t load() const;
            bool is_divisible() const;
            size_t edgeBegin() const;
            size_t edgeEnd() const;
            ColourBufferConstIterator selectEdge(size_t i) const;
        protected:
            ColourBufferConstIterator _p;
            std::shared_ptr<sdsl::rank_support_v5<1,1>>    rp;
            std::shared_ptr<sdsl::select_support_mcl<1,1>> sp;
            size_t _begin, _end, _grainsize, _edgeOffset;
        };

    public:
        // graph structure params

        uint8_t  k;
        uint64_t maxColour,
                 numGenomes;
        uint8_t  colourBitWidth,
                 baseIdBitWidth;

        // graph structure data

        sdsl::bit_vector               last;
        sdsl::rank_support_v5<1,1>     lastRankSupport;
        sdsl::select_support_mcl<1,1>  lastSelectSupport;
        sdsl::wt_int<>                 W;
        ColourBufferRegistry           colourBuffers;
        std::array<size_t, 5>          F,
                                       C;

    public:
        inline size_t node_count() const { return lastRankSupport(last.size()); }
        inline size_t edge_count() const { return W.size(); }
    };

    void makeLastBufferOnDisk(sdsl::bit_vector &wplus, StaticGraphDiskFileConfig const &cfg);
    void makeEdgesBufferOnDisk(sdsl::int_vector<4> &edges, StaticGraphDiskFileConfig const &cfg);
    void makeStaticGraphOnDisk(BaseGraph &g, fs::path const &databaseDir);
    void makeStaticGraphOnDisk(WriteableGraph &g, fs::path const &databaseDir);

    // template definitions ------------------------------------------------------

    namespace detail
    {

        template <typename Cursor_, typename Unwind_, typename Kmers_>
        void suffixwiseFill(Cursor_ &crs_,
                            std::vector<range_t<Unwind_>> const &data,
                            Unwind_ Apply,
                            std::vector<suffix::SuffixTable> const &blocks,
                            Kmers_ &kmers,
                            Kmers_ &temp,
                            terminals::TerminalRange &tmls,
                            suffix::SmallRollingNuclSeq suffix_,
                            uint8_t s_,
                            uint32_t threads)
        {
            if (suffix_.size() == s_)
            {
                LOG(INFO) << "inserting suffix " << suffix_.toString();

                extractAndSortSuffix(data, Apply, blocks, suffix_, kmers, temp, threads);
                terminals::TerminalRange tx = tmls.endsWith(suffix_);

                flush::flush(crs_, kmers, tx,
                             BaseGraph::LOCK_SAMPLE_RATE,
                             static_cast<uint8_t>(mers::dna4ToDna5(suffix_.msb())));
            } else {
                insertPartialSuffix(crs_, tmls, suffix_);

                /* A */ suffixwiseFill(crs_, data, Apply, blocks, kmers, temp, tmls,
                                       suffix_ + uint64_t(0b00), s_, threads);
                /* C */ suffixwiseFill(crs_, data, Apply, blocks, kmers, temp, tmls,
                                       suffix_ + uint64_t(0b01), s_, threads);
                /* G */ suffixwiseFill(crs_, data, Apply, blocks, kmers, temp, tmls,
                                       suffix_ + uint64_t(0b10), s_, threads);
                /* T */ suffixwiseFill(crs_, data, Apply, blocks, kmers, temp, tmls,
                                       suffix_ + uint64_t(0b11), s_, threads);
            }
        }

        template <typename Unwind_, typename Kmers_>
        void extractAndSortSuffix(std::vector<range_t<Unwind_>> const &data,
                                  Unwind_ Apply,
                                  std::vector<suffix::SuffixTable> const &counts,
                                  suffix::SmallRollingNuclSeq suffix_,
                                  Kmers_ &kmers,
                                  Kmers_ &temp,
                                  uint32_t threads)
        {
            // block start for each genome

            auto blocks = accumulatePartials(counts, [&](suffix::SuffixTable const &st)->size_t
                { return st[suffix_]; });
            size_t nRecs = blocks.back();

            LOG(INFO) << "inserting " << nRecs << " k-mers with this suffix";

            // extract and sort k-mers with this suffix

            kmers.resize(nRecs);
            temp.resize(nRecs);
            kmers.insertKmers(data, Apply, blocks, suffix_);
            if (nRecs > 1)
                kmers.sort(&temp, threads);
        }

        template <typename Cursor_>
        void insertPartialSuffix(Cursor_ &crs_,
                                 terminals::TerminalRange &tmls,
                                 suffix::SmallRollingNuclSeq suffix_)
        {
            auto vals = tmls.retrieve(suffix_);
            if (!vals.empty())
            {
                flush::flush(crs_, vals, vals.size());
            }
        }

    } // namespace detail

} // namespace graph
