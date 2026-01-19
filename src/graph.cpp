#include <maki/graph.hpp>

#include <iostream>
#include <sdsl/util.hpp>
#include <oneapi/tbb/parallel_reduce.h>

// De Bruijn graph ==========================================================================

namespace detail
{

    size_t calculateNumEdgesBound(const Dna4GenomeVector &genomes)
    {
        LOG(INFO) << "calculating bound on number of edges";

        // each nucleotide in the sequence contributes up to 1 edge, and also
        // need to account for terminal character in each sequence.
        return oneapi::tbb::parallel_reduce(
            oneapi::tbb::blocked_range(genomes.cbegin(), genomes.cend()),
            (size_t)0,
            [](const oneapi::tbb::blocked_range<Dna4GenomeVector::BaseGenomeVectorConstIter> &r, size_t total) -> size_t
            {
                for (const Dna4Genome &g : r)
                {
                    total += g.length() + g.numFragments();
                }
                return total;
            },
            std::plus{});
    }

    size_t calculateNumEdgesBound(const ChunkedDna4Genome &genome)
    {
        return genome.numEdges();
    }

    size_t calculateNumEdgesBound(const std::vector<reads::ReadChunks> &reads)
    {
        return oneapi::tbb::parallel_reduce(
            oneapi::tbb::blocked_range(reads.cbegin(), reads.cend()),
            (size_t)0,
            [](const oneapi::tbb::blocked_range<typename std::vector<reads::ReadChunks>::const_iterator> &r, size_t total) -> size_t
            {
                for (const reads::ReadChunks &rset : r)
                {
                    total += reads::detail::numEdges(rset);
                }
                return total;
            },
            std::plus{});
    }

    size_t calculateNumColourValuesBound(const Dna4GenomeVector &genomes, uint8_t k)
    {
        LOG(INFO) << "calculating bound on number of colour values";

        // also add starting terminals
        return oneapi::tbb::parallel_reduce(
            oneapi::tbb::blocked_range(genomes.cbegin(), genomes.cend()),
            (size_t)0,
            [k](const oneapi::tbb::blocked_range<Dna4GenomeVector::BaseGenomeVectorConstIter> &r, size_t total) -> size_t
            {
                for (const Dna4Genome &g : r)
                {
                    total += g.numKmers(k) + (k * g.numFragments());
                }
                return total;
            },
            std::plus{});
    }

    size_t bytesUsed(const Dna4GenomeVector &vec)
    {
        using GenomeRange = oneapi::tbb::blocked_range<Dna4GenomeVector::BaseGenomeVectorConstIter>;

        auto rangeSizeof = [](const GenomeRange &r, size_t currentSize) -> size_t {
            for (const Dna4Genome &g : r)
                currentSize += g.rss();
            return currentSize;
        };

        return oneapi::tbb::parallel_reduce(GenomeRange(vec.cbegin(), vec.cend()),
                                            (size_t)0,
                                            rangeSizeof,
                                            std::plus{});
    }

    size_t bytesUsed(const ChunkedDna4Genome &genome)
    {
        return genome.sequenceRSS() + genome.chunksRSS();
    }

    size_t bytesUsed(const std::vector<reads::ReadChunks> &reads)
    {
        using ReadChunksRange = oneapi::tbb::blocked_range<std::vector<reads::ReadChunks>::const_iterator>;
        return oneapi::tbb::parallel_reduce(
            ReadChunksRange(reads.cbegin(), reads.cend()),
            (size_t)0,
            [](ReadChunksRange const &r, size_t currentSize)->size_t {
                for (auto it = r.begin(); it != r.end(); ++it)
                    currentSize += it->rss();
                return currentSize;
            },
            std::plus{});
    }

} // namespace detail

namespace graph
{

    GraphDiskFileConfig::GraphDiskFileConfig(std::string _base)
        : base(_base)
        , edges(base / fs::path("edges.dat"))
        , wplus(base / fs::path("wplus.dat"))
        , buffers(base / fs::path("colours.h5"))
        , meta(base / fs::path("meta.dat"))
        , genomesids(base / fs::path("genome_ids.txt.gz"))
        , annotids(base / fs::path("annotation_ids.txt.gz"))
    {
    }

    bool GraphDiskFileConfig::validate() const
    { return fs::exists(edges) && fs::exists(wplus) && fs::exists(meta); }

    StaticGraphDiskFileConfig::StaticGraphDiskFileConfig(std::string _base)
        : base(_base)
        , last(base / "last.dat")
        , lastRank(base / "last_rank.dat")
        , lastSelect(base / "last_select.dat")
        , W(base / "W.dat")
        , buffers(base / "colours.h5")
        , meta(base / "meta.dat")
        , genomesids(base / "genome_ids.txt.gz")
        , annotids(base / "annotation_ids.txt.gz")
        , tempW(base / "tempW.dat")
        , tempDir(base / ".temp")
    {
    }

    bool StaticGraphDiskFileConfig::validate() const
    { return fs::exists(last) && fs::exists(W) && fs::exists(meta); }

    bool ColourBufferIterator::operator==(const ColourBufferIterator &rhs) const
    { return (colour == rhs.colour) && (boundary == rhs.boundary); }

    bool ColourBufferIterator::operator!=(const ColourBufferIterator &rhs) const
    { return (colour != rhs.colour) || (boundary != rhs.boundary); }

    ColourBufferIterator::difference_type ColourBufferIterator::operator-(const ColourBufferIterator &rhs) const
    { return boundary - rhs.boundary; }

    ColourBufferIterator& ColourBufferIterator::operator+=(size_t dx) {
        colour += dx;
        boundary += dx;
        return *this;
    }

    ColourBufferIterator& ColourBufferIterator::operator++() {
        ++colour;
        ++boundary;
        return *this;
    }

    bool ColourBufferConstIterator::operator==(const ColourBufferConstIterator &rhs) const
    { return (colour == rhs.colour) && (boundary == rhs.boundary); }

    bool ColourBufferConstIterator::operator!=(const ColourBufferConstIterator &rhs) const
    { return (colour != rhs.colour) || (boundary != rhs.boundary); }

    ColourBufferConstIterator::difference_type ColourBufferConstIterator::operator-(const ColourBufferConstIterator &rhs) const
    { return boundary - rhs.boundary; }

    ColourBufferConstIterator& ColourBufferConstIterator::operator+=(size_t dx) {
        colour += dx;
        boundary += dx;
        return *this;
    }

    ColourBufferConstIterator& ColourBufferConstIterator::operator++() {
        ++colour;
        ++boundary;
        return *this;
    }

    ColourBufferConstIterator ColourBufferConstIterator::operator+(size_t dx) const
    { return ColourBufferConstIterator(colour + dx, boundary + dx); }

    void ColourBuffer::reserve(size_t size) {
        colours.reserve(size);
        boundaries.reserve(size);
    }

    void ColourBuffer::resize(size_t newSize) {
        colours.resize(newSize);
        boundaries.resize(newSize);
    }

    void ColourBuffer::clear() {
        colours.clear();
        boundaries.clear();
    }

    void ColourBuffer::bit_compress()
    { sdsl::util::bit_compress(colours); }

    void ColourBuffer::set_width(size_t width)
    { colours.width(width); }

    size_t ColourBuffer::size() const noexcept
    { return colours.size(); }

    bool ColourBuffer::empty() const
    { return colours.empty() && boundaries.empty(); }

    ColourBuffer::ColourBuffer(uint8_t colourWidth, size_t numEntries)
        : colours(numEntries, 0, colourWidth)
        , boundaries(numEntries)
    {
    }

    ColourBuffer::ColourBuffer(uint8_t colourWidth)
        : colours()
        , boundaries()
    {
        set_width(colourWidth);
    }

    ColourBufferIterator ColourBuffer::begin()
    { return {colours.begin(), boundaries.begin()}; }

    ColourBufferIterator ColourBuffer::end()
    { return {colours.end(), boundaries.end()}; }

    ColourBufferConstIterator ColourBuffer::cbegin() const
    { return {colours.cbegin(), boundaries.cbegin()}; }

    ColourBufferConstIterator ColourBuffer::cend() const
    { return {colours.cend(), boundaries.cend()}; }

    void _BaseCursor::operator++()
    { ++_i; }

    void _BaseCursor::operator+=(size_t di_)
    { _i += di_; }

    void _BaseCursor::operator=(size_t rhs_)
    { _i = rhs_; }

    size_t _BaseCursor::operator-(const _BaseCursor &other) const noexcept
    { return _i - other._i; }

    size_t _BaseCursor::operator-(size_t other) const noexcept
    { return _i - other; }

    size_t _BaseCursor::position() const noexcept
    { return _i; }

    _BaseCursor::operator size_t() const noexcept
    { return _i; }

    _BaseCursor::operator int64_t() const
    { return static_cast<int64_t>(_i); }
    
    void _BaseCursor::setToCurrent() noexcept
    { currentBlock = static_cast<int64_t>(_i); }

    void _BaseCursor::setToCurrent(uint8_t e) noexcept
    { lastBlock[e] = currentBlock; }

    bool _BaseCursor::notCurrentBlock(uint8_t e) const noexcept
    { return lastBlock[e] != currentBlock; }

    BaseWriteCursor BaseWriteCursor::moveTo(size_t i_) const
    { return BaseWriteCursor(*this, i_); }

    void BaseWriteCursor::emplace(BufferValue &&v)
    { _tmpValues.emplace_back(mers::dna4ToDna5(v.edge())); }

    uint8_t BaseWriteCursor::k() const noexcept
    { return __k; }

    GraphWriteCursor GraphWriteCursor::moveTo(size_t i_, size_t reserve, ZstdCompressor *cx)
    { return GraphWriteCursor(*this, i_, reserve, cx); }

    GraphWriteCursor GraphWriteCursor::moveTo(size_t i_)
    { return GraphWriteCursor(*this, i_); }

    void GraphWriteCursor::compressColourBuffer() {
        colours.bit_compress();
    }

    void GraphWriteCursor::resetColourBuffer() {
        colours.clear();
        colours.set_width(colourBitWidth);
    }

    void GraphWriteCursor::setCompressor(ZstdCompressor *zcx_)
    { zcx = zcx_; }

    void GraphWriteCursor::emplace(BufferValue &&v)
    { _tmpValues.emplace_back(v.colour(), mers::dna4ToDna5(v.edge())); }

    size_t GraphWriteCursor::getColourPosition() const
    { return colours.size(); }

    uint8_t GraphWriteCursor::k() const noexcept
    { return __k; }

    BufferIdPair::BufferIdPair(size_t ePos_, size_t colId_, size_t bndId_, size_t len_, uint8_t w_)
        : edgePosition(ePos_),
            colourBufferId(colId_),
            boundaryBufferId(bndId_),
            bufferLength(len_),
            colourWidth(w_)
    {
    }

    std::string BufferIdPair::getColourBufferId() const
    { return std::to_string(colourBufferId); }

    std::string BufferIdPair::getBoundaryBufferId() const
    { return std::to_string(boundaryBufferId); }

    void PrintTo(const BufferIdPair &obj, std::ostream *os) {
        (*os) << obj;
        return;
    }

    std::map<size_t, BufferIdPair>::iterator ColourBufferRegistry::begin()
    { return buffers.begin(); }

    std::map<size_t, BufferIdPair>::iterator ColourBufferRegistry::end()
    { return buffers.end(); }

    std::map<size_t, BufferIdPair>::const_iterator ColourBufferRegistry::begin() const
    { return buffers.begin(); }

    std::map<size_t, BufferIdPair>::const_iterator ColourBufferRegistry::end() const
    { return buffers.end(); }

    std::map<size_t, BufferIdPair>::const_iterator ColourBufferRegistry::cbegin() const
    { return buffers.cbegin(); }

    std::map<size_t, BufferIdPair>::const_iterator ColourBufferRegistry::cend() const
    { return buffers.cend(); }

    std::reverse_iterator<std::map<size_t, graph::BufferIdPair>::const_iterator> ColourBufferRegistry::rbegin() const noexcept
    { return buffers.rbegin(); }

    std::reverse_iterator<std::map<size_t, graph::BufferIdPair>::const_iterator> ColourBufferRegistry::rend() const noexcept
    { return buffers.rend(); }

    std::string_view ColourBufferRegistry::filePath() const noexcept
    { return h5file; }

    std::map<size_t, BufferIdPair>::const_iterator ColourBufferRegistry::getKey(size_t id) const
    { return buffers.find(id); }

    size_t ColourBufferRegistry::numBuffers() const noexcept
    { return buffers.size(); }

    WriteableGraph initialiseEmptyGraph(const Dna4GenomeVector &genomes, uint8_t k_, size_t edges_, fs::path bufferPath_) {
        return WriteableGraph(k_, genomes.getMaxFeatureId(), genomes.numGenomes(), genomes.getBaseFeatureWidth(), edges_, bufferPath_);
    }

    BaseGraph initialiseEmptyGraph(uint8_t k_, size_t edges_) {
        return BaseGraph(k_, 1, edges_);
    }

    void saveToDisk(const WriteableGraph &graph, std::string databaseDir) {
        GraphDiskFileConfig config(databaseDir);

        // serialise metadata

        saveSmallBuffers(config.meta, graph);

        // serialise data buffers

        LOG(INFO) << "serialising edges to " << config.edges;
        fileutils::store_to_file(graph.edges, config.edges);

        LOG(INFO) << "serialising wplus to " << config.wplus;
        fileutils::store_to_file(graph.wplus, config.wplus);

        LOG(INFO) << "database successfully serialised to disk";
    }

    void saveToDisk(const TraversableGraph &graph, std::string databaseDir)
    {
        StaticGraphDiskFileConfig config(databaseDir);

        // serialise metadata

        saveSmallBuffers(config.meta, graph);

        // serialise data buffers

        LOG(INFO) << "serialising last to " << config.last;
        fileutils::store_to_file(graph.last, config.last);

        LOG(INFO) << "serialising lastRankSupport to " << config.lastRank;
        sdsl::store_to_file(graph.lastRankSupport, config.lastRank);

        LOG(INFO) << "serialising lastSelectSupport to " << config.lastSelect;
        sdsl::store_to_file(graph.lastSelectSupport, config.lastSelect);

        LOG(INFO) << "serialising W to " << config.W;
        sdsl::store_to_file(graph.W, config.W);
    }

    void loadFromDisk(std::string databaseDir, WriteableGraph &graph)
    {
        GraphDiskFileConfig config(databaseDir);
        loadFromDisk(config, graph);
    }

    void loadFromDisk(GraphDiskFileConfig const &config, WriteableGraph &graph)
    {
        // load metadata

        loadSmallBuffers(config.meta, graph);

        // load data buffers

        LOG(INFO) << "reading edges from " << config.edges;
        fileutils::load_from_file(graph.edges, config.edges);

        LOG(INFO) << "reading wplus from " << config.wplus;
        fileutils::load_from_file(graph.wplus, config.wplus);

        graph.locks.realloc(graph.edges.size(), WriteableGraph::LOCK_SAMPLE_RATE);

        LOG(INFO) << "graph successfully loaded";
    }

    void loadFromDisk(std::string databaseDir, TraversableGraph &graph)
    {
        StaticGraphDiskFileConfig config(databaseDir);
        loadFromDisk(config, graph);
    }

    void loadFromDisk(StaticGraphDiskFileConfig const &config, TraversableGraph &graph)
    {
        // load metadata

        loadSmallBuffers(config.meta, graph);

        // load data buffers

        LOG(INFO) << "loading last from " << config.last;
        fileutils::load_from_file(graph.last, config.last);

        LOG(INFO) << "loading lastRankSupport from " << config.lastRank;
        sdsl::load_from_file(graph.lastRankSupport, config.lastRank);
        graph.lastRankSupport.set_vector(&graph.last);

        LOG(INFO) << "loading lastSelectSupport from " << config.lastSelect;
        sdsl::load_from_file(graph.lastSelectSupport, config.lastSelect);
        graph.lastSelectSupport.set_vector(&graph.last);

        LOG(INFO) << "loading W from " << config.W;
        sdsl::load_from_file(graph.W, config.W);

        LOG(INFO) << "graph successfully loaded";
    }

    std::ostream &operator<<(std::ostream &os, const BufferIdPair &obj)
    {
        return os
               << "BufferIdPair("
               << obj.colourBufferId << ", "
               << obj.boundaryBufferId
               << ')';
    }

    BufferIdPair ColourBufferRegistry::registerBuffer(size_t edgeIndex, size_t numElements, uint8_t elementWidth)
    {
        // ensure we can create IDs using current schema
        assert(edgeIndex < std::numeric_limits<size_t>::max() - 1u);

        size_t
            colId = edgeIndex << 1,
            bndId = colId | 1u;

        // -- START THREAD LOCK --
        {
            const std::lock_guard<std::mutex> lock(registerLock);

            // save data

            auto result = buffers.try_emplace(/*key*/ edgeIndex,
                                              /*val*/ edgeIndex, colId, bndId, numElements, elementWidth);
            if (!result.second)
            {
                throw std::runtime_error("key already exists: " + std::to_string(edgeIndex));
            }

            return result.first->second;
        }
        // -- END THREAD LOCK --
    }

    void ColourBufferRegistry::saveBuffer(const std::string &bufferId, const CompressedIntVector &vec)
    {
        std::array<size_t, 1> dims = {vec.size()};

        // allocate datasets

        typedef std::remove_cv_t<std::remove_pointer_t<decltype(vec.data())>> Byte;

        // -- START THREAD LOCK --
        {
            const std::lock_guard<std::mutex> lock(registerLock);
            HighFive::DataSet dset = h5p->createDataSet<Byte>(bufferId, HighFive::DataSpace(dims));

            // dataset attributes

            dset.createAttribute("oelem", vec.info.numElements);
            dset.createAttribute("ow", vec.info.elementWidth);
            dset.createAttribute("osize", vec.info.numWords);
            dset.createAttribute("cxsize", vec.cxSize);

            dset.write_raw(vec.data());
        }
        // -- END THREAD LOCK --
    }

    void ColourBufferRegistry::readBuffer(const std::string &bufferId, CompressedIntVector &vec) const
    {
        // load dataset and attributes

        HighFive::DataSet dset = h5p->getDataSet(bufferId);

        dset.getAttribute("oelem").read(vec.info.numElements);
        dset.getAttribute("ow").read(vec.info.elementWidth);
        dset.getAttribute("osize").read(vec.info.numWords);
        dset.getAttribute("cxsize").read(vec.cxSize);

        vec.resize(vec.cxSize);

        // read compressed data

        dset.read_raw(vec.data());
    }

    void ColourBufferRegistry::saveColours(size_t edgeIndex, ColourBuffer const &data, ZstdCompressor &zstd)
    {
        BufferIdPair ids = registerBuffer(edgeIndex, data.colours.size(), data.colours.width());

        zstd.compress(data.colours);
        saveBuffer(ids.getColourBufferId(), zstd.rsl);

        zstd.compress(data.boundaries);
        saveBuffer(ids.getBoundaryBufferId(), zstd.rsl);
    }

    void ColourBufferRegistry::readColours(ColourBuffer &out, const BufferIdPair &ids, ZstdDecompressor &zstd) const
    {
        readBuffer(ids.getColourBufferId(), zstd.src);
        zstd.decompress(out.colours);

        readBuffer(ids.getBoundaryBufferId(), zstd.src);
        zstd.decompress(out.boundaries);
    }

    ColourBuffer ColourBufferRegistry::readColours(const BufferIdPair &ids, ZstdDecompressor &zstd) const
    {
        ColourBuffer data;
        readColours(data, ids, zstd);
        return data;
    }

    void ColourBufferRegistry::readBoundaries(sdsl::bit_vector &out, const BufferIdPair &ids, ZstdDecompressor &zstd) const
    {
        readBuffer(ids.getBoundaryBufferId(), zstd.src);
        zstd.decompress(out);
    }

    std::ostream &operator<<(std::ostream &os, const ColourBufferRegistry &obj)
    {
        os
            << "ColourBufferRegistry("
            << obj.h5file
            << "; { ";

        for (const auto &[k, v] : obj.buffers)
        {
            os
                << "{ "
                << k
                << ", "
                << v
                << " } ";
        }

        return os
               << "})";
    }

    size_t ColourBufferRegistry::getTotalSize() const
    {
        return std::accumulate(cbegin(), cend(), (size_t)0,
                               [](size_t total, auto const &kp)
                               {
                                   return total + kp.second.bufferLength;
                               });
    }

    uint8_t ColourBufferRegistry::getMaxColourWidth() const
    {
        const auto &buffer = *std::max_element(cbegin(), cend(), [](const auto &lhs, const auto &rhs)
            { return lhs.second.colourWidth < rhs.second.colourWidth; });
        return buffer.second.colourWidth;
    }

    ColourBuffer ColourBufferRegistry::combineAllBuffers() const
    {
        ZstdDecompressor zstd;

        // output buffer

        uint8_t colourWidth = getMaxColourWidth();
        size_t numElements = getTotalSize();
        ColourBuffer arr(colourWidth, numElements);

        // order of buffers

        std::vector<size_t> sortedKeys;
        sortedKeys.reserve(numBuffers());
        for (auto const &rec : buffers)
            sortedKeys.push_back(rec.first);
        std::sort(sortedKeys.begin(), sortedKeys.end());

        // write to output buffer

        auto colOut = arr.colours.begin();
        auto bndOut = arr.boundaries.begin();

        for (size_t const &key : sortedKeys)
        {
            ColourBuffer src = readColours(buffers.at(key), zstd);
            colOut = std::copy(src.colours.cbegin(), src.colours.cend(), colOut);
            bndOut = std::copy(src.boundaries.cbegin(), src.boundaries.cend(), bndOut);
        }

        return arr;
    }

    std::mutex registerLock;

    _BaseCursor::_BaseCursor(size_t i_)
        : currentBlock(i_),
          lastBlock({ -1, -1, -1, -1, -1 }),
          F({ 0, 0, 0, 0, 0 }),
          C({ 0, 0, 0, 0, 0 }),
          _i(i_),
          _tmpValues()
    {
    }

    void _BaseCursor::resetBaseCursor(size_t i_) {
        currentBlock = i_;
        lastBlock = { -1, -1, -1, -1, -1 };
        F = { 0, 0, 0, 0, 0 };
        C = { 0, 0, 0, 0, 0 };
        _i = i_;
        _tmpValues.clear();
    }

    BaseWriteCursor::BaseWriteCursor(BaseGraph &graph, size_t i_)
        : _BaseCursor(i_)
        , locks(graph.locks)
        , edges(graph.edges)
        , wplus(graph.wplus)
        , __k(graph.kmer_size())
    {
    }

    BaseWriteCursor::BaseWriteCursor(BaseGraph &graph)
        : _BaseCursor()
        , locks(graph.locks)
        , edges(graph.edges)
        , wplus(graph.wplus)
        , __k(graph.kmer_size())
    {
    }

    BaseWriteCursor::BaseWriteCursor(BaseWriteCursor const &rhs, size_t i)
        : _BaseCursor(i)
        , locks(rhs.locks, pv::SpinRegionManager::Accessor::shallow_copy_tag)
        , edges(rhs.edges)
        , wplus(rhs.wplus)
        , __k(rhs.__k)
    {
    }

    auto BaseWriteCursor::sortTempValues()
    {
        std::sort(_tmpValues.begin(),
                  _tmpValues.end(),
                  [](const BufferValue &lhs, const BufferValue &rhs) {
                      return lhs.edge() < rhs.edge();
                  });
        return std::unique(_tmpValues.begin(), _tmpValues.end());
    }

    void BaseWriteCursor::writeNode(uint8_t msb)
    {
        auto tmpEnd = sortTempValues();
        auto it = _tmpValues.cbegin();

        // if there are valid edges other than `$`, then `$` is not required
        if ((it != tmpEnd) && (it->edge() == 0b000) && (std::next(it) != tmpEnd))
            ++it;

        uint64_t previousEdge = 0b110; // impossible value
        while (it != tmpEnd)
        {
            uint64_t nextEdge = it->edge();

            if (nextEdge != previousEdge)
            {
                ++F[msb];

                locks.access(_i);
                if (notCurrentBlock(nextEdge))
                {
                    edges[_i] = nextEdge | 0b1000;
                    setToCurrent(nextEdge);
                } else {
                    edges[_i] = nextEdge;
                }
                ++_i;
                previousEdge = nextEdge;
            }

            ++it;
        }

        locks.access(_i-1);
        wplus[_i-1] = 1;

        _tmpValues.clear();
        ++C[msb];
    }

    GraphWriteCursor::GraphWriteCursor(WriteableGraph &graph, size_t i_, uint8_t colWidth_)
        : _BaseCursor(i_)
        , locks(graph.locks)
        , edges(graph.edges)
        , wplus(graph.wplus)
        , colourRegistry(graph.colourBuffers)
        , __k(graph.kmer_size())
        , colours()
        , zcx()
        , colourBitWidth(colWidth_)
    {
    }

    GraphWriteCursor::GraphWriteCursor(WriteableGraph &graph, size_t i_, uint8_t colWidth_, size_t reserve, ZstdCompressor *cx)
        : _BaseCursor(i_)
        , locks(graph.locks)
        , edges(graph.edges)
        , wplus(graph.wplus)
        , colourRegistry(graph.colourBuffers)
        , __k(graph.kmer_size())
        , colours(colWidth_)
        , zcx(cx)
        , colourBitWidth(colWidth_)
    {
        colours.reserve(reserve);
    }

    GraphWriteCursor::GraphWriteCursor(WriteableGraph &graph)
        : _BaseCursor()
        , locks(graph.locks)
        , edges(graph.edges)
        , wplus(graph.wplus)
        , colourRegistry(graph.colourBuffers)
        , __k(graph.kmer_size())
        , colours()
        , zcx()
        , colourBitWidth(0)
    {
    }

    GraphWriteCursor::GraphWriteCursor(GraphWriteCursor const &rhs, size_t i, size_t reserve, ZstdCompressor *cx)
        : _BaseCursor(i)
        , locks(rhs.locks, pv::SpinRegionManager::Accessor::shallow_copy_tag)
        , edges(rhs.edges)
        , wplus(rhs.wplus)
        , colourRegistry(rhs.colourRegistry)
        , __k(rhs.__k)
        , colours(rhs.colourBitWidth)
        , zcx(cx)
        , colourBitWidth(rhs.colourBitWidth)
    {
        colours.reserve(reserve);
    }

    GraphWriteCursor::GraphWriteCursor(GraphWriteCursor const &rhs, size_t i)
        : _BaseCursor(i)
        , locks(rhs.locks, pv::SpinRegionManager::Accessor::shallow_copy_tag)
        , edges(rhs.edges)
        , wplus(rhs.wplus)
        , colourRegistry(rhs.colourRegistry)
        , __k(rhs.__k)
        , colours()
        , zcx()
        , colourBitWidth(rhs.colourBitWidth)
    {
    }

    GraphWriteCursor::~GraphWriteCursor() {
        flushColours();
    }

    auto GraphWriteCursor::sortTempValues() {
        std::sort(
            _tmpValues.begin(),
            _tmpValues.end(),
            [](const BufferValue &lhs, const BufferValue &rhs) {
                return ((lhs.edge() << 61) | lhs.colour()) < ((rhs.edge() << 61) | rhs.colour());
            });
        return std::unique(_tmpValues.begin(), _tmpValues.end());
    }

    void GraphWriteCursor::flushColours() {
        if (!colours.empty()) {
            locks.releaseLock();

            compressColourBuffer();
            colourRegistry.saveColours(position(), colours, *zcx);
            resetColourBuffer();
        }
    }

    void GraphWriteCursor::writeNode(uint8_t msb) {
        auto tmpEnd = sortTempValues();
        auto it = _tmpValues.cbegin();

        // if there are valid edges other than `$`, then `$` is not required

        if ((it != tmpEnd) && (it->edge() == 0b000) && (std::next(it) != tmpEnd))
            ++it;

        // insert edge and colour values
        
        uint64_t previousEdge = 0b110; // impossible value
        while (it != tmpEnd)
        {
            uint64_t nextEdge = it->edge();

            if (nextEdge != previousEdge)
            {
                ++F[msb];

                locks.access(_i);
                if (notCurrentBlock(nextEdge))
                {
                    edges[_i] = nextEdge | 0b1000;
                    setToCurrent(nextEdge);
                } else {
                    edges[_i] = nextEdge;
                }
                ++_i;

                colours.boundaries.push_back(1);
                previousEdge = nextEdge;
            } else colours.boundaries.push_back(0);

            colours.colours.push_back(it->colour());
            ++it;
        }

        locks.access(_i-1);
        wplus[_i-1] = 1;
        _tmpValues.clear();
        ++(C[msb]);

        // check size of colour buffer and flush if too large

        if (colours.size() > MAX_BUFFER_SIZE)
            flushColours();
    }

    sdsl::int_vector<0> const& GraphWriteCursor::getColours() const noexcept {
        return colours.colours;
    }

    sdsl::bit_vector const& GraphWriteCursor::getBoundaries() const noexcept {
        return colours.boundaries;
    }

    ColourBufferRegistry::ColourBufferRegistry(std::string_view h5file_)
        : h5file(h5file_)
        , h5p(std::make_unique<HighFive::File>(
              h5file,
              HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate))
        , buffers()
    {
    }

    BaseGraph::BaseGraph(size_t _k, uint64_t _numGenomes, size_t _reserve)
        : k(_k)
        , numGenomes(_numGenomes)
        , locks(_reserve, LOCK_SAMPLE_RATE)
        , edges(_reserve)
        , wplus(_reserve)
    {
    }

    BaseGraph::BaseGraph(uint8_t k_,
                         uint64_t _numGenomes,
                         sdsl::int_vector<4> &&edges_,
                         sdsl::bit_vector &&wplus_,
                         std::array<size_t,5> F_,
                         std::array<size_t,5> C_)
        : k(k_)
        , numGenomes(_numGenomes)
        , locks(edges_.size(), LOCK_SAMPLE_RATE)
        , edges(std::move(edges_))
        , wplus(std::move(wplus_))
        , F(F_)
        , C(C_)
    {
    }

    WriteableGraph::WriteableGraph(size_t _k,
                                   uint64_t _maxColour,
                                   uint64_t _numGenomes,
                                   uint8_t baseColourWidth_,
                                   size_t _reserve)
        : BaseGraph(_k, _numGenomes, _reserve)
        , maxColour(_maxColour)
        , colourBitWidth(required_bits(_maxColour))
        , baseIdBitWidth(baseColourWidth_)
        , colourBuffers(ColourBufferRegistry::temp_flag)
    {
    }

    WriteableGraph::WriteableGraph(size_t _k,
                                   uint64_t _maxColour,
                                   uint64_t _numGenomes,
                                   uint8_t baseColourWidth_,
                                   size_t _reserve,
                                   std::string _bufferPath)
        : BaseGraph(_k, _numGenomes, _reserve)
        , maxColour(_maxColour)
        , colourBitWidth(required_bits(_maxColour))
        , baseIdBitWidth(baseColourWidth_)
        , colourBuffers(_bufferPath)
    {
    }

    TraversableGraph::TraversableGraph(BaseGraph &&g)
        : k(g.k)
        , maxColour(1)
        , numGenomes(g.numGenomes)
        , colourBitWidth(1)
        , baseIdBitWidth(1)
        , last(std::move(g.wplus))
        , lastRankSupport(&last)
        , lastSelectSupport(&last)
        , W()
        , colourBuffers()
        , F(g.F)
        , C(g.C)
    {
        sdsl::construct_im(W, std::move(g.edges));
    }

    TraversableGraph::TraversableGraph(BaseGraph &&g, fs::path const &tempFile)
        : k(g.k)
        , maxColour(1)
        , numGenomes(g.numGenomes)
        , colourBitWidth(1)
        , baseIdBitWidth(1)
        , last(std::move(g.wplus))
        , lastRankSupport(&last)
        , lastSelectSupport(&last)
        , W()
        , colourBuffers()
        , F(g.F)
        , C(g.C)
    {
        // serialise array to disk
        sdsl::store_to_file(std::move(g.edges), tempFile);

        // construct wavelet tree
        sdsl::memory_monitor::start();
        sdsl::construct(W, tempFile, 0);
        sdsl::memory_monitor::stop();
        sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(LOG(INFO));

        // remove temporary file
        fs::remove(tempFile);
    }

    TraversableGraph::TraversableGraph(graph::WriteableGraph &&g)
        : k(g.k)
        , maxColour(g.maxColour)
        , numGenomes(g.numGenomes)
        , colourBitWidth(g.colourBitWidth)
        , baseIdBitWidth(g.baseIdBitWidth)
        , last(std::move(g.wplus))
        , lastRankSupport(&last)
        , lastSelectSupport(&last)
        , W()
        , colourBuffers(std::move(g.colourBuffers))
        , F(g.F)
        , C(g.C)
    {
        sdsl::construct_im(W, std::move(g.edges));
    }

    TraversableGraph::TraversableGraph(WriteableGraph &&g, fs::path const &tempFile)
        : k(g.k)
        , maxColour(g.maxColour)
        , numGenomes(g.numGenomes)
        , colourBitWidth(g.colourBitWidth)
        , baseIdBitWidth(g.baseIdBitWidth)
        , last(std::move(g.wplus))
        , lastRankSupport(&last)
        , lastSelectSupport(&last)
        , W()
        , colourBuffers(std::move(g.colourBuffers))
        , F(g.F)
        , C(g.C)
    {
        // serialise array to disk
        sdsl::store_to_file(std::move(g.edges), tempFile);

        // construct wavelet tree
        sdsl::memory_monitor::start();
        sdsl::construct(W, tempFile, 0);
        sdsl::memory_monitor::stop();
        sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(LOG(INFO));

        // remove temporary file
        fs::remove(tempFile);
    }

    TraversableGraph::TraversableGraph(uint8_t k,
                                       uint64_t maxColour,
                                       uint64_t numGenomes,
                                       uint8_t colourBitWidth,
                                       uint8_t baseIdBitWidth,
                                       ColourBufferRegistry &&colourBuffers,
                                       std::array<size_t, 5> const &F,
                                       std::array<size_t, 5> const &C)
        : k(k)
        , maxColour(maxColour)
        , numGenomes(numGenomes)
        , colourBitWidth(colourBitWidth)
        , baseIdBitWidth(baseIdBitWidth)
        , colourBuffers(std::move(colourBuffers))
        , F(F)
        , C(C)
    {
    }

    TraversableGraph::ColourBufferRange::ColourBufferRange(const ColourBuffer &buffer, size_t grainsize, size_t edgeEnd)
        : _p(buffer.cbegin())
        , rp(std::make_shared<sdsl::rank_support_v5<1,1>>(&buffer.boundaries))
        , sp(std::make_shared<sdsl::select_support_mcl<1,1>>(&buffer.boundaries))
        , _begin(0)
        , _end(buffer.size())
        , _grainsize(grainsize)
        , _edgeOffset(edgeEnd - load())
    {
    }

    TraversableGraph::ColourBufferRange::ColourBufferRange(const ColourBuffer &buffer, size_t grainsize)
        : _p(buffer.cbegin())
        , rp(std::make_shared<sdsl::rank_support_v5<1,1>>(&buffer.boundaries))
        , sp(std::make_shared<sdsl::select_support_mcl<1,1>>(&buffer.boundaries))
        , _begin(0)
        , _end(buffer.size())
        , _grainsize(grainsize)
        , _edgeOffset(0)
    {}

    TraversableGraph::ColourBufferRange::ColourBufferRange(const ColourBufferRange &rhs, size_t begin, size_t end)
        : _p(rhs._p)
        , rp(rhs.rp)
        , sp(rhs.sp)
        , _begin(begin)
        , _end(end)
        , _grainsize(rhs._grainsize)
        , _edgeOffset(rhs._edgeOffset)
    {
    }

    TraversableGraph::ColourBufferRange::ColourBufferRange(ColourBufferRange &lhs, oneapi::tbb::proportional_split p)
        : ColourBufferRange(lhs)
    {
        size_t leftEdges = lhs.load() * p.left() / (p.left() + p.right());
        _begin = sp->select(rp->rank(_begin) + leftEdges);
        lhs._end = _begin;
    }

    TraversableGraph::ColourBufferRange::ColourBufferRange(ColourBufferRange &lhs, oneapi::tbb::split)
        : ColourBufferRange(lhs, oneapi::tbb::proportional_split(1, 1))
    {}

    TraversableGraph::ColourBufferRange
    TraversableGraph::ColourBufferRange::subrange(size_t begin, size_t end) const {
        size_t rnkOffset = rp->rank(_begin);
        // `begin` and `end` need to be converted to 1-based indices
        return ColourBufferRange(*this, sp->select(rnkOffset + begin + 1u), sp->select(rnkOffset + end + 1u));
    }

    ColourBufferConstIterator TraversableGraph::ColourBufferRange::begin() const
    { return _p + _begin; }

    ColourBufferConstIterator TraversableGraph::ColourBufferRange::end() const
    { return _p + _end; }

    bool TraversableGraph::ColourBufferRange::empty() const
    { return _begin == _end; }

    size_t TraversableGraph::ColourBufferRange::size() const
    { return _end - _begin; }

    size_t TraversableGraph::ColourBufferRange::load() const
    { return rp->rank(_end) - rp->rank(_begin); }

    bool TraversableGraph::ColourBufferRange::is_divisible() const
    { return load() >= (2 * _grainsize); }

    size_t TraversableGraph::ColourBufferRange::edgeBegin() const
    { return _edgeOffset + rp->rank(_begin); }

    size_t TraversableGraph::ColourBufferRange::edgeEnd() const
    { return _edgeOffset + rp->rank(_end); }

    ColourBufferConstIterator TraversableGraph::ColourBufferRange::selectEdge(size_t i) const
    { return _p + sp->select(i - _edgeOffset); }

    // k-mer sorting ---------------------------------------------------------------------------

    // bulk fill
    // ---------

    void fillGraphBulk(GraphWriteCursor &crs, const Dna4GenomeVector &genomes_, uint32_t threads_)
    {
        uint8_t
            k = crs.k(),
            valueBytes = _mkimem_details::value_size(genomes_.getFeatureWidth());

        // extract k-mers

        auto blocks = accumulatePartials(genomes_.buffer(), [&](Dna4Genome const &g)->size_t
            { return g.numKmers(k); });
        
        mers::KmerBuffer kmers(blocks.back(), valueBytes, k);
        kmers.insertKmers(genomes_, UnwindGenome{ k }, blocks);
        kmers.sort(threads_);

        // extract terminals

        terminals::TerminalBuffer tmnl = terminals::extractTerminalsAndSort(genomes_, UnwindGenome{}, k);
        auto tmnlRange = tmnl.asRange();

        flush::flush(crs, kmers, tmnlRange, WriteableGraph::LOCK_SAMPLE_RATE);
    }

    void fillGraphBulk(BaseWriteCursor &crs, ChunkedDna4Genome const &genome_, uint32_t threads_)
    {
        uint8_t k = crs.k();

        // extract k-mers

        auto blocks = accumulatePartials<SequenceRange>(genome_, [&](SequenceRange const &rng)->size_t
            { return rng.numKmers(k); });
        
        mers::UncolouredKmerBuffer kmers(blocks.back(), k);
        kmers.insertKmers(genome_, UnwindChunkedGenome{}, blocks);
        kmers.sort(threads_);

        // extract terminals

        terminals::TerminalBuffer tmnl = terminals::extractTerminalsAndSort(genome_, UnwindChunkedGenome{}, k);
        auto tmnlRange = tmnl.asRange();

        flush::flush(crs, kmers, tmnlRange, BaseGraph::LOCK_SAMPLE_RATE);
    }

    void fillGraphBulk(BaseWriteCursor &crs, std::vector<reads::ReadChunks> const &reads_, uint32_t threads_)
    {
        uint8_t k = crs.k();

        // extract k-mers and terminals
        
        auto blocks = accumulatePartials(reads_, reads::detail::numEdges);

        terminals::TerminalBuffer kmers(blocks.back(), k, terminals::TerminalBuffer::autofit_tag);
        kmers.insertKmers(reads_, reads::UnwindReads{}, blocks);
        kmers.sort(threads_);

        flush::flush(crs, kmers, BaseGraph::LOCK_SAMPLE_RATE);
    }

    // suffix-wise fill
    // ----------------

    void fillGraphBySuffix(GraphWriteCursor &crs, const Dna4GenomeVector &genomes_, uint8_t s_, uint32_t threads_)
    {
        uint8_t
            k = crs.k(),
            k_eff = k - s_,
            valueBytes = _mkimem_details::value_size(genomes_.getFeatureWidth());

        UnwindGenome Unwind{ k };

        // count suffixes in each genome

        auto genomeSuffixCounts = suffix::countSuffixes(genomes_, Unwind, s_, k - s_);
        suffix::SuffixTable suffixTotals = suffix::accumulate(genomeSuffixCounts);
        size_t requiredBufferSize = suffixTotals.maxValue();

        // allocate buffer for k-mer extraction

        mers::KmerBuffer
            kmers(requiredBufferSize, valueBytes, k, k_eff),
            temp(requiredBufferSize, valueBytes, k, k_eff);

        // extract sequence terminals

        terminals::TerminalBuffer tmnl = terminals::extractTerminalsAndSort(genomes_, UnwindGenome{}, k);
        auto tmnlRange = tmnl.asRange();

        // suffix-wise insertion

        detail::suffixwiseFill(crs,
                               genomes_, Unwind, genomeSuffixCounts,
                               kmers, temp,
                               tmnlRange,
                               suffix::SmallRollingNuclSeq(0), s_,
                               threads_);
    }

    void fillGraphBySuffix(BaseWriteCursor &crs, ChunkedDna4Genome const &genome_, uint8_t s_, uint32_t threads_)
    {
        uint8_t
            k = crs.k(),
            k_eff = k - s_;

        // count suffixes in each genome

        auto regionSuffixCounts = suffix::countSuffixes(genome_, UnwindChunkedGenome{}, s_, k - s_);
        suffix::SuffixTable suffixTotals = suffix::accumulate(regionSuffixCounts);
        size_t requiredBufferSize = suffixTotals.maxValue();

        // allocate buffer for k-mer extraction

        mers::UncolouredKmerBuffer
            kmers(requiredBufferSize, k, k_eff),
            temp(requiredBufferSize, k, k_eff);

        // extract sequence terminals

        terminals::TerminalBuffer tmnl = terminals::extractTerminalsAndSort(genome_, UnwindChunkedGenome{}, k);
        auto tmnlRange = tmnl.asRange();

        // suffix-wise insertion

        detail::suffixwiseFill(crs,
                               genome_, UnwindChunkedGenome{}, regionSuffixCounts,
                               kmers, temp,
                               tmnlRange,
                               suffix::SmallRollingNuclSeq(0), s_,
                               threads_);
    }

    void fillGraphBySuffix(BaseWriteCursor &crs, std::vector<reads::ReadChunks> const &reads_, uint8_t s_, uint32_t threads_)
    {
        assert(s_ > 0);

        uint8_t
            k = crs.k(),
            k_eff = k - s_;

        // count suffixes in each chunk of reads
        
        auto chunkSuffixCounts = suffix::countSuffixes(reads_, reads::UnwindReads{}, s_);
        suffix::SuffixTable suffixTotals = suffix::accumulate(chunkSuffixCounts);
        size_t requiredBufferSize = suffixTotals.maxValue();

        // allocate buffer for k-mer extraction

        terminals::TerminalBuffer
            kmers(requiredBufferSize, k, k_eff, terminals::TerminalBuffer::autofit_tag),
            temp(requiredBufferSize, k, k_eff, terminals::TerminalBuffer::autofit_tag);

        // extract short terminals

        terminals::TerminalBuffer tmnl = terminals::extractTerminalsAndSort(reads_, reads::UnwindReads{}, s_);
        auto tmnlRange = tmnl.asRange();

        // suffix-wise insertion

        detail::suffixwiseFill(crs,
                               reads_, reads::UnwindReads{}, chunkSuffixCounts,
                               kmers, temp,
                               tmnlRange,
                               suffix::SmallRollingNuclSeq(0), s_,
                               threads_);
    }

    void makeLastBufferOnDisk(sdsl::bit_vector &wplus, StaticGraphDiskFileConfig const &cfg)
    {
        // last - Rank support

        {
            LOG(INFO) << "serialising lastRankSupport to " << cfg.lastRank;
            sdsl::rank_support_v5<1,1> lRnk;
            sdsl::util::init_support(lRnk, &wplus);
            sdsl::store_to_file(std::move(lRnk), cfg.lastRank);
        }

        // last - Select support

        {
            LOG(INFO) << "serialising lastSelectSupport to " << cfg.lastSelect;
            sdsl::select_support_mcl<1,1> lSel;
            sdsl::util::init_support(lSel, &wplus);
            sdsl::store_to_file(std::move(lSel), cfg.lastSelect);
        }

        // last

        LOG(INFO) << "serialising last to " << cfg.last;
        fileutils::store_to_file(std::move(wplus), cfg.last);
    }

    void makeEdgesBufferOnDisk(sdsl::int_vector<4> &edges, StaticGraphDiskFileConfig const &cfg)
    {
        LOG(INFO) << "constructing wavelet tree for W";

        // size of original data structure
        LOG(INFO)
            << "original size of Edges: "
            << static_cast<long double>(sdsl::size_in_bytes(edges)) / 1.0e9l
            << " Gb";

        // serialise array to disk
        sdsl::store_to_file(std::move(edges), cfg.tempW);

        sdsl::wt_int<> W;

        // construct wavelet tree
        //sdsl::memory_monitor::start();
        sdsl::construct(W, cfg.tempW, 0);
        //sdsl::memory_monitor::stop();
        //sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(LOG(INFO));

        // remove temporary file
        fs::remove(cfg.tempW);

        // size of original data structure
        LOG(INFO)
            << "size of wavelet tree W: "
            << static_cast<long double>(sdsl::size_in_bytes(W)) / 1.0e9l
            << " Gb";

        // write to disk
        LOG(INFO) << "serialising W to " << cfg.W;
        sdsl::store_to_file(W, cfg.W);
    }

    void makeStaticGraphOnDisk(BaseGraph &g, fs::path const &databaseDir)
    {
        StaticGraphDiskFileConfig config(databaseDir);

        makeLastBufferOnDisk(g.wplus, config);
        makeEdgesBufferOnDisk(g.edges, config);

        // serialise graph metadata

        {
            TraversableGraph graph(g.k, 1, g.numGenomes, 1, 1, {}, g.F, g.C);
            saveSmallBuffers(config.meta, graph);
        }
    }

    void makeStaticGraphOnDisk(WriteableGraph &g, fs::path const &databaseDir)
    {
        StaticGraphDiskFileConfig config(databaseDir);

        makeLastBufferOnDisk(g.wplus, config);
        makeEdgesBufferOnDisk(g.edges, config);

        // serialise graph metadata

        {
            TraversableGraph graph(g.k, g.maxColour, g.numGenomes, g.colourBitWidth, g.baseIdBitWidth, std::move(g.colourBuffers), g.F, g.C);
            saveSmallBuffers(config.meta, graph);
        }
    }

    uint8_t TraversableGraph::block(size_t i) const
    {
        if      (F[1] > i) return 0u|0b1000u;
        else if (F[2] > i) return 1u|0b1000u;
        else if (F[3] > i) return 2u|0b1000u;
        else if (F[4] > i) return 3u|0b1000u;
        else               return 4u|0b1000u;
    }

    uint8_t TraversableGraph::edge(size_t i) const
    {
        return W[i] & 0b0111;
    }

    Dna4Sequence TraversableGraph::kmer(size_t i) const
    {
        size_t remaining = k+1;
        Dna4Sequence s;
        s.reserve(remaining);
        std::optional i_ = i;
        // check first edge, it may be $
        if (auto e_ = edge(*i_); e_ > 0)
            s.push_back(parsing::dna5ToDna4(e_));
        i_ = bwd(*i_);
        --remaining;
        // continue
        while ((i_) && (remaining)) {
            s.push_back(parsing::dna5ToDna4(edge(*i_)));
            i_ = bwd(*i_);
            --remaining;
        }
        auto rs = s | std::views::reverse;
        return Dna4Sequence(rs);
    }

    // WARNING: no bounds check on input
    std::optional<size_t> TraversableGraph::fwd(size_t i) const
    {
        uint8_t e = edge(i);
        if (e == 0u)
            return std::nullopt;
        size_t rnk = W.rank(i + 1, e|0b1000u);
        return std::make_optional(lastSelectSupport(lastRankSupport(F[e]) + rnk));
    }

    // WARNING: no bounds check on input
    std::optional<size_t> TraversableGraph::bwd(size_t i) const
    {
        uint8_t c = block(i);
        if ((c^0b1000u) == 0u)
            return std::nullopt;
        size_t r1 = (lastRankSupport(i)+1u), r2 = lastRankSupport(F[c^0b1000u]);
        if (r1 < r2)
            LOG(ERROR) << "r1<r2...i=" << i << "/" << last.size()-1u << ", c=" << (size_t)c << ", F[c]=" << F[c^0b1000u];
        return W.select(r1-r2, c);
    }

    std::optional<size_t> TraversableGraph::outgoing(size_t i, uint8_t e) const
    {
        IndexRange r = getNode(i);

        size_t j = WPred(r.end, e | 0b1000u);
        if (j >= r.begin)
        {
            return std::make_optional(j);
        }
        else
        {
            j = WPred(r.end, e & 0b0111u);
            if (j >= r.begin)
                return std::make_optional(j);
            else
                return std::nullopt;
        }
    }

    sdsl::bit_vector TraversableGraph::getTerminals() const
    {
        // construct a bit-vector representing position of terminals

        // number of nodes in the graph
        size_t m = sdsl::util::cnt_one_bits(last);

        // bit-vector where 1 designates a terminal

        sdsl::bit_vector bv(m);

        // fill bit vector

        std::vector<_IndexRangeDepth> pfxs;
        pfxs.reserve(4);

        IndexRange root = getRoot();
        pfxs.emplace_back(root.begin, root.end, 0u);

        while ((!pfxs.empty()) && (pfxs.size() < 4))
        {
            _IndexRangeDepth r = pfxs.front();
            pfxs.erase(pfxs.begin());

            bv[lastRankSupport(r.end) - 1u] = 1;

            if (r.depth < (k - 1u))
            {
                while (!r.empty())
                {
                    IndexRange chld = getNode(fwd(r.begin).value());
                    pfxs.emplace_back(chld.begin, chld.end, r.depth + 1);
                    ++r.begin;
                }
            }
        }

        oneapi::tbb::parallel_for(oneapi::tbb::blocked_range(std::make_move_iterator(pfxs.begin()),
                                                             std::make_move_iterator(pfxs.end())),
                                  FillTerminalBitVector(this, &bv));

        LOG(INFO) << "found " << sdsl::util::cnt_one_bits(bv) << " terminals out of " << m << " nodes";

        return bv;
    }

    std::mutex lockForBV;

    void
    TraversableGraph::FillTerminalBitVector::operator()(_IndexRangeDepth &&b_) const
    {
        std::vector<size_t> cache;

        _IndexRangeDepth_Stack pfxs;
        pfxs.emplace(std::move(b_));

        while (!pfxs.empty())
        {
            _IndexRangeDepth &rng = pfxs.top();

            if ((rng.depth < (g->k - 1u)) && (!rng.empty()))
            {
                IndexRange chld = g->getNode(g->fwd(rng.begin++).value());
                pfxs.emplace(chld.begin, chld.end, rng.depth + 1);
            }
            else
            {
                cache.push_back(g->lastRankSupport(rng.end) - 1u);
                if (cache.size() == maxCacheSize)
                    flushCache(cache);
                pfxs.pop();
            }
        }

        flushCache(cache);
    }

    void
    TraversableGraph::FillTerminalBitVector::flushCache(std::vector<size_t> &cache) const
    {
        std::sort(cache.begin(), cache.end());
        {
            std::lock_guard<std::mutex> guard(lockForBV);
            for (size_t const &i : cache)
                (*bv)[i] = 1;
        }
        cache.clear();
    }

    TraversableGraph::SuccinctTerminalIndex::SuccinctTerminalIndex(const TraversableGraph &g_)
        : sdv(g_.getTerminals()),
          sdvRank(&sdv),
          sdvSelect(&sdv)
    {
    }

    TraversableGraph::KmerRangeGenerator::KmerRangeGenerator(TraversableGraph const &g_)
        : sdv(std::make_shared<SuccinctTerminalIndex>(g_)),
          lastRank(g_.lastRankSupport),
          lastSelect(g_.lastSelectSupport),
          i(0u),
          endI(g_.last.size()),
          tI(1u),
          maxTerminal(sdv->rank(sdv->size()))
    {
    }

    TraversableGraph::KmerRangeGenerator::KmerRangeGenerator(KmerRangeGenerator const &obj_, size_t begin, size_t end)
        : sdv(obj_.sdv),
          lastRank(obj_.lastRank),
          lastSelect(obj_.lastSelect),
          i(begin),
          endI(end),
          tI(sdv->rank(lastRank(begin)) + 1u),
          maxTerminal(sdv->rank(lastRank(end)))
    {
    }

    std::optional<TraversableGraph::IndexRange>
    TraversableGraph::KmerRangeGenerator::next()
    {
        while (i < endI)
        {
            if (tI <= maxTerminal)
            {
                size_t
                    tNI = sdv->select(tI),
                    tEI = std::min<size_t>(tNI == 0u ? 0u : lastSelect(tNI) + 1u,
                                           endI),
                    tEJ = lastSelect(tNI + 1u) + 1u;

                auto r = IndexRange(i, tEI);

                ++tI;
                i = tEJ;

                if (!r.empty())
                    return r;
            }
            else
            {
                auto r = IndexRange(i, endI);
                i = endI;
                return r;
            }
        }
        return std::nullopt;
    }

} // namespace graph
