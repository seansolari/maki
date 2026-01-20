#include <syncstream>
#include <zstr.hpp>
#include <maki/dist.hpp>
#include <maki/timing.hpp>

namespace dist
{

    namespace pairwise
    {

        // counting shared kmers (matrix-based) --------------------------------------------------------------------------- //

        MatrixCountColourPairs::MatrixCountColourPairs(Matrix *m_, bitmask proj_)
            : m(m_)
            , proj(proj_)
        {
        }

        void MatrixCountColourPairs::operator()(const ColourBufferRange &r) const
        {
            ColourVector tmp;

            auto it = r.begin(),
                 end = r.end();
            while (it != end) {
                do {
                    tmp.push_back(proj(*it.colour));
                    ++it;
                } while ((it != end) && (*it.boundary == 0));

                std::sort(tmp.begin(), tmp.end());
                auto end = std::unique(tmp.begin(), tmp.end());
                tmp.resize(end - tmp.begin());

                flushColours(tmp);
            }
        }

        void MatrixCountColourPairs::flushColours(ColourVector &colours) const
        {
            // write colours to temp buffer

            for (size_t i = 1; i < colours.size(); ++i)
            {
                for (size_t j = 0; j < i; ++j)
                {
                    // values are sorted: j < i => colours[j] < colours[i]
                    ++(*m)(colours[i], colours[j]);
                }
            }

            // clear input buffer

            colours.clear();
        }

        Matrix countPairsDense(const graph::TraversableGraph &g)
        {
            size_t numColours = g.getNumGenomes() + 1;
            Matrix counts(numColours);
            ApplyAcrossColours(g.colourBuffers, 1, MatrixCountColourPairs(&counts, bitmask(g.getBaseColourBitWidth())));
            return counts;
        }

        SampleColourPairs::SampleColourPairs(ColourSampler *obj_, uint8_t w_) : obj(obj_) , keyWidth(w_), msk(keyWidth) {}

        SampleColourPairs::SampleColourPairs(SampleColourPairs &rhs, oneapi::tbb::split) : SampleColourPairs(rhs) {}

        void SampleColourPairs::operator()(const ColourBufferRange &r) const {
            ColourSampler::Accessor data(*obj);
            ColourVector tmp;
            auto it = r.begin(), end = r.end();
            size_t aBig, aSmall, bBig, bSmall;
            while (it != end) {
                do {
                    tmp.push_back(*it.colour);
                    ++it;
                } while ((it != end) && (*it.boundary == 0));
                for (size_t i = 1; i < tmp.size(); ++i) {
                    aBig = msk(tmp[i]);
                    aSmall = tmp[i] >> keyWidth;
                    if (aSmall) {
                        for (size_t j = 0; j < i; ++j) {
                            bBig = msk(tmp[j]);
                            bSmall = tmp[j] >> keyWidth;
                            if (bSmall) {
                                data.increment(bBig-1, aBig-1, bSmall|(aSmall<<32));
                                data.increment(aBig-1, bBig-1, aSmall|(bSmall<<32));
                            }
                        }
                    }
                }
                tmp.clear();
            }
        }

        std::unordered_set<size_t> ColourSampler::gatherKeys(uint8_t baseWidth) {
            assert(baseWidth < 32u);
            std::unordered_set<size_t> keys;
            forEach([&](size_t gid1, size_t gid2, const KeyValuePair &kv){
                keys.insert(gid1|((kv.key&0xFFFFFFFFULL)<<baseWidth));
                keys.insert(gid2|((kv.key&~0xFFFFFFFFULL)>>(32u-baseWidth)));
            });
            return keys;
        }

        ColourSampler samplePairsDense(const graph::TraversableGraph &g) {
            ColourSampler counts(g.getNumGenomes());
            ApplyAcrossColours(g.colourBuffers, 1, SampleColourPairs(&counts, g.getBaseColourBitWidth()));
            return counts;
        }

        void writeToCompressed(const Matrix &m, std::string filepath)
        {
            zstr::ofstream fh(filepath);

            size_t N = m.size();
            for (size_t r = 1; r < N; ++r) {
                for (size_t c = 0; c < r; ++c) {
                    auto const &v = m(r, c);
                    if (v)
                        fh << r << '\t' << c << '\t' << (uint32_t)v << '\n';
                }
            }

            fh.close();
            LOG(INFO) << "pairwise counts written to " << filepath;
        }

        void writeToTsv(const Matrix &m, std::string filepath)
        {
            std::ofstream fh(filepath);

            size_t N = m.size();
            for (size_t r = 1; r < N; ++r) {
                for (size_t c = 0; c < r; ++c) {
                    auto const &v = m(r, c);
                    if (v)
                        fh << r << '\t' << c << '\t' << (uint32_t)v << '\n';
                }
            }

            fh.close();
            LOG(INFO) << "pairwise counts written to " << filepath;
        }

    } // namespace pairwise

    namespace query
    {

        // counting shared k-mers (dict) ------------------------------------------------------------------- //

        size_t QuerySet::size() const noexcept
        { return keys.size(); }

        void QuerySet::insert(uint64_t k_) {
            int64_t v = static_cast<int64_t>(size());
            keys.try_emplace(k_, v);
        }

        int64_t QuerySet::index(uint64_t k_) const {
            if (auto it = keys.find(k_); it != keys.end()) {
                return it->second;
            } else return -1;
        }

        void ConcurrentCounter::incrementRange(uint64_t lhs, sdsl::int_vector<0>::const_iterator rhsIt, sdsl::int_vector<0>::const_iterator rhsEnd) {
            while (rhsIt != rhsEnd) {
                uint64_t rhs = *rhsIt++;
                if (lhs != rhs) {
                    auto [it, allc] = data.emplace(copyable_pair<uint64_t>(lhs, rhs), (uint32_t)1u);
                    if (!allc) ++it->second;
                }
            }
        }

        ConcurrentBuffers::ConcurrentBuffers(size_t n_) : data(n_) {}

        size_t ConcurrentBuffers::size() const noexcept
        { return data.size(); }

        std::vector<ConcurrentCounter>::const_iterator ConcurrentBuffers::begin() const
        { return data.begin(); }

        std::vector<ConcurrentCounter>::const_iterator ConcurrentBuffers::end() const
        { return data.end(); }

        std::unordered_set<uint64_t> ConcurrentBuffers::gatherKeys() const {
            std::unordered_set<size_t> keys;
            for (const auto &buffer : data) {
                for (const auto &[kp, v] : buffer.data) {
                    keys.insert(kp.first);
                    keys.insert(kp.second);
                }
            }
            return keys;
        }

        CountQueries::CountQueries(const QuerySet &q_, ConcurrentBuffers *b_, uint8_t w_)
            : queries(q_)
            , buffers(b_)
            , msk(w_)
        {}

        void CountQueries::operator()(const ColourBufferRange &r) const {
            auto prev = r.begin(),
                 current = prev,
                 end = r.end();

            while (current != end)
            {
                do ++current; while ((current != end) && (*current.boundary == 0));
                searchQueries(prev.colour, current.colour);
                prev = current;
            }
        }

        void CountQueries::searchQueries(const sdsl::int_vector<0>::const_iterator &begin, const sdsl::int_vector<0>::const_iterator &end) const {
            for (auto it = begin; it != end; ++it) {
                int64_t bufferIndex = queries.index(msk(*it));
                if (bufferIndex != -1) {
                    ConcurrentCounter &buffer = (*buffers)[bufferIndex];
                    buffer.incrementRange(*it, begin, end);
                }
            }
        }

        ConcurrentBuffers searchGenomes(const QuerySet &queries, const graph::TraversableGraph &g) {
            ConcurrentBuffers buffers(queries.size());
            CountQueries Op(queries, &buffers, g.getBaseColourBitWidth());
            ApplyAcrossColours(g.colourBuffers, 1, Op);
            return buffers;
        }
        
    }  // namespace query
    
    namespace unique
    {

        // counting unique kmers --------------------------------------------------------------------------- //
        
        void writeToDisk(const Map &m, std::string filepath)
        {
            if (filepath.ends_with(".gz"))
                writeToCompressed(m, filepath);
            else
                writeToTsv(m, filepath);
        }

        void writeToCompressed(const Map &m, std::string filepath)
        {
            zstr::ofstream fh(filepath);

            for (size_t i = 0; i < m.size(); ++i)
                if ((uint32_t)m[i] > 0)
                    fh << i << '\t' << (uint32_t)m[i] << '\n';

            fh.close();
            LOG(INFO) << "counts written to " << filepath;
        }

        void writeToTsv(const Map &m, std::string filepath)
        {
            std::ofstream fh(filepath);

            for (size_t i = 0; i < m.size(); ++i)
                if ((uint32_t)m[i] > 0)
                    fh << i << '\t' << (uint32_t)m[i] << '\n';

            fh.close();
            LOG(INFO) << "counts written to " << filepath;
        }

    } // namespace unique

    namespace cooc
    {
        
        size_t Xlist::Hash::operator()(const size_t &index) const {
            size_t seed = arr._itemWidth;
            for (size_t i = index * arr._itemWidth; i < (index + 1) * arr._itemWidth; ++i)
                seed ^= std::hash<typename sdsl::int_vector<0>::value_type>{}(arr._keys[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }

        bool Xlist::Compare::operator()(const size_t &lhs, const size_t &rhs) const {
            for (size_t i = lhs * arr._itemWidth, j = rhs * arr._itemWidth;
                 i < (lhs + 1) * arr._itemWidth;
                 ++i, ++j)
                if (arr._keys[i] != arr._keys[j]) return false;
            return true;
        }

        Xlist::Xlist(size_t bitWidth_, size_t itemWidth_) : Xlist(bitWidth_, itemWidth_, 0) {}

        Xlist::Xlist(size_t bitWidth_, size_t itemWidth_, size_t reserve_)
            : _items(0)
            , _itemWidth(itemWidth_)
            , _keys(0, 0, bitWidth_)
            , _hash(*this)
            , _compare(*this)
            , _map(std::max(reserve_, (size_t)10), _hash, _compare)
        {
            if (reserve_)
                _keys.reserve(_itemWidth * reserve_);
        }

        bool Xlist::update()
        { return ++_map[_items] == (uint32_t)1; }

        RankMap::RankMap(size_t bitWidth_)
            : _bypk((bitWidth_ + 7) / 8)
            , _data()
        {}

        size_t RankMap::ByteVectorHash::operator()(const Key& v) const noexcept {
            std::string_view sv(reinterpret_cast<const char*>(v.data()), v.size());
            return std::hash<std::string_view>{}(sv);
        }

        size_t RankMap::ByteVectorHash::operator()(std::string_view sv) const noexcept {
            return std::hash<std::string_view>{}(sv);
        }

        bool RankMap::ByteVectorEq::operator()(const Key& a, const Key& b) const noexcept {
            return a == b;
        }

        bool RankMap::ByteVectorEq::operator()(const Key& a, std::string_view b) const noexcept {
            std::string_view sva(reinterpret_cast<const char*>(a.data()), a.size());
            return sva == b;
        }

        RankMap::Key RankMap::compressU64(const ColourVector& src) const {
            Key out(src.size() * _bypk);
            const uint8_t* p = reinterpret_cast<const uint8_t*>(src.data());
            uint8_t* o = out.data();
            for (size_t i = 0, stride = 8; i < src.size(); ++i, o += _bypk) {
                std::memcpy(o, p + i * stride, _bypk);
            }
            return out;
        }

        ColourVector RankMap::decompressU64(const Key& bytes) const {
            const size_t n = bytes.size() / _bypk;
            ColourVector out(n, 0);
            const uint8_t* src = bytes.data();
            uint8_t* p = reinterpret_cast<uint8_t*>(out.data());
            for (size_t i = 0, stride = 8; i < n; ++i) {
                std::memcpy(p + i * stride, src, _bypk);
                src += _bypk;
            }
            return out;
        }

        void RankMap::increment(const ColourVector &v) {
            auto key = compressU64(v);
            std::string_view keyView(reinterpret_cast<const char*>(key.data()), key.size());
            _data.lazy_emplace_l(
                keyView,
                [](Map::value_type& kv) { ++kv.second; },
                [&](const Map::constructor& ctor) { ctor(std::move(key), /*initial*/ 1u); }
            );
        }

        void RankMap::set(const ColourVector &v, uint32_t z) {
            auto key = compressU64(v);
            std::string_view keyView(reinterpret_cast<const char*>(key.data()), key.size());
            _data.lazy_emplace_l(
                keyView,
                [&](Map::value_type& kv) { kv.second = z; },
                [&](const Map::constructor& ctor) { ctor(std::move(key), /*initial*/ z); }
            );
        }

        RankMap::Map::const_iterator RankMap::find(const ColourVector &v) const {
            auto key = compressU64(v);
            std::string_view keyView(reinterpret_cast<const char*>(key.data()), key.size());
            return _data.find(key);
        }

        RankMap::Map::const_iterator RankMap::end() const
        { return _data.end(); }

        void RankMap::flush(std::ostream &bos) const {
            oneapi::tbb::parallel_for((size_t)0, _data.subcnt(), (size_t)1, [&](size_t submapIndex)->void {
                std::osyncstream os(bos);
                _data.with_submap(submapIndex, [&](const Map::EmbeddedSet& set)->void {
                    for (const auto &[key, count] : set) {
                        auto bigKey = decompressU64(key);
                        auto it_ = bigKey.cbegin(), end_ = bigKey.cend();
                        os << *it_++;
                        while (it_ != end_)
                            os << ',' << *it_++;
                        os << '\t' << count << '\n';
                    }
                });
            });
        }

        bool RankMap::operator==(const RankMap &rhs) const
        { return _data == rhs._data; }

        BaseColourCounter::BaseColourCounter(RankMap *m_) : m(m_) {}

        CountColourLists::CountColourLists(RankMap *m_, uint8_t w_) : BaseColourCounter(m_), msk(w_) {}

        void CountColourLists::operator()(const ColourBufferRange &r) const {
            ColourVector tmp;
            auto it = r.begin(),
                 end = r.end();
            while (it != end) {
                do {
                    auto c = *it.colour;
                    ++it;
                    if (c != msk(c)) // only count annotations
                        tmp.push_back(c);
                } while ((it != end) && (*it.boundary == 0));
                if (tmp.size() > 1) {
                    m->increment(tmp);
                }
                tmp.clear();
            }
        }

        void CountAllLists::operator()(const dist::ColourBufferRange &r) const {
            ColourVector tmp;
            auto it = r.begin(),
                 end = r.end();
            while (it != end) {
                do {
                    auto c = *it.colour;
                    tmp.push_back(c);
                    ++it;
                } while ((it != end) && (*it.boundary == 0));
                m->increment(tmp);
                tmp.clear();
            }
        }

        RankMap countPatterns(graph::TraversableGraph const &g, size_t grainsize) {
            RankMap map(g.colourBitWidth);
            ApplyAcrossColours(g.colourBuffers, grainsize, CountColourLists(&map, g.getBaseColourBitWidth()));
            return map;
        }

        RankMap countRaw(graph::TraversableGraph const &g, size_t grainsize) {
            RankMap map(g.colourBitWidth);
            ApplyAcrossColours(g.colourBuffers, grainsize, CountAllLists(&map));
            return map;
        }

        RetrieveColouredKmers::RetrieveColouredKmers(const graph::TraversableGraph *g_, const RankMap *m_, uint8_t w_, std::ostream *os_)
            : g(g_)
            , qry(m_)
            , msk(w_)
            , os(os_)
        {}

        void RetrieveColouredKmers::operator()(const ColourBufferRange &r) const {
            ColourVector tmp;
            size_t e = 0;
            auto it = r.begin(),
                 end = r.end();
            while (it != end) {
                do {
                    auto c = *it.colour;
                    ++it;
                    if (c != msk(c)) // only count annotations
                        tmp.push_back(c);
                } while ((it != end) && (*it.boundary == 0));
                if (tmp.size() > 1) {
                    if (auto q = qry->find(tmp); q != qry->end()) {
                        size_t edgeIndex = r.edgeBegin() + e;
                        if (edgeIndex >= g->edge_count())
                            LOG(ERROR) << "requested edge " << edgeIndex << " [e=" << e << ",offset=" << r.edgeBegin() << "] out of " << g->edge_count() << " edges";
                        printKmer(edgeIndex, q->second);
                    }
                }
                tmp.clear();
                ++e;
            }
        }

        void RetrieveColouredKmers::printKmer(size_t edgeIndex, uint32_t value) const
        { std::osyncstream(*os) << toString(g->kmer(edgeIndex)) << '\t' << value << '\n'; }

        void printColours(graph::TraversableGraph const &g, const RankMap &qry, std::ostream &os) {
            ApplyAcrossColours(g.colourBuffers, 1u, RetrieveColouredKmers(&g, &qry, g.getBaseColourBitWidth(), &os));
        }

    } // namespace cooc

} // namespace dist
