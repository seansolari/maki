#pragma once

#include <array>
#include <atomic>
#include <iostream>
#include <memory>
#include <queue>
#include <utility>
#include <vector>
#include <brent/brent.hpp>
#include <sdsl/vectors.hpp>
#include <maki/bitvectors.hpp>
#include <maki/dist.hpp>
#include <maki/graph.hpp>
#include <maki/tree.hpp>
#include <maki/utils.hpp>

namespace pv = bv::threaded;

struct GenomeToken {
    std::string name;
    size_t sequenceLength;
    std::string source;

    GenomeToken(const std::string &name_, size_t length_, const std::string &source_) : name(name_) , sequenceLength(length_), source(source_) {}
    GenomeToken(std::string &&name_, size_t &&length_, std::string &&source_) : name(std::move(name_)) , sequenceLength(std::move(length_)), source(std::move(source_)) {}

    GenomeToken(const GenomeToken&) =default;
    GenomeToken(GenomeToken&&) =default;
    GenomeToken& operator=(const GenomeToken&) =default;
    GenomeToken& operator=(GenomeToken&&) =default;
};

namespace classify {

    using ResultVector = std::vector<int64_t>;

    ResultVector classify(graph::TraversableGraph const &qry, graph::TraversableGraph const &ref, size_t grainsize);

    auto partitionResults(ResultVector &r, const graph::ColourBufferRegistry &buffers);
    dist::cooc::RankMap retrieveColours(ResultVector &&r, graph::TraversableGraph const &ref);

    using Input         = std::array<int64_t,5>;
    using Blocks        = std::array<int64_t,5>;
    using Output        = std::array<int64_t,5>;
    using Vector        = std::vector<int64_t>;
    using Accessor      = pv::SpinRegionManager::Accessor;

    class PartialInterleaving {
        static constexpr size_t LOCK_SAMPLE_RATE = 8 * 1024;

        graph::TraversableGraph const &qry, &ref;
        size_t grainsize, nQueries;
        std::array<int64_t,6> rC;
        std::array<Vector,2> data;
        sdsl::bit_vector qryB, qryBp, refB, refBp;
        pv::SpinRegionManager locks;
        Vector refBPos, refBpPos;
        sdsl::bit_vector::rank_1_type qryBpRank, refBpRank;
        sdsl::bit_vector::select_1_type qryBpSelect, refBpSelect;
        uint8_t h;

    public:
        PartialInterleaving(graph::TraversableGraph const &qry, graph::TraversableGraph const &ref, size_t grainsize_);

        void interleave();
        ResultVector fetch();
    
        inline static auto makeVector(size_t size)
        { return Vector(size); }

        int64_t getPrevQryBlock(int64_t zPos) const;
        int64_t getPrevRefBlock(int64_t zPos) const;

        // get index corresponding to first edge of node at index `nodeIndex` in query graph
        inline int64_t qryEdge(int64_t nodeIndex)
        { return (nodeIndex == 0) ? 0 : qry.lastSelectSupport.select(nodeIndex) + 1; }

        // get index corresponding to first edge of node at index `nodeIndex` in reference graph
        inline int64_t refEdge(int64_t nodeIndex)
        { return (nodeIndex == 0) ? 0 : ref.lastSelectSupport.select(nodeIndex) + 1; }

        int64_t refFindBetween(int64_t l, int64_t r, uint8_t c) const;

        // output buffer for current interleaving iteration
        inline auto& current()
        { return data[h%2]; }

        // input buffer for current interleaving iteration
        inline auto& previous()
        { return data[(h+1)%2]; }

        inline auto& refOverlap() const
        { return refBp; }

        inline auto& qryOverlap() const
        { return qryBp; }

        inline auto& refBlocks() const
        { return refBpPos; }

        Input from(int64_t eqry);
        void initPrevBlocks(int64_t, int64_t&, int64_t&, Blocks&, Output&);
        void interleaveRange(int64_t, int64_t);
        void doOneInterleave();
        void nextIteration();
        void fetchRange(int64_t, int64_t, ResultVector&);
    };

} // namespace classify

namespace stats
{

    namespace detail
    {

        // probability functions
        
        double normalPPF(double p);
        double normalCDF(double z);

        // convert between mutated k-mers and base-pairs

        double r1ToQ(double r1, double k);
        double qToR1(double q, double k);

        // confidence interval for Nmut

        double varNMut(double _L, double _k, double r1);
        double varNMut(double _L, double _k, double r1, double q);

        double nLow(double _L, double _k, double _q, double _z);
        double nHigh(double _L, double _k, double _q, double _z);

        double nLowDerivative(double _L, double _k, double _q, double _a);
        double nHighDerivative(double _L, double _k, double _q, double _a);

        struct OptimBase {
            double _L, _k, _z, _nMut;
            OptimBase(double L, double k, double z, double nMut) : _L(L), _k(k), _z(z), _nMut(nMut) {}
        };

        struct NHighOptim : public brent::func_base
                          , public OptimBase
        {
            using OptimBase::OptimBase;

            virtual double operator()(double q_) override final
            { return nHigh(_L, _k, q_, _z) - _nMut; }
        };

        struct NLowOptim : public brent::func_base
                         , public OptimBase
        {
            using OptimBase::OptimBase;

            virtual double operator()(double q_) override final
            { return nLow(_L, _k, q_, _z) - _nMut; }
        };

        double qFromNMutHigh(double _L, double _k, double _nMut, double _z, bool check = true);
        double qFromNMutLow(double _L, double _k, double _nMut, double _z, bool check = true);

    } // namespace detail

    struct JaccardResult {
        size_t gid1, gid2,      // IDs in newick tree
               l1, l2,
               k,               // k-mer size
               shared, c1, c2;
        double jcd, pval, r1Low, r1High;
        std::string lcaName;

        void calculatePvals(double z);
        friend std::ostream& operator<<(std::ostream&, const JaccardResult&);
        void clear();
    protected:
        void calculatePoissonPValue();
        void calculateR1(double z, double minL = 100.0);
    };

    struct SamdistResult {
        size_t fid1, fid2,      // feature IDs
               l1, l2,
               k,               // k-mer size
               shared, c1, c2;
        double jcd;
        std::string lcaName;
        const std::string *seed1, *seed2;
        friend std::ostream& operator<<(std::ostream&, const SamdistResult&);
        void clear();
    };

} // namespace stats

namespace std {

    template <>
    struct hash<stats::JaccardResult> {
        size_t operator()(const stats::JaccardResult& obj) const {
            size_t h_ = 0;
            h_ ^= hash<size_t>{}(obj.gid1) + 0x9e3779b9;
            h_ ^= hash<size_t>{}(obj.gid2) + 0x9e3779b9 + (h_ << 6) + (h_ >> 2);
            h_ ^= hash<size_t>{}(obj.k) + 0x9e3779b9 + (h_ << 6) + (h_ >> 2);
            return h_;
        }
    };

} // namespace std

namespace stats
{
    
    template <typename T>
    struct HashMaxHeap {
        struct Item {
            T      value;
            size_t hash;

            inline bool operator<(const Item &rhs) const
            { return hash < rhs.hash; }
        };

    protected:
        class _ObjLess {
            std::vector<Item> &_v;

        public:
            _ObjLess(std::vector<Item> &v_) : _v(v_) {}

            inline bool operator()(size_t lhs, size_t rhs) const
            { return _v[lhs] < _v[rhs]; }
        };

        size_t N;   // size of heap
        std::vector<Item> objects;
        std::hash<T> hash;
        std::priority_queue<size_t,std::vector<size_t>,_ObjLess> maxHeap; // max-heap of object indices

    public:
        HashMaxHeap(size_t N_)
            : N(N_)
            , objects()
            , hash()
            , maxHeap(_ObjLess(objects))
        {
            objects.reserve(N);
        }

        void push(const T& t_) {
            size_t t_h = hash(t_);
            if (maxHeap.size() == N) {
                size_t maxObjectIndex = maxHeap.top();
                if (objects[maxObjectIndex].hash > t_h) {
                    maxHeap.pop();
                    objects[maxObjectIndex].value = t_;
                    objects[maxObjectIndex].hash = t_h;
                    maxHeap.push(maxObjectIndex);
                }
            } else {
                objects.emplace_back(t_, t_h);
                maxHeap.push(objects.size()-1);
            }
        }

        inline size_t size() const noexcept
        { return objects.size(); }

        inline std::vector<Item>::iterator begin()
        { return objects.begin(); }

        inline std::vector<Item>::iterator end()
        { return objects.end(); }
    };

    class JaccardResultsSampler : public stream::OutputHandler<JaccardResult> {
        size_t N;
        double z;
        std::unordered_map<std::string,HashMaxHeap<JaccardResult>> results;
    
    public:
        JaccardResultsSampler(size_t N_, double z, std::ostream &os_);
        virtual ~JaccardResultsSampler() override;
        virtual void write(JaccardResult &rec) override;
    };

    std::unique_ptr<stream::OutputHandler<JaccardResult>> makeJaccardSampler(size_t take_, double z_, std::ostream &os_);

} // namespace stats

// specialise direct output for jaccard result

template <>
struct stream::DirectOutput<stats::JaccardResult> : public OutputHandler<stats::JaccardResult> {
protected:
    double z;

public:
    DirectOutput(double z_, std::ostream &os_) : OutputHandler<stats::JaccardResult>(os_), z(z_) {}

    virtual void write(stats::JaccardResult &rec) override {
        rec.calculatePvals(z);
        this->os << rec << '\n';
    }
};
