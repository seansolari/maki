#include <math.h>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <gsl/gsl_cdf.h>
#include <glog/logging.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_invoke.h>
#include <maki/classify.hpp>

namespace bmp = boost::multiprecision;

maki_tree::Tree loadAndReduce(const std::string &nwk_file, const std::unordered_map<size_t, GenomeToken> &genomes) {
    LOG(INFO) << "parsing Newick taxonomy from " << nwk_file;
    maki_tree::Tree tree;
    tree.load(nwk_file.c_str());

    LOG(INFO) << "reducing tree to selected genomes";
    size_t missing = 0;
    std::unordered_set<std::string> currentGenomes;
    for (const auto &[_, tk] : genomes) {
        if (tree.id(tk.name) != tree.end())
            currentGenomes.insert(tk.name);
        else ++missing;
    }
    if (missing > 0)
        LOG(ERROR) << missing << " genomes could not be found in provided phylogeny";
    return tree.reduce(currentGenomes);
}

namespace classify {

    int64_t PartialInterleaving::refFindBetween(int64_t l, int64_t r, uint8_t c) const {
        auto lrank = ref.W.rank(l, c),
             rrank = ref.W.rank(r, c);
        if (lrank < rrank) {
            return ref.W.select(lrank + 1, c);
        } else {
            c ^= 0b1000;
            lrank = ref.W.rank(l, c);
            rrank = ref.W.rank(r, c);
            return (lrank < rrank)
                ? ref.W.select(lrank + 1, c)
                : -1;
        }
    }

    int64_t PartialInterleaving::getPrevQryBlock(int64_t zPos) const {
        auto zRank = qryBpRank(zPos);
        if (zRank)
            return qryBpSelect(zRank);
        else return 0;
    }

    int64_t PartialInterleaving::getPrevRefBlock(int64_t zPos) const {
        auto zRank = refBpRank(zPos);
        if (zRank)
            return refBpSelect(zRank);
        else return 0;
    }

    Input PartialInterleaving::from(int64_t eqry) {
        return Input {
            0 /* dummy value */,
            static_cast<int64_t>(qry.C[1] + qry.W.rank(eqry, 0b1001)),
            static_cast<int64_t>(qry.C[2] + qry.W.rank(eqry, 0b1010)),
            static_cast<int64_t>(qry.C[3] + qry.W.rank(eqry, 0b1011)),
            static_cast<int64_t>(qry.C[4] + qry.W.rank(eqry, 0b1100))
        };
    }

    void PartialInterleaving::initPrevBlocks(int64_t zqry,
                                      int64_t &prevBlock,
                                      int64_t &prevRefBlockStart,
                                      Blocks &blocks,
                                      Output &O)
    {
        int64_t prevQryBlock = getPrevQryBlock(zqry+1),
                prevRefBlock = getPrevRefBlock(zqry+1);

        if ((prevQryBlock == zqry) || (prevRefBlock == zqry))
            return;
        else if (prevRefBlock > prevQryBlock) { // in block started by reference
            prevBlock = prevRefBlock;
            prevRefBlockStart = refEdge(refBpPos[prevBlock]);
        } else {                                // in block started by query
            prevBlock = prevQryBlock;
            prevRefBlockStart = -1;
        }

        auto eref = refEdge(previous()[zqry]);
        for (uint64_t c = 1; c < 5; ++c) {
            auto cRankStart = qry.W.rank(qryEdge(prevBlock), c|0b1000),
                 cRankEnd   = qry.W.rank(qryEdge(zqry),      c|0b1000);
            if (cRankStart < cRankEnd) {
                blocks[c] = prevBlock;
                O[c] = rC[c] + ref.W.rank(eref, c|0b1000);
            }
        }
    }

    // `zqry` is the start position of the interleaving, meaning it is the index of a node.
    // Note it is a 0-based number, while select is a 1-based operation, i.e. selecting the
    // first item must call `select(1)`.
    void PartialInterleaving::interleaveRange(int64_t zqry, int64_t zqryEnd) {
        Vector  &Zp = previous(),
                &Z = current();
        int64_t eqry = qryEdge(zqry),       // current edge in query graph
                block = -1,                 // current block
                refBlockStart = -1;
        Input   I = from(eqry);             // output position for each character in interleaving
        Blocks  blockId = {-1,-1,-1,-1,-1}; // block from which each character was last written
        Output  O;                          // cache current position for each character in reference graph
        initPrevBlocks(zqry, block, refBlockStart, blockId, O);

        Accessor a(locks);
        while (zqry < zqryEnd) {
            if (qryBp[zqry]) {
                block = zqry;
                refBlockStart = -1;
            } else if (refBp[zqry]) {
                block = zqry;
                refBlockStart = refEdge(refBpPos[zqry]);
            }

            auto eref = refEdge(Zp[zqry++]);
            do
            {
                uint8_t edge = qry.W[eqry],
                        c = edge & 0b0111;
                if ((edge & 0b1000) && c) {
                    int64_t outz = I[c]++;
                    a.access(outz);
                    if (blockId[c] != block) {
                        blockId[c] = block;
                        auto cRankAtEnd = ref.W.rank(eref, c|0b1000);
                        if (auto cRankAtStart = (refBlockStart == -1)
                                ? cRankAtEnd
                                : ref.W.rank(refBlockStart, c|0b1000);
                            cRankAtStart < cRankAtEnd)
                        {
                            refB[outz] = 1;
                            refBPos[outz] = rC[c] + cRankAtStart;
                        } else {
                            qryB[outz] = 1;
                        }
                        O[c] = rC[c] + cRankAtEnd;
                    }
                    Z[outz] = O[c];
                }
            } while (qry.last[eqry++] == 0);
        }
    }

    void PartialInterleaving::nextIteration() {
        oneapi::tbb::parallel_invoke(
            [&]{ qryBp |= qryB; },
            [&]{ refBp |= refB; },
            [&]{
                oneapi::tbb::parallel_for(
                    (size_t)0, nQueries,
                    [&](size_t const &i)->void {
                        int64_t newv = refBPos[i];
                        if (newv) {
                            refBpPos[i] = newv;
                            refBPos[i] = 0;
                        }
                    }
                );
            }
        );

        oneapi::tbb::parallel_invoke(
            [&]{ sdsl::util::set_to_value(qryB, 0);             },
            [&]{ sdsl::util::set_to_value(refB, 0);             },
            [&]{ sdsl::util::init_support(qryBpRank, &qryBp);    },
            [&]{ sdsl::util::init_support(qryBpSelect, &qryBp);  },
            [&]{ sdsl::util::init_support(refBpRank, &refBp);    },
            [&]{ sdsl::util::init_support(refBpSelect, &refBp);  }
        );

        ++h;
    }

    void PartialInterleaving::doOneInterleave() {
        oneapi::tbb::parallel_for(
            oneapi::tbb::blocked_range<size_t>(0, nQueries, grainsize),
            [&](oneapi::tbb::blocked_range<size_t> const &r)->void {
                interleaveRange(r.begin(), r.end());
            }
        );
        nextIteration();
    }

    void PartialInterleaving::interleave() {
        while (h < ref.k)
            doOneInterleave();
    }

    void PartialInterleaving::fetchRange(int64_t zqry, int64_t zqryEnd, ResultVector &out) {
        Vector &Z = previous();
        int64_t qryIt = qryEdge(zqry),
                qryEnd;
        while (zqry < zqryEnd) {
            qryEnd = qryEdge(zqry + 1);
            if (qryBp[zqry] == 0) {
                int64_t refBegin = refEdge(Z[zqry] - 1),
                        refEnd = refEdge(Z[zqry]);
                while (qryIt < qryEnd) {
                    uint8_t c = qry.W[qryIt];
                    out[qryIt] = refFindBetween(refBegin, refEnd, c);
                    ++qryIt;
                }
            } else {
                qryIt = qryEnd;
            }
            ++zqry;
        }
    }

    ResultVector PartialInterleaving::fetch() {
        ResultVector result(qry.W.size(), -1);
        oneapi::tbb::parallel_for(
            oneapi::tbb::blocked_range<size_t>(0, nQueries, grainsize),
            [&](oneapi::tbb::blocked_range<size_t> const &r)->void {
                fetchRange(r.begin(), r.end(), result);
            }
        );
        return result;
    }

    PartialInterleaving::PartialInterleaving(graph::TraversableGraph const &qry_, graph::TraversableGraph const &ref_, size_t grainsize_)
        : qry(qry_)
        , ref(ref_)
        , grainsize(grainsize_)
        , nQueries(qry.node_count())
        , rC({
            /* <N */ static_cast<int64_t>(ref.C[0]),
            /* <A */ static_cast<int64_t>(ref.C[1]),
            /* <C */ static_cast<int64_t>(ref.C[2]),
            /* <G */ static_cast<int64_t>(ref.C[3]),
            /* <T */ static_cast<int64_t>(ref.C[4]),
            /* <. */ static_cast<int64_t>(ref.node_count())
        })
        , data({ makeVector(nQueries), makeVector(nQueries) })
        , qryB(nQueries, 1024)
        , qryBp(nQueries, 1024)
        , refB(nQueries, 1024)
        , refBp(nQueries, 1024)
        , locks(nQueries, LOCK_SAMPLE_RATE)
        , refBPos(makeVector(nQueries))
        , refBpPos(makeVector(nQueries))
        , qryBpRank()
        , refBpRank()
        , qryBpSelect()
        , refBpSelect()
        , h(0)
    {
        assert(qry.k == ref.k);

        // h=1 interleaving

        std::array<size_t,6> qC = {
            /* <N */ qry.C[0],
            /* <A */ qry.C[1],
            /* <C */ qry.C[2],
            /* <G */ qry.C[3],
            /* <T */ qry.C[4],
            /* <. */ nQueries
        };

        Vector &Z = current(),
               &Zp = previous();
        size_t z = 0,
               zend;
        for (size_t c = 0; c < 5; ++c) {
            size_t qryNodes = qC[c+1]-qC[c],
                   refNodes = rC[c+1]-rC[c];
            if (refNodes) {
                refB[z] = 1;
                refBPos[z] = rC[c];
            }
            else if (qryNodes)
                qryB[z] = 1;

            zend = z + qryNodes;
            while (z < zend) {
                if (c == 0) Zp[z] = rC[c+1];
                Z[z++] = rC[c+1];
            }
        }

        // update rank-select structures

        nextIteration();
    }

    ResultVector classify(const graph::TraversableGraph &qry, const graph::TraversableGraph &ref, size_t grainsize) {
        PartialInterleaving buffers(qry, ref, grainsize);
        buffers.interleave();
        return buffers.fetch();
    }

    auto partitionResults(ResultVector &r, const graph::ColourBufferRegistry &buffers) {
        // (1) check end condition
        if (buffers.rbegin()->first <= static_cast<size_t>(r.back()))
            throw std::runtime_error("partitioning is impossible, out of bounds edge request");
        // start partitioning
        std::unordered_map<size_t,std::pair<size_t,size_t>> partitions;
        size_t lhs = 0, lhsNext = 0;
        auto rhs = buffers.cbegin();
        while (lhs < r.size()) {
            // find the first buffer whose end contains the requested node,
            // valid bounds on `rhs` is guaranteed by (1)
            while (rhs->first <= static_cast<size_t>(r[lhs]))
                ++rhs;
            // find values contained within this buffer
            do ++lhsNext;
            while ((lhsNext < r.size()) && (static_cast<size_t>(r[lhsNext]) < rhs->first));
            partitions.try_emplace(rhs->first, lhs, lhsNext);
            lhs = lhsNext;
        }
        return partitions;
    }

    dist::cooc::RankMap retrieveColours(ResultVector &&r, const graph::TraversableGraph &ref) {
        // filter null matches
        r.erase(std::remove_if(r.begin(), r.end(), [](int64_t x) { return x == -1; }),
                r.end());
        // partition results by buffer ranges
        auto chunks = partitionResults(r, ref.colourBuffers);
        // retrieve required buffers
        dist::cooc::RankMap arr(ref.getColourBitWidth());
        std::vector<std::pair<const size_t, graph::BufferIdPair>> buffers;
        for (const auto &[bufferId, _] : chunks) buffers.push_back(*ref.colourBuffers.getKey(bufferId));
        dist::ApplyAcrossColours(
            ref.colourBuffers,
            1u,
            [&](const graph::TraversableGraph::ColourBufferRange &rng)->void {
                dist::ColourVector key;
                auto [begin, end] = chunks[rng.edgeEnd()];
                while (begin < end) {
                    // collect colour associated with edge
                    auto col = rng.selectEdge(r[begin++]+1u);
                    do {
                        key.push_back(*col.colour);
                        ++col;
                    } while ((col != rng.end()) && (*col.boundary == 0));
                    // increment count for colour
                    arr.increment(key);
                    key.clear();
                }
            },
            buffers
        );
        return arr;
    }


} // namespace classify

namespace stats
{

    namespace detail
    {
        
        // Normal p-value -> Z-score
        double normalPPF(double p)
        { return gsl_cdf_ugaussian_Pinv(p); }

        // Normal Z-score -> p-value
        double normalCDF(double z)
        { return gsl_cdf_ugaussian_P(z); }

        double r1ToQ(double r1, double k) {
            bmp::cpp_bin_float_50
                _r1 = r1,
                _q = 1ULL - bmp::pow(1ULL-_r1, k);
            double q = _q.convert_to<double>();
            return q;
        }

        double qToR1(double q, double k) {
            bmp::cpp_bin_float_50
                _q = q,
                _r1 = 1ULL - bmp::pow(1ULL-_q, 1.0 / k);
            double r1 = _r1.convert_to<double>();
            return r1;
        }

        double varNMut(double _L, double _k, double r1)
        { return r1 == 0.0 ? 0.0 : varNMut(_L, _k, r1, r1ToQ(r1, _k)); }

        double varNMut(double _L, double _k, double r1, double q) {
            bmp::cpp_bin_float_50 _r1 = r1, _q = q;

            bmp::cpp_bin_float_50 varN =
                _L * (1.0 - _q) * (_q * (2.0 * _k + (2.0 / _r1) - 1.0) - 2.0 * _k)
                    + _k * (_k - 1.0) * bmp::pow(1.0 - _q, 2.0)
                    + (2.0 * (1.0 - _q) / bmp::pow(_r1, 2.0)) * ((1.0 + (_k - 1.0) * (1.0 - _q)) * _r1 - _q);
            assert(varN >= 0.0);
            return varN.convert_to<double>();
        }

        double nLow(double _L, double _k, double _q, double _z)
        { return _L * _q - _z * sqrt(varNMut(_L, _k, qToR1(_q, _k), _q)); }

        double nHigh(double _L, double _k, double _q, double _z)
        { return _L * _q + _z * sqrt(varNMut(_L, _k, qToR1(_q, _k), _q)); }

        double nLowDerivative(double _L, double _k, double _q, double _a) {
            double _r1 = qToR1(_q, _k), _z = normalPPF(1.0 - _a / 2.0);
            
            double
                _middle = _r1 == 1.0 ? _k * _L : (_k * _L * (1.0 - _q)) / (1.0 - _r1),
                _offset =
                    ((1.0-_q)*(4.0*_q+_r1*(-6.0+2.0*_k+2.0*_r1-2.0*_k*_r1+_L*(-2.0-2.0*(-1.0+_k)*_r1+_k*std::pow(_r1,2.0))+2.0*(1.0-_q)*(3.0-3.0*_k*(1.0-_r1)-_r1-std::pow(_k,3.0)*std::pow(_r1,2.0)+std::pow(_k,2.0)*(-2.0*_r1+std::pow(_r1,2.0))+_L*(1.0-_r1+_k*(2.0-_r1)*_r1+2.0*std::pow(_k,2.0)*std::pow(_r1,2.0))))))
                    / (2.0*(1.0-_r1)*std::pow(_r1,2.0)*sqrt((1.0-_q)*(-2.0*(1.0-_r1)+_L*(2.0-_r1)*_r1)+std::pow(1.0-_q,2.0)*(2.0*(1.0-_r1)+_k*(2.0-_r1)*_r1+std::pow(_k,2.0)*std::pow(_r1,2.0)+_L*(-2.0*_r1+std::pow(_r1,2.0)-2.0*_k*std::pow(_r1,2.0)))));

            return _middle - _z * _offset;
        }

        double nHighDerivative(double _L, double _k, double _q, double _a) {
            double _r1 = qToR1(_q, _k), _z = normalPPF(1.0 - _a / 2.0);
            double
                _middle = _r1 == 1.0 ? _k * _L : (_k * _L * (1.0 - _q)) / (1.0 - _r1),
                _offset =
                    ((1.0-_q)*(4.0*_q+_r1*(-6.0+2.0*_k+2.0*_r1-2.0*_k*_r1+_L*(-2.0-2.0*(-1.0+_k)*_r1+_k*std::pow(_r1,2.0))+2.0*(1.0-_q)*(3.0-3.0*_k*(1.0-_r1)-_r1-std::pow(_k,3.0)*std::pow(_r1,2.0)+std::pow(_k,2.0)*(-2.0*_r1+std::pow(_r1,2.0))+_L*(1.0-_r1+_k*(2.0-_r1)*_r1+2.0*std::pow(_k,2.0)*std::pow(_r1,2.0))))))
                    / (2.0*(1.0-_r1)*std::pow(_r1,2.0)*sqrt((1.0-_q)*(-2.0*(1.0-_r1)+_L*(2.0-_r1)*_r1)+std::pow(1.0-_q,2.0)*(2.0*(1.0-_r1)+_k*(2.0-_r1)*_r1+std::pow(_k,2.0)*std::pow(_r1,2.0)+_L*(-2.0*_r1+std::pow(_r1,2.0)-2.0*_k*std::pow(_r1,2.0)))));

            return _middle + _z * _offset;
        }

        /**
         * q_for_n_mutated_high-- find q s.t. nHigh(q) == nMut
         *
         * Note: nMut==0 is a special case. When q=0 the formula for variance has a zero
         * in the denominator and thus fails to compute. However, the limit of that
         * formula as q goes to zero is zero (and in fact, it is easy to see that the
         * variance is truly zero when q=0). This means the formulas for nLow and nHigh,
         * e.g. L*q-z*sqrt(varN), are zero when q=0. Thus if nMut=0, 0 is the q for
         * which nHigh(q) == nMut.
         *     
         * nMut==L is another special case. There are two solutions in this case, one
         * of which is q=1. We are interested in the other solutions
         */
        double qFromNMutHigh(double _L, double _k, double _nMut, double _z, bool check) {
            if (_nMut == 0.0)
                return 0.0;

            NHighOptim F = { _L, _k, _z, _nMut };

            double qLeft = 1.0E-5,
                   qRight = _nMut < _L ? 1.0 : 1.0 - 1.0E-5;
            for (size_t attempt = 0; attempt < 10; ++attempt)
            {
                if (F(qLeft) >= 0) qLeft /= 2.0;
                else break;
            }
            
            // at this point,
            //    n_high(L,k,qLeft,z)  - nMut < 0
            //    n_high(L,k,qRight,z) - nMut > 0
            // so we can use the Brent's method to find the solution in the bracketed
            // interval
            double qSoln = brent::zero(qLeft, qRight, /*tol*/ 2.0E-12, F);

            if (check && (qSoln != 1.0)) {
                // limit(dNHigh) as q->1 appears to be non-negative
                double dNHigh = nHighDerivative(_L, _k, qSoln, 2.0 * (1.0 - normalCDF(_z)));
                if (dNHigh <= 0.0)
                    LOG(WARNING)
                        << "solution of nHigh(q)=" << static_cast<size_t>(_nMut)
                        << " (for L=" << static_cast<size_t>(_L)
                        << ", k=" << static_cast<size_t>(_k)
                        << ") fails derivative test...d(nHigh)/dr at q=" << qSoln
                        << " is " << dNHigh;
            }

            return qSoln;
        }

        /**
         * q_for_n_mutated_low-- find q s.t. nLow(q) == nMut
         *
         * Note: nMut==L is a special case. When q=1 all kmers are mutated and thus
         * var(nMut)=0. The formula for nLow then simplifies as
         *    nLow = Lq - z*sqrt(varN) = L*q = L
         * Moreover, if q<1 then varN>0 and nLow < Lq < L.  Thus if nMut=L, 1 is the q
         * for which nLow(q) == nMut.
         */
        double qFromNMutLow(double _L, double _k, double _nMut, double _z, bool check) {
            if (_nMut == _L)
                return 1.0;

            NLowOptim F = { _L, _k, _z, _nMut };

            double qLeft = 1.0E-5,
                   qRight = 1.0L;
            for (size_t attempt = 0; attempt < 10; ++attempt)
            {
                if (F(qLeft) >= 0) qLeft /= 2.0;
                else break;
            }
            
            // at this point,
            //    n_low(L,k,qLeft,z)  - nMut < 0
            //    n_low(L,k,qRight,z) - nMut > 0
            // so we can use the Brent's method to find the solution in the bracketed
            // interval
            double qSoln = brent::zero(qLeft, qRight, /*tol*/ 2.0E-12, F);

            if (check) {
                // limit(dNHigh) as q->1 appears to be non-negative
                double dNLow = nLowDerivative(_L, _k, qSoln, 2.0 * (1.0 - normalCDF(_z)));
                if (dNLow <= 0.0)
                    LOG(WARNING)
                        << "solution of nLow(q)=" << static_cast<size_t>(_nMut)
                        << " (for L=" << static_cast<size_t>(_L)
                        << ", k=" << static_cast<size_t>(_k)
                        << ") fails derivative test...d(nLow)/dr at q=" << qSoln
                        << " is " << dNLow;
            }

            return qSoln;
        }

    } // namespace detail

    // Poisson approximation to the binomial distribution with lambda = np
    void JaccardResult::calculatePoissonPValue() {
        double lambda = std::pow(2.0, std::log2(l1) + std::log2(l2) - static_cast<double>(2 * k));
        pval = gsl_cdf_poisson_Q(shared, lambda);
    }

    void JaccardResult::calculateR1(double z, double minL) {
        double lAve = static_cast<double>(l1 + l2) / 2.0,
               nMut = lAve * (1.0 - jcd) / (1.0 + jcd);
        if (lAve > minL) {
            r1Low = detail::qToR1(detail::qFromNMutHigh(lAve, k, nMut, z), k);
            r1High = detail::qToR1(detail::qFromNMutLow(lAve, k, nMut, z), k);
        }
    }
    
    void JaccardResult::calculatePvals(double z) {
        calculatePoissonPValue();
        calculateR1(z);
    }

    std::ostream& operator<<(std::ostream &os, const JaccardResult &rec) {
        return os
            << rec.gid1 << '\t'
            << rec.gid2 << '\t'
            << rec.shared << '\t'
            << rec.c1 << '\t'
            << rec.c2 << '\t'
            << rec.jcd << '\t'
            << rec.pval << '\t'
            << rec.r1Low << '\t'
            << rec.r1High << '\t'
            << rec.lcaName;
    }

    void JaccardResult::clear() {
        gid1 = gid2 = shared = c1 = c2 = l1 = l2 = k = 0;
        jcd = pval = r1Low = r1High = 0.0;
        lcaName.clear();
    }

    JaccardResultsSampler::JaccardResultsSampler(size_t N_, double z_, std::ostream &os_)
        : stream::OutputHandler<JaccardResult>(os_)
        , N(N_)
        , z(z_)
        , results()
    {}

    JaccardResultsSampler::~JaccardResultsSampler() {
        this->os << std::scientific;
        for (auto &[key, grp] : results) {
            LOG(INFO) << "key " << key << " has " << grp.size() << " results";
            for (auto &rec : grp) {
                rec.value.calculatePvals(z);
                this->os << rec.value << '\n';
            }
        }
    }

    void JaccardResultsSampler::write(JaccardResult &rec) {
        size_t keyEnd = rec.lcaName.find("__");
        if (keyEnd != std::string::npos) {
            auto [p, _] = results.try_emplace(rec.lcaName.substr(0, keyEnd), N);
            p->second.push(rec);
        }
    }

    std::unique_ptr<stream::OutputHandler<JaccardResult>> makeJaccardSampler(size_t take_, double z_, std::ostream &os_) {
        if (!take_)
            return std::make_unique<stream::DirectOutput<JaccardResult>>(z_, os_);
        else
            return std::make_unique<JaccardResultsSampler>(take_, z_, os_);
    }

    std::ostream& operator<<(std::ostream &os, const SamdistResult &rec) {
        return os
            << rec.fid1 << '\t'
            << rec.fid2 << '\t'
            << rec.shared << '\t'
            << rec.c1 << '\t'
            << rec.c2 << '\t'
            << rec.jcd << '\t'
            << rec.lcaName << '\t'
            << ((rec.seed1 != nullptr) ? *rec.seed1 : "") << '\t' 
            << ((rec.seed2 != nullptr) ? *rec.seed2 : "");
    }

    void SamdistResult::clear() {
        fid1 = fid2 = shared = c1 = c2 = l1 = l2 = k = 0;
        jcd = 0.0;
        lcaName.clear();
        seed1 = seed2 = nullptr;
    }

} // namespace stats

