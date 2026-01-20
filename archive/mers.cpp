#include <maki/mers.hpp>

#include <indicators/progress_bar.hpp>

namespace mers
{

    namespace detail
    {

        void _Kmers::print(std::ostream &os) const
        {
            for (_Kmer const &p : kmers)
            {
                print_i(os, p);
                os << '\n';
            }
        }

        void _Kmers::print_i(std::ostream &os, _Kmer const &p_) const
        {
            for (size_t i = 0; i < k - p_.sequence.size(); ++i)
                os << '$';
            for (auto const &c : p_.sequence)
                os << seqan3::to_char(c);
            os << ' ' << p_.edge;
        }

        void _Kmers::push_back(Dna4Sequence const &seq)
        {
            auto
                l = seq.cbegin(),
                r = l;
            for (; r < l + k; ++r)
            {
                kmers.emplace_back(l, r, seqan3::to_char(*r));
            }
            while (r < seq.cend())
            {
                kmers.emplace_back(l, r, seqan3::to_char(*r));
                ++l;
                ++r;
            }
            kmers.emplace_back(l, r, '$');
        }

        void _Kmers::sort()
        {
            std::stable_sort(kmers.begin(), kmers.end(), _Kmers_LessThan{});
        }

        void _Kmers::uniq()
        {
            auto
                begin = kmers.begin(),
                end = std::unique(begin, kmers.end(), _Kmers_Eq{});
            kmers.resize(end - begin);
        }

        void _Kmers::uniq_nodes()
        {
            auto
                begin = kmers.begin(),
                end = std::unique(begin, kmers.end(), _Nodes_Eq{});
            kmers.resize(end - begin);
        }

        void _Kmers::remove_null_terminals()
        {
            for (size_t i = 1; i < kmers.size(); ++i)
            {
                if ((kmers[i - 1].sequence == kmers[i].sequence) && (kmers[i - 1].edge == '$'))
                {
                    kmers.erase(kmers.begin() + (i - 1));
                }
            }
        }

        void _Kmers::force_remove_terminals()
        {
            kmers.erase(std::remove_if(kmers.begin(),
                                       kmers.end(),
                                       [&](_Kmer const &v) -> bool
                                       { return (v.sequence.size() < k) || (v.edge == '$'); }),
                        kmers.end());
        }

        size_t
        _pairwise(mers::detail::_Kmers const &lhs,
                  mers::detail::_Kmers const &rhs)
        {
            size_t count = 0u;
            _CountIterator out(count);

            std::set_intersection(lhs.begin(), lhs.end(),
                                  rhs.begin(), rhs.end(),
                                  out,
                                  _Kmers_LessThan{});

            return count;
        }

        std::map<std::pair<size_t, size_t>, size_t>
        _pairwise(std::vector<mers::detail::_Kmers> const &kmers)
        {
            std::map<std::pair<size_t, size_t>, size_t> res;

            for (size_t i = 1; i < kmers.size(); ++i)
            {
                for (size_t j = 0; j < i; ++j)
                {
                    size_t c = _pairwise(kmers[j], kmers[i]);
                    if (c > 0)
                        res[std::make_pair(j + 1, i + 1)] = c;
                }
            }

            return res;
        }

        bool _Kmers_LessThan::operator()(_Kmer const &p_, _Kmer const &q_) const
        {
            auto
                prBegin = std::make_reverse_iterator(p_.sequence.end()),
                prEnd = std::make_reverse_iterator(p_.sequence.begin()),
                qrBegin = std::make_reverse_iterator(q_.sequence.end()),
                qrEnd = std::make_reverse_iterator(q_.sequence.begin());

            if (std::equal(prBegin, prEnd,
                           qrBegin, qrEnd))
            {
                return p_.edge < q_.edge; // Note '$' < 'A' < 'C' ...
            }
            else
            {
                return std::lexicographical_compare(prBegin, prEnd,
                                                    qrBegin, qrEnd);
            }
        }

        bool _Nodes_Eq::operator()(_Kmer const &p_, _Kmer const &q_) const
        {
            return std::equal(p_.sequence.begin(), p_.sequence.end(),
                              q_.sequence.begin(), q_.sequence.end());
        }

        bool _Kmers_Eq::operator()(_Kmer const &p_, _Kmer const &q_) const
        {
            return std::equal(p_.sequence.begin(), p_.sequence.end(),
                              q_.sequence.begin(), q_.sequence.end()) &&
                   (p_.edge == q_.edge);
        }

        bool _Kmers_FwdEq::operator()(_Kmer const &l_, _Kmer const &r_) const
        {
            assert(r_.sequence.size() > 0);
            size_t rangeSize = r_.sequence.size() - 1;
            auto
                l_end = l_.sequence.cend(),
                r_begin = r_.sequence.cbegin();
            return std::equal(l_end - rangeSize, l_end,
                              r_begin, r_begin + rangeSize) &&
                   (l_.edge == seqan3::to_char(r_.sequence.back()));
        }

    } // namespace detail

    size_t estimateUniqueKmers(size_t totalSequenceLength, uint8_t k)
    {
        long double totalSequenceLengthFl = static_cast<long double>(totalSequenceLength);
        long double
            dens = totalSequenceLengthFl / std::pow(4.0l, k),
            res = totalSequenceLengthFl / (1.0l + dens);
        return static_cast<size_t>(std::ceil(res));
    }

    std::string ContiguousRollingNuclSeq::toString() const
    {
        std::string sequence;
        sequence.reserve(_num_nts);
        for (size_t i = 0; i < _num_nts; ++i)
        {
            uint8_t rank = (_data[i / 4] >> ((i % 4) * 2)) & 0b11;
            switch (rank)
            {
            case 0:
                sequence.push_back('A');
                break;
            case 1:
                sequence.push_back('C');
                break;
            case 2:
                sequence.push_back('G');
                break;
            case 3:
                sequence.push_back('T');
                break;
            default:
                break;
            }
        }
        return sequence;
    }

    void ContiguousRollingNuclSeq::rollBack(uint8_t nt)
    {
        // proceed whole bytes

        uint8_t
            i = 0,
            _e;
        for (; i < _num_cells - 1; ++i)
        {
            _e = _data[i + 1] & 0b11;
            _data[i] = (_data[i] >> 2) | (_e << 6);
        }

        // finalise partial byte

        _data[i] = (_data[i] >> 2) | (nt << _edge_bit_offset);
    }

    KmerBuffer::iterator
    KmerBuffer::insert(iterator it,
                       Dna4SequenceConstIter seqIter,
                       Dna4SequenceConstIter endIter,
                       uint64_t colour,
                       bool endIsTerminal)
    {
        ContiguousRollingNuclSeq kmer(k, seqIter);
        seqIter += k;

        uint8_t edge;

        while (seqIter != endIter)
        {
            edge = dna4ToShort(*seqIter++);

            it.writeKmer(kmer);
            it.writeValue(BufferValue(colour, edge));
            ++it;

            kmer.rollBack(edge);
        }

        if (endIsTerminal)
        {
            it.writeKmer(kmer);
            it.writeValue(BufferValue((uint64_t)0u, terminalEdge));
            ++it;
        }

        return it;
    }

    KmerBuffer::iterator
    KmerBuffer::insert(iterator it,
                       Dna4SequenceConstIter seqIter,
                       Dna4SequenceConstIter endIter,
                       uint64_t colour,
                       bool endIsTerminal,
                       suffix::SmallRollingNuclSeq key)
    {
        uint8_t s = key.size();
        seqIter += k - s;

        // initialise
        suffix::SmallRollingNuclSeq suffix(s, seqIter);
        seqIter += s;

        uint64_t edge;

        // fill
        while (seqIter != endIter)
        {
            edge = dna4ToLong(*seqIter);

            if (suffix == key)
            {
                it.writeKmer(seqIter - k, k_eff);
                it.writeValue(BufferValue(colour, edge));
                ++it;
            }

            suffix.rollBack(edge);
            ++seqIter;
        }

        // terminal edge
        if ((suffix == key) && endIsTerminal)
        {
            it.writeKmer(seqIter - k, k_eff);
            it.writeValue(BufferValue((uint64_t)0u, terminalEdge));
            ++it;
        }

        return it;
    }

    void KmerBuffer::sort(KmerBuffer *temp, uint32_t threads_)
    {
        KmerBuffer *result = sort::sort(this, temp, threads_);
        if (result != this)
            _data.swap(result->_data);
    }

    void KmerBuffer::sort(uint32_t threads_)
    {
        KmerBuffer temp(_num_records, _value_bytes, k, k_eff);
        sort(&temp, threads_);
    }

    UncolouredKmerBuffer::iterator
    UncolouredKmerBuffer::insert(iterator it,
                                 Dna4SequenceConstIter seqIter,
                                 Dna4SequenceConstIter endIter,
                                 bool endIsTerminal)
    {
        ContiguousRollingNuclSeq kmer(k, seqIter);
        seqIter += k;

        uint8_t edge;

        while (seqIter != endIter)
        {
            edge = dna4ToShort(*seqIter++);

            it.writeKmer(kmer);
            it.writeValue(BufferValue(edge));
            ++it;

            kmer.rollBack(edge);
        }

        if (endIsTerminal)
        {
            it.writeKmer(kmer);
            it.writeValue(BufferValue(terminalEdge));
            ++it;
        }

        return it;
    }

    UncolouredKmerBuffer::iterator
    UncolouredKmerBuffer::insert(iterator it,
                                 Dna4SequenceConstIter seqIter,
                                 Dna4SequenceConstIter endIter,
                                 bool endIsTerminal,
                                 suffix::SmallRollingNuclSeq key)
    {
        uint8_t s = key.size();
        seqIter += k - s;

        // initialise
        suffix::SmallRollingNuclSeq suffix {s, seqIter};
        seqIter += s;

        uint64_t edge;

        // fill
        while (seqIter != endIter)
        {
            edge = dna4ToLong(*seqIter);

            if (suffix == key)
            {
                it.writeKmer(seqIter - k, k_eff);
                it.writeValue(BufferValue(edge));
                ++it;
            }

            suffix.rollBack(edge);
            ++seqIter;
        }

        // terminal edge
        if ((suffix == key) && endIsTerminal)
        {
            it.writeKmer(seqIter - k, k_eff);
            it.writeValue(BufferValue(terminalEdge));
            ++it;
        }

        return it;
    }

    void UncolouredKmerBuffer::sort(UncolouredKmerBuffer *temp, uint32_t threads_)
    {
        UncolouredKmerBuffer *result = sort::sort(this, temp, threads_);
        if (result != this)
            _data.swap(result->_data);
    }

    void UncolouredKmerBuffer::sort(uint32_t threads_)
    {
        UncolouredKmerBuffer temp(_num_records, k, k_eff);
        sort(&temp, threads_);
    }

    KmerDiffClass
    KmerDiff::operator()(const uint8_t *_lhs, const uint8_t *_rhs) const
    {
        // - check most significant bytes
        if ((_seqWidth > 1) && (!std::equal(_lhs + 1, _lhs + _seqWidth,
                                            _rhs + 1)))
        {
            return BW_0_K;
        }
        else
        {
            // - check least significant byte
            uint8_t
                _lhs_v = *_lhs,
                _rhs_v = *_rhs;

            if (_lhs_v == _rhs_v)
                return IS_0;
            else if ((_lhs_v >> 2) == (_rhs_v >> 2))
                return IS_K;
            else
                return BW_0_K;
        }
    }

    KmerDiffClass
    MaskedBytesDiff::operator()(const uint8_t *l_, const uint8_t *r_) const
    {
        uint8_t msk = _msb_mask;

        for (int i = _width - 1; i >= 0; --i)
        {
            uint8_t
                l_val = *(l_ + i),
                r_val = *(r_ + i);

            l_val &= msk;

            if (l_val != r_val)
            {
                if (i == 0)
                {
                    if ((l_val >> 2) == (r_val >> 2))
                    {
                        return IS_K;
                    }
                }
                return BW_0_K;
            }

            msk = 0xFFu;
        }

        return IS_0;
    }

} // namespace mers

mers::KmerOverlapVector adjacentDifference(mers::KmerBuffer &kmers_)
{
    return mers::adjacentDifference(kmers_.begin(),
                                    kmers_.end(),
                                    mers::KmerDiff{ kmers_.keyBytes() });
}

mers::KmerOverlapVector adjacentDifference(mers::UncolouredKmerBuffer &kmers_)
{
    return mers::adjacentDifference(kmers_.begin(),
                                    kmers_.end(),
                                    mers::KmerDiff{ kmers_.keyBytes() });
}
