#include <maki/suffix.hpp>

namespace suffix
{

    std::string SmallRollingNuclSeq::toString() const
    {
        std::string nucleotideSequence;
        nucleotideSequence.reserve(_s);
        for (size_t i = 0; i < _s; ++i)
        {
            SmallRollingNuclSeq::window_t rank = (_data >> (2 * i)) & 0b11;
            switch (rank)
            {
            case 0:
                nucleotideSequence.push_back('A');
                break;
            case 1:
                nucleotideSequence.push_back('C');
                break;
            case 2:
                nucleotideSequence.push_back('G');
                break;
            case 3:
                nucleotideSequence.push_back('T');
                break;
            default:
                break;
            }
        }
        return nucleotideSequence;
    }

    void
    SuffixTable::count(Dna4SequenceConstIter it, Dna4SequenceConstIter end)
    {
        SmallRollingNuclSeq suffix { s };

        // initialise suffix

        for (size_t i = 0; i < s; ++i)
            suffix.rollBack(parsing::dna4ToLong(*it++));
        ++(*this)[suffix];
        
        // count rest of sequence

        while (it != end)
        {
            suffix.rollBack(parsing::dna4ToLong(*it++));
            ++(*this)[suffix];
        }
    }

    size_t SuffixTable::maxValue() const {
        return *std::max_element(_data.cbegin(), _data.cend());
    }

    SuffixTable& inplaceAdd(SuffixTable &dest, SuffixTable const &rhs)
    {
        assert(dest.suffixSize() == rhs.suffixSize());
        auto it = dest.begin();
        for (uint64_t v : rhs)
        {
            *it += v;
            ++it;
        }
        return dest;
    }

    SuffixTable accumulate(std::vector<SuffixTable> const &tables)
    {
        assert(!tables.empty());
        SuffixTable result(tables[0].suffixSize());
        std::accumulate(tables.cbegin(), tables.cend(), std::ref(result), inplaceAdd);
        return result;
    }

    SuffixBV::SuffixBV(size_t s)
        : _data((SuffixBV::offset(s + 1) + _w - 1) / _w, (size_t)0)
    {
    }

    size_t SuffixBV::offset(size_t s)
    {
        /*
        Geometric series:
            0 -> 0
            1 -> 4^0
            2 -> 4^0 + 4^1
              ...
            s -> 4^0 + 4^1 + ... + 4^(s-1) = (4^s - 1) / 3;
        */
        size_t x = 1u;
        return ((x << (2 * s)) - 1u) / 3u;
    }

    size_t SuffixBV::bitIndex(SmallRollingNuclSeq const &sfx)
    {
        return SuffixBV::offset(sfx.size()) + (size_t)sfx.data();
    }

    void SuffixBV::insert(SmallRollingNuclSeq const &sfx)
    {
        size_t i = bitIndex(sfx);
        _data[i / _w] |= (size_t)1ul << (i % _w);
    }

    bool SuffixBV::contains(SmallRollingNuclSeq const &sfx) const
    {
        size_t i = bitIndex(sfx);
        return static_cast<size_t>(_data[i / _w]) & ((size_t)1ul << (i % _w));
    }

} // namespace suffix
