#pragma once
#include <algorithm>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <oneapi/tbb.h>
#include <oneapi/tbb/concurrent_vector.h>
#include <sdsl/bit_vectors.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <maki/bitvectors.hpp>
#include <maki/colourEncoder.hpp>
#include <maki/fasta.hpp>
#include <maki/suffix.hpp>
#include <maki/matrix.hpp>
#include <maki/radixSort.hpp>
#include <maki/utils.hpp>
#include <maki/maki.h>

namespace pv = bv::threaded;

namespace mers
{

    using parsing::dna4ToDna5;
    using parsing::dna4ToLong;
    using parsing::dna4ToRank;
    using parsing::dna4ToShort;

    namespace detail
    {

        struct _Kmer
        {
            _Kmer() = default;

            _Kmer(Dna4SequenceConstIter begin, Dna4SequenceConstIter end, char e_)
                : sequence(begin, end),
                  edge(e_)
            {
            }

            _Kmer(Dna4Sequence &&sequence, char edge)
                : sequence(std::move(sequence)),
                  edge(edge)
            {
            }

            Dna4Sequence sequence;
            char edge;
        };

        // Raw k-mer pointer, used for testing
        struct _Kmers
        {
            _Kmers(size_t k_) : k(k_), kmers() {}

            void push_back(Dna4Sequence const &seq);

            void sort();

            void uniq();

            void remove_null_terminals();

            void force_remove_terminals();

            void uniq_nodes();

            void print(std::ostream &os) const;

            inline std::ostream &operator<<(std::ostream &os) const
            {
                print(os);
                return os;
            }

            inline size_t size() const noexcept { return kmers.size(); }

            void print_i(std::ostream &os, _Kmer const &p_) const;

            inline auto begin() const { return kmers.begin(); }

            inline auto end() const { return kmers.end(); }

            inline _Kmer const &operator[](size_t i) const noexcept { return kmers[i]; }

            size_t k;
            std::vector<_Kmer> kmers;
        };

        class _CountReference
        {
        public:
            _CountReference(size_t &target) : x(target) {}

            template <class... Args>
            inline _CountReference &operator=([[maybe_unused]] Args &&...args)
            {
                ++x;
                return *this;
            }

        protected:
            size_t &x;
        };

        class _CountIterator
        {
        public:
            _CountIterator(size_t &target) : x(target) {}

            inline _CountReference operator*() { return _CountReference(x); }

            inline void operator++() const { return; }

            inline _CountIterator operator++(int) { return _CountIterator(x); }

        protected:
            size_t &x;
        };

        size_t _pairwise(mers::detail::_Kmers const &, mers::detail::_Kmers const &);

        std::map<std::pair<size_t, size_t>, size_t> _pairwise(std::vector<mers::detail::_Kmers> const &);

        struct _Kmers_LessThan
        {
            bool operator()(_Kmer const &p_, _Kmer const &q_) const;
        };

        struct _Nodes_Eq
        {
            bool operator()(_Kmer const &p_, _Kmer const &q_) const;
        };

        struct _Kmers_Eq
        {
            bool operator()(_Kmer const &p_, _Kmer const &q_) const;
        };

        struct _Kmers_FwdEq
        {
            bool operator()(_Kmer const &l_, _Kmer const &r_) const;
        };

    } // namespace detail

    // estimate number of unique k-mers

    size_t estimateUniqueKmers(size_t totalSequenceLength, uint8_t k);

    /**
     * Byte-wise orientation of DNA
     * ----------------------------
     *
     *  - A byte is made of 8 bits:
     *      ------------------------------------------------------------
     *      | b1 b2 b3 b4 b5 b6 b7 b8 | b9 b10 b11 b12 b13 b14 b15 b16 |
     *      ------------------------------------------------------------
     *
     *  - Keys for RADULS must be presented with orientation of least significant
     *      byte on the left (in this case b1-b8).
     *
     *  - When constructing the bit-packed sequence representation, we pushed
     *      bases consecutively, i.e. the sequence
     *          ( dna1 dna2 dna3 dna4 )
     *      was converted to
     *          ( [dna1b1 dna1b2] [dna2b1 dna2b2] [dna3b1 dna3b2] [dna4b1 dna4b2] ).
     *
     *  - These are pushed into the bit-packed vector in the following order
     *      ---------------------------------------------------------------
     *      | dna4b1 dna4b2 dna3b1 dna3b2 dna2b1 dna2b2 dna1b1 dna1b2 | ...
     *      ---------------------------------------------------------------
     *          such that dna1 is the least significant character of the bit-packed byte.
     *
     *  - Sorting bit-packed k-mers in this orientation, where last dna of the k-mer gives the
     *      most sigificant bits, means that k-mers will be sorted in reverse
     *      lexicographical order.
     *
     *  - Consider a 6-mer, which requires 2 bytes of storage. Their orientation is
     *      ---------------------------------------------
     *      | dna4 dna3 dna2 dna1 | NULL NULL dna6 dna5 |
     *      ---------------------------------------------
     *      | ======== b0 ======= | ======= b1 ======== |
     *      ---------------------------------------------
     *      where b0 is the least significant byte, and b1 is the most significant byte.
     *      This is the orientation required by RADULS.
     */

    /**
     * K-mer Buffer
     */

    // Buffer Key -------------------------------------------------------------------------

    /**
     * writeMerTo
     * ------------
     *
     * Encode base-pairs from input iterator into values pointed
     * at by output iterator. E.g. encodes 10-mer into bytes of
     * the following orientation:
     *
     *      `| 4 3 2 1 | 8 7 6 5 | _ _ 10 9 |`
     *
     */
    template <
        typename OutputPtr,
        typename O = std::remove_pointer_t<OutputPtr>>
    inline OutputPtr writeMerTo(Dna4SequenceConstIter it, OutputPtr out, uint8_t k)
    {
        constexpr size_t bpPerRecord = sizeof(O) * 4; // base-pairs per output record

        // write complete cells
        for (size_t i = 0; i < k / bpPerRecord; ++i)
        {
            O &val = *out++;
            for (size_t j = 0; j < bpPerRecord; ++j)
            {
                val |= dna4ToRank<O>(*it++) << (2 * j);
            }
        }

        // write partial cell
        O &val = *out;
        for (size_t j = 0; j < k % bpPerRecord; ++j)
        {
            val |= dna4ToRank<O>(*it++) << (2 * j);
        }
        return out;
    }

    template <typename OutputPtr,
              typename O = std::remove_pointer_t<OutputPtr>>
        requires std::random_access_iterator<OutputPtr>
    inline void reverseWriteMerTo(Dna4SequenceConstIter it, OutputPtr out, uint8_t k, uint8_t n)
    {
        constexpr size_t bpPerRecord = sizeof(O) * 4; // base-pairs per output record
        for (long i = k-1, end = k-1-n;
             i > end;
             --i)
            *(out + (i / bpPerRecord)) |= dna4ToRank<O>(*it--) << (2 * (i % bpPerRecord));
    }

    // suffix routines ================================================================

    /**
     * Represent nt using the 2-bit encoding `[A, C, G, T] <=> [0, 1, 2, 3]`,
     * also allowing for a terminal edge `$ < A`. Bytes are ordered in vector
     * from least significant to most significant. This is useful for
     * efficiently encoding consecutive k-mers without needing to re-encode
     * the previous letters.
     */
    class ContiguousRollingNuclSeq
    {
    public:
        using ConstUShortSpan = ndim::Span<const uint8_t *>;
        using UShortVector = std::vector<uint8_t>;
        using UShortVectorConstRevIter = typename UShortVector::const_reverse_iterator;

    public:
        // allocate space and initialise edge
        ContiguousRollingNuclSeq(size_t num_nts)
            : _num_nts(num_nts)
            , _num_cells(_mkimem_details::key_size(num_nts))
            , _edge_bit_offset(((num_nts - 1) % 4) * 2)
            , _data(_num_cells, 0)
        {
        }

        ContiguousRollingNuclSeq(size_t _length, Dna4SequenceConstIter _it) : ContiguousRollingNuclSeq(_length)
        {
            writeMerTo(_it, _data.data(), _length);
        }

        ContiguousRollingNuclSeq(Dna4SequenceConstIter begin, Dna4SequenceConstIter end) : ContiguousRollingNuclSeq(std::distance(begin, end))
        {
            writeMerTo(begin, _data.data(), _num_nts);
        }

    public:
        size_t _num_nts;
        uint8_t
            _num_cells,
            _edge_bit_offset;
        UShortVector _data;

    public:
        inline size_t size() const noexcept { return _num_nts; }

        // Returns pointer to least significant byte.
        inline const uint8_t *data() const { return _data.data(); }

        // Returns iterator to most significant byte.
        inline UShortVectorConstRevIter rbegin() const { return _data.rbegin(); }

        inline UShortVectorConstRevIter rend() const { return _data.rend(); }

        std::string toString() const;

        void rollBack(uint8_t nt);

        inline ConstUShortSpan view() const
        {
            return ConstUShortSpan(_data.data(), static_cast<size_t>(_num_cells));
        }
    };

    template <typename T, typename Derived>
    class BaseKmerRandomIterator : public ndim::StridedIteratorBase<T,Derived>
    {
        using BaseIterType = ndim::StridedIteratorBase<T,Derived>;
        using difference_type = std::ptrdiff_t;

    protected:
        uint32_t
            _key_bytes,
            _value_offset;

    public:
        BaseKmerRandomIterator(uint32_t key_bytes, uint32_t value_offset)
            : BaseIterType()
            , _key_bytes(key_bytes)
            , _value_offset(value_offset)
        {
        }

        BaseKmerRandomIterator(T *p, size_t w, uint32_t key_bytes, uint32_t value_offset)
            : BaseIterType(p, w)
            , _key_bytes(key_bytes)
            , _value_offset(value_offset)
        {
        }

        constexpr inline void writeKmer(Dna4SequenceConstIter in, uint8_t kEff) const
            requires(!std::is_const_v<T>)
        {
            writeMerTo(in, this->_data, kEff);
            return;
        }

        constexpr inline void writeKmer(const ContiguousRollingNuclSeq &kmer) const
            requires(!std::is_const_v<T>)
        {
            std::memcpy(this->_data, kmer.data(), _key_bytes);
            return;
        }

        inline uint8_t readEdge() const { return *(this->_data + _value_offset) & 0b111; }

        // most significant byte of key
        inline uint8_t keyMSB() const { return *(this->_data + _key_bytes - 1); }
    };

    template <typename T>
    class KmerBufferRandomAccessIterator : public BaseKmerRandomIterator<T, KmerBufferRandomAccessIterator<T>>
    {
        using BaseIterType = BaseKmerRandomIterator<T, KmerBufferRandomAccessIterator<T>>;

    protected:
        uint32_t _value_bytes;

    public:
        KmerBufferRandomAccessIterator(uint32_t key_bytes, uint32_t value_bytes, uint32_t value_offset)
            : BaseIterType(key_bytes, value_offset)
            , _value_bytes(value_bytes)
        {
        }

        KmerBufferRandomAccessIterator(T *p, size_t w, uint32_t key_bytes, uint32_t value_bytes, uint32_t value_offset)
            : BaseIterType(p, w, key_bytes, value_offset)
            , _value_bytes(value_bytes)
        {
        }

        constexpr inline void writeValue(BufferValue val) const
            requires(!std::is_const_v<T>)
        {
            val.flush(this->_data + this->_value_offset, _value_bytes);
            return;
        }

        inline BufferValue readValue() const { return BufferValue(this->_data + this->_value_offset, _value_bytes); }
    };

    template <typename T>
    class UncolouredKmerBufferRandomAccessIterator : public BaseKmerRandomIterator<T, UncolouredKmerBufferRandomAccessIterator<T>>
    {
        using BaseIterType = BaseKmerRandomIterator<T, UncolouredKmerBufferRandomAccessIterator<T>>;

    public:
        UncolouredKmerBufferRandomAccessIterator(uint32_t key_bytes, uint32_t value_offset)
            : BaseIterType(key_bytes, value_offset)
        {
        }

        UncolouredKmerBufferRandomAccessIterator(T *p, size_t w, uint32_t key_bytes, uint32_t value_offset)
            : BaseIterType(p, w, key_bytes, value_offset)
        {
        }

        constexpr inline void writeValue(uint8_t val) const
            requires(!std::is_const_v<T>)
        {
            *(this->_data + this->_value_offset) = val;
            return;
        }

        inline uint8_t readValue() const { return *(this->_data + this->_value_offset); }
    };

    template <template <typename T> class Iter>
    class BaseKmerBuffer
    {
        using MyMatrix = ndim::Matrix<Iter>;

    public:
        static constexpr uint8_t terminalEdge = 0b00000111;

    public:
        using size_type = typename MyMatrix::size_type;
        using iterator = typename MyMatrix::iterator;
        using const_iterator = typename MyMatrix::const_iterator;
        using sentinel_type = typename MyMatrix::sentinel_type;

    protected:
        BaseKmerBuffer(uint64_t num_records, uint32_t value_bytes, uint8_t _k, uint8_t _k_eff, uint32_t key_padding)
            : _num_records(num_records)
            , _record_capacity(num_records)
            , k(_k)
            , k_eff(_k_eff)
            , _key_bytes(_mkimem_details::key_size(k_eff) + key_padding)
            , _rec_size(_mkimem_details::record_size(_key_bytes, value_bytes))
            , _value_offset(_rec_size - value_bytes)
            , _data(_num_records, _rec_size, sort::ALIGNMENT)
        {
        }

        BaseKmerBuffer(uint32_t value_bytes,
                       uint8_t _k,
                       uint8_t _k_eff,
                       uint32_t key_padding,
                       std::initializer_list<std::initializer_list<uint8_t>> data)
            : _num_records(data.size())
            , _record_capacity(data.size())
            , k(_k)
            , k_eff(_k_eff)
            , _key_bytes(_mkimem_details::key_size(k_eff) + key_padding)
            , _rec_size(_mkimem_details::record_size(_key_bytes, value_bytes))
            , _value_offset(_rec_size - value_bytes)
            , _data(data)
        {
        }

    public:
        BaseKmerBuffer(uint64_t num_records, uint32_t value_bytes, uint8_t _k, uint8_t _k_eff)
            : _num_records(num_records)
            , _record_capacity(num_records)
            , k(_k)
            , k_eff(_k_eff)
            , _key_bytes(_mkimem_details::key_size(k_eff))
            , _rec_size(_mkimem_details::record_size(_key_bytes, value_bytes))
            , _value_offset(_rec_size - value_bytes)
            , _data(_num_records, _rec_size, sort::ALIGNMENT)
        {
        }

        BaseKmerBuffer(uint64_t num_records, uint32_t value_bytes, uint8_t k_eff)
            : BaseKmerBuffer(num_records, value_bytes, k_eff, k_eff)
        {
        }

        BaseKmerBuffer(uint32_t value_bytes,
                       uint8_t k,
                       uint8_t k_eff,
                       std::initializer_list<std::initializer_list<uint8_t>> data)
            : _num_records(data.size())
            , _record_capacity(data.size())
            , k(k)
            , k_eff(k_eff)
            , _key_bytes(_mkimem_details::key_size(k_eff))
            , _rec_size(_mkimem_details::record_size(_key_bytes, value_bytes))
            , _value_offset(_rec_size - value_bytes)
            , _data(data)
        {
        }

        BaseKmerBuffer(const BaseKmerBuffer &obj) = delete;

        BaseKmerBuffer(BaseKmerBuffer &&obj) = default;

    protected:
        uint64_t
            _num_records,
            _record_capacity;
        uint8_t
            k,
            k_eff;
        uint32_t
            _key_bytes,
            _rec_size,
            _value_offset;
        MyMatrix _data;

    public:
        inline void touch(uint32_t threads)
        {
            sort::touchArray(_data.data(), _num_records, _rec_size, threads);
            return;
        }

    public:
        inline size_t size() const noexcept { return _num_records; }

        // Will not reallocate - can only resize up to the initially declared capacity.
        // Sets all values in this range to 0.
        inline void resize(size_t new_size_) noexcept
        {
            assert(new_size_ <= _record_capacity);
            _num_records = new_size_;
            std::fill_n(_data.data(), new_size_ * _rec_size, 0);
        }

        inline size_t keyLeastSigByteOffset() const noexcept { return (sizeof(uint8_t) * 4) % (k_eff - 1); }

        inline uint8_t getK() const noexcept { return k; }

        inline uint8_t getEffK() const noexcept { return k_eff; }

        inline uint32_t recordBytes() const noexcept { return _rec_size; }

        inline uint32_t keyBytes() const noexcept { return _key_bytes; }

        inline uint32_t valueOffset() const noexcept { return _value_offset; }

        template <class... Args>
        inline uint8_t* data(Args&& ...args) { return _data.data(std::forward<Args>(args)...); }

        template <class... Args>
        inline uint8_t const *cdata(Args&& ...args) const { return _data.cdata(std::forward<Args>(args)...); }

        template <class... Args>
        inline auto operator[](Args &&...args) { return _data.operator[](std::forward<Args>(args)...); }

        inline const MyMatrix &memoryview() const { return _data; }
    };
    
    class KmerBuffer : public BaseKmerBuffer<KmerBufferRandomAccessIterator>
    {
        using BaseBufferType = BaseKmerBuffer<KmerBufferRandomAccessIterator>;

    public:
        KmerBuffer(uint64_t num_records, uint32_t value_bytes, uint8_t _k, uint8_t _k_eff)
            : BaseBufferType(num_records, value_bytes, _k, _k_eff)
            , _value_bytes(value_bytes)
        {
        }

        KmerBuffer(uint64_t num_records, uint32_t value_bytes, uint8_t k_eff)
            : BaseBufferType(num_records, value_bytes, k_eff, k_eff)
            , _value_bytes(value_bytes)
        {
        }

        KmerBuffer(uint32_t value_bytes,
                   uint8_t k,
                   uint8_t k_eff,
                   std::initializer_list<std::initializer_list<uint8_t>> data)
            : BaseBufferType(value_bytes, k, k_eff, data)
            , _value_bytes(value_bytes)
        {
        }

    protected:
        uint32_t _value_bytes;

    public:
        inline iterator at(size_type idx) { return _data.at(idx, _key_bytes, _value_bytes, _value_offset); }

        inline const_iterator constAt(size_type idx) const { return _data.at(idx, _key_bytes, _value_bytes, _value_offset); }

        inline const_iterator at(size_type idx) const { return _data.at(idx, _key_bytes, _value_bytes, _value_offset); }

        inline uint32_t valueBytes() const noexcept { return _value_bytes; }

        inline iterator begin() { return _data.begin(_key_bytes, _value_bytes, _value_offset); }

        inline sentinel_type end()
        {
            return _num_records == _record_capacity
                ? _data.end(_key_bytes, _value_bytes, _value_offset)
                : _data.begin(_key_bytes, _value_bytes, _value_offset) + _num_records;
        }

    public:
        iterator insert(iterator it,
                        Dna4SequenceConstIter begin,
                        Dna4SequenceConstIter end,
                        uint64_t colour,
                        bool endIsTerminal);

        iterator insert(iterator it,
                        Dna4SequenceConstIter begin,
                        Dna4SequenceConstIter end,
                        uint64_t colour,
                        bool endIsTerminal,
                        suffix::SmallRollingNuclSeq key);

        template <typename Unwind_, class ...Args>
        void
        insertKmers(std::vector<range_t<Unwind_>> const &data,
                    Unwind_ Apply,
                    std::vector<size_t> const &blocks,
                    Args ...args)
        {
            assert(data.size() == blocks.size());
            oneapi::tbb::parallel_for(
                oneapi::tbb::blocked_range<size_t>(0, blocks.size()),
                [&](oneapi::tbb::blocked_range<size_t> const &r)->void {
                    for (size_t i = r.begin(); i != r.end(); ++i)
                    {
                        iterator it = at(i == 0 ? 0 : blocks[i - 1]);
                        Apply.forEach(data[i], [&](Dna4SequenceConstIter begin, Dna4SequenceConstIter end, uint64_t colour, bool endIsTerminal)->void {
                            it = insert(it, begin, end, colour, endIsTerminal, args...); });
                    }
                }
            );
        }

    public:
        void sort(KmerBuffer *temp, uint32_t threads_);
        void sort(uint32_t threads_);
    };

    class UncolouredKmerBuffer : public BaseKmerBuffer<UncolouredKmerBufferRandomAccessIterator>
    {
        using BaseBufferType = BaseKmerBuffer<UncolouredKmerBufferRandomAccessIterator>;

    public:
        UncolouredKmerBuffer(uint64_t num_records, uint8_t _k, uint8_t _k_eff)
            : BaseBufferType(num_records, 1u, _k, _k_eff)
        {
        }

        UncolouredKmerBuffer(uint64_t num_records, uint8_t k_eff)
            : BaseBufferType(num_records, 1u, k_eff, k_eff)
        {
        }

        UncolouredKmerBuffer(uint8_t k,
                             uint8_t k_eff,
                             std::initializer_list<std::initializer_list<uint8_t>> data)
            : BaseBufferType(1u, k, k_eff, data)
        {
        }

    public:
        inline iterator at(size_type idx) { return _data.at(idx, _key_bytes, _value_offset); }

        inline const_iterator constAt(size_type idx) const { return _data.at(idx, _key_bytes, _value_offset); }

        inline const_iterator at(size_type idx) const { return _data.at(idx, _key_bytes, _value_offset); }

        static inline uint32_t valueBytes() noexcept { return 1u; }

        inline iterator begin() { return _data.begin(_key_bytes, _value_offset); }

        inline sentinel_type end()
        {
            return _num_records == _record_capacity
                ? _data.end(_key_bytes, _value_offset)
                : _data.begin(_key_bytes, _value_offset) + _num_records;
        }

    public:
        iterator insert(iterator it,
                        Dna4SequenceConstIter begin,
                        Dna4SequenceConstIter end,
                        bool endIsTerminal);

        iterator insert(iterator it,
                        Dna4SequenceConstIter begin,
                        Dna4SequenceConstIter end,
                        bool endIsTerminal,
                        suffix::SmallRollingNuclSeq key);

        template <typename Unwind_, class ...Args>
        void
        insertKmers(std::vector<range_t<Unwind_>> const &data,
                    Unwind_ Apply,
                    std::vector<size_t> const &blocks,
                    Args ...args)
        {
            assert(data.size() == blocks.size());
            oneapi::tbb::parallel_for(
                oneapi::tbb::blocked_range<size_t>(0, blocks.size()),
                [&](oneapi::tbb::blocked_range<size_t> const &r)->void {
                    for (size_t i = r.begin(); i != r.end(); ++i)
                    {
                        iterator it = at(i == 0 ? 0 : blocks[i - 1]);
                        Apply.forEach(data[i], [&](Dna4SequenceConstIter begin, Dna4SequenceConstIter end, bool endIsTerminal)->void {
                            it = insert(it, begin, end, endIsTerminal, args...); });
                    }
                }
            );
        }

    public:
        void sort(UncolouredKmerBuffer *temp, uint32_t threads_);
        void sort(uint32_t threads_);
    };

    template <class Buffer>
    std::string toString(Buffer const &buffer)
    {
        return std::accumulate(
            buffer.begin(),
            buffer.end(),
            std::string{},
            [](std::string lhs, ndim::Span<const uint8_t*> rhs) -> std::string
            {
                std::string outer = std::accumulate(
                    rhs.begin(),
                    rhs.end(),
                    std::string{},
                    [](std::string inner, uint8_t x) -> std::string
                    {
                        std::string inner_str = "0b" + std::bitset<8>(x).to_string();
                        return inner.size() == 0 ? inner_str : inner + ", " + inner_str;
                    });
                return lhs.size() == 0 ? outer : lhs + "\n" + outer;
            });
    }

    // k-mer comparison ------------------------------------------------------------

    enum KmerDiffClass : uint8_t
    {
        IS_0   = 0b00u,  // corresponds to pattern `+  +  + ... +`
        IS_K   = 0b01u,  // corresponds to pattern `-  +  + ... +`
        BW_0_K = 0b10u // corresponds to pattern `* ... - ... *`
    };

    struct KmerDiff
    {
        static constexpr uint8_t result_width = 2;

        KmerDiff(uint32_t key_bytes) : _seqWidth(key_bytes) {}

        // `_lhs` and `_rhs` point to the least significant bytes of each k-mer
        KmerDiffClass operator()(const uint8_t *_lhs, const uint8_t *_rhs) const;

    private:
        uint32_t _seqWidth;
    };

    struct MaskedBytesDiff
    {
        MaskedBytesDiff(uint8_t k_)
            : _width(_mkimem_details::key_size(k_)),
              _msb_mask(k_ % 4 == 0 ? 0xFFu : (0xFFu >> (2 * (4 - (k_ % 4)))))
        {
        }

    protected:
        size_t _width;
        uint8_t _msb_mask;

    public:
        /**
         * Compare bytes of two arrays, where `l_` and `r_` point to the least significant
         * bytes. A mask is applied to the most significant byte of `l_`.
         */
        KmerDiffClass operator()(const uint8_t *l_, const uint8_t *r_) const;
    };

    using KmerOverlapVector = sdsl::int_vector<2>;

    template <typename RandomIt_, typename Diff_>
    sdsl::int_vector<2> adjacentDifference(RandomIt_ begin, RandomIt_ end, Diff_ D)
    {
        size_t numElements = end - begin;
        KmerOverlapVector data(numElements, 0);
        if (numElements == 0)
            return data;
        
        data[0] = BW_0_K;

        // move data into thread-safe container

        pv::IntVector<2> parArr(std::move(data), 512);

        tbb::parallel_for(
            tbb::blocked_range<size_t>(1, numElements, KMER_OVERLAP_GRAINSIZE),
            [&](tbb::blocked_range<size_t> const &r)-> void
            {
                pv::IntVector<2>::Accessor outCursor(parArr);
                size_t i = r.begin() - 1;

                RandomIt_ it = begin + i;
                uint8_t const *prev = std::to_address(it),
                              *curr;

                while (++i != r.end())
                {
                    curr = std::to_address(++it);
                    outCursor[i] = D(curr, prev);
                    prev = curr;
                }
            });

        // move data out of container

        KmerOverlapVector result = std::move(parArr).data();
        return result;
    }

    template <KmerDiffClass Xchr, typename RandomIt_, typename Diff_>
    std::pair<sdsl::int_vector<2>,
              sdsl::bit_vector>
    indexedAdjacentDifference(RandomIt_ begin, RandomIt_ end, Diff_ D)
    {
        size_t numElements = end - begin;
        sdsl::int_vector<2> data(numElements, 0);
        sdsl::bit_vector    blocks(numElements, 0);

        if (numElements == 0)
            return std::make_pair(std::move(data), std::move(blocks));

        data[0] = BW_0_K;
        if constexpr (Xchr == BW_0_K)
            blocks[0] = 1;

        // thread-safe views

        static constexpr size_t blocksize = 8 * 512;
        pv::SpinRegionManager lockManager(numElements, blocksize);

        tbb::parallel_for(
            tbb::blocked_range<size_t>(1, numElements, KMER_OVERLAP_GRAINSIZE),
            [&](tbb::blocked_range<size_t> const &r)-> void
            {
                pv::SpinRegionManager::Accessor locks(lockManager);

                size_t i = r.begin() - 1;
                RandomIt_ it = begin + i;
                uint8_t const *prev = std::to_address(it),
                              *curr;

                while (++i != r.end())
                {
                    curr = std::to_address(++it);
                    KmerDiffClass dx = D(curr, prev);

                    // lock is acquired for region containing `i`
                    locks.access(i);
                    data[i] = dx;
                    if (dx == Xchr)
                        blocks[i] = 1;

                    prev = curr;
                }
            });

        // move data out of thread-safe views

        auto result = std::make_pair(std::move(data), std::move(blocks));
        return result;
    }

} // namespace mers

mers::KmerOverlapVector adjacentDifference(mers::KmerBuffer &);

mers::KmerOverlapVector adjacentDifference(mers::UncolouredKmerBuffer &);

template <mers::KmerDiffClass Xchr>
auto indexedAdjacentDifference(mers::KmerBuffer &kmers_) {
    return mers::indexedAdjacentDifference<Xchr>(
        kmers_.begin(),
        kmers_.end(),
        mers::KmerDiff{ kmers_.keyBytes() });
}

template <mers::KmerDiffClass Xchr>
auto indexedAdjacentDifference(mers::UncolouredKmerBuffer &kmers_) {
    return mers::indexedAdjacentDifference<Xchr>(
        kmers_.begin(),
        kmers_.end(),
        mers::KmerDiff{ kmers_.keyBytes() });
}
