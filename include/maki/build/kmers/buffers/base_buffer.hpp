#pragma once
#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <oneapi/tbb/parallel_for.h>
#include <sdsl/int_vector.hpp>

#include "maki/build/utils/bits.hpp"
#include "maki/build/utils/locks.hpp"
#include "maki/maki.h"
#include "sort.hpp"

namespace ndim {

template <typename T> struct ReversePointer {
  using size_type = std::size_t;
  using element_type = T;
  using difference_type = std::ptrdiff_t;
  using value_type = std::remove_cv_t<T>;
  using pointer = T *;
  using reference = T &;
  using iterator_category = std::contiguous_iterator_tag;

public:
  ReversePointer() : _p(nullptr) {}

  ReversePointer(T *p_) : _p(p_) {}

  ReversePointer(const ReversePointer &obj) = default;

  ReversePointer(ReversePointer &&obj) = default;

protected:
  T *_p;

public:
  inline ReversePointer &operator=(ReversePointer &&rhs) {
    _p = rhs._p;
    return *this;
  }

  inline ReversePointer &operator=(const ReversePointer &rhs) {
    _p = rhs._p;
    return *this;
  }

  inline ReversePointer &operator=(T *rhs) {
    _p = rhs;
    return *this;
  }

  inline ReversePointer &operator+=(difference_type rhs) {
    _p -= rhs;
    return *this;
  }

  inline ReversePointer &operator-=(difference_type rhs) {
    _p += rhs;
    return *this;
  }

  inline operator ReversePointer<const T>() const
    requires(!std::is_const_v<T>)
  {
    return ReversePointer((const T *)_p);
  }

  inline reference operator*() const { return *_p; }

  inline pointer operator->() { return _p; }

  inline reference operator[](difference_type rhs) const { return *(_p - rhs); }

  inline ReversePointer &operator++() {
    --_p;
    return *this;
  }

  inline ReversePointer &operator--() {
    ++_p;
    return *this;
  }

  inline ReversePointer operator++(int) {
    ReversePointer tmp(*this);
    --_p;
    return tmp;
  }

  inline ReversePointer operator--(int) {
    ReversePointer tmp(*this);
    ++_p;
    return tmp;
  }

  inline difference_type operator-(const ReversePointer &rhs) const {
    return rhs._p - _p;
  }

  inline ReversePointer operator+(difference_type rhs) const {
    return ReversePointer(_p - rhs);
  }

  inline ReversePointer operator-(difference_type rhs) const {
    return ReversePointer(_p + rhs);
  }

  friend inline ReversePointer operator+(difference_type lhs,
                                         const ReversePointer &rhs) {
    return reverse_pointer(rhs._p - lhs);
  }

  inline bool operator==(const ReversePointer &rhs) const {
    return _p == rhs._p;
  }

  inline bool operator!=(const ReversePointer &rhs) const {
    return _p != rhs._p;
  }

  inline bool operator>(const ReversePointer &rhs) const { return _p < rhs._p; }

  inline bool operator<(const ReversePointer &rhs) const { return _p > rhs._p; }

  inline bool operator>=(const ReversePointer &rhs) const {
    return _p <= rhs._p;
  }

  inline bool operator<=(const ReversePointer &rhs) const {
    return _p >= rhs._p;
  }

  inline bool operator==(const T *rhs) const { return _p == rhs; }

  inline bool operator!=(const T *rhs) const { return _p != rhs; }

  operator T *() const { return _p; }

  inline auto asConst() const noexcept { return ReversePointer<const T>(_p); }
};

template <typename T> struct Span {};

template <typename T> struct Span<T *> {
public:
  using iterator = T *;
  using const_iterator = const T *;
  using reverse_iterator = ReversePointer<T>;
  using element_type = T;
  using value_type = std::remove_cv_t<T>;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using pointer = T *;
  using const_pointer = const T *;
  using reference = T &;

public:
  constexpr Span() = default;

  constexpr Span(T *p_, size_type w_) : _p(p_), _w(w_) {}

  constexpr Span(const Span &obj) = default;

  constexpr Span(Span &&obj) = default;

protected:
  T *_p;
  size_type _w;

public:
  inline Span &operator=(Span &&rhs) {
    _p = rhs._p;
    _w = rhs._w;
    return *this;
  }

  inline Span &operator=(const Span &rhs) {
    _p = rhs._p;
    _w = rhs._w;
    return *this;
  }

  inline Span &operator=(T *rhs) {
    _p = rhs;
    return *this;
  }

  inline bool operator==(const T val) const {
    for (size_type i = 0; i < _w; ++i) {
      if (*(_p + i) != val)
        return false;
    }
    return true;
  }

  inline void write(Span<value_type *> const &rhs) const {
    std::copy(rhs.cbegin(), rhs.cend(), begin());
  }

  inline void write(Span<const value_type *> const &rhs) const {
    std::copy(rhs.cbegin(), rhs.cend(), begin());
  }

  inline bool operator==(Span<value_type *> const &rhs) const {
    return std::equal(cbegin(), cend(), rhs.cbegin());
  }

  inline bool operator==(Span<const value_type *> const &rhs) const {
    return std::equal(cbegin(), cend(), rhs.cbegin());
  }

  inline bool operator!=(Span<value_type *> const &rhs) const {
    return !std::equal(cbegin(), cend(), rhs.cbegin());
  }

  inline bool operator!=(Span<const value_type *> const &rhs) const {
    return !std::equal(cbegin(), cend(), rhs.cbegin());
  }

  inline operator Span<const T *>() const
    requires(!std::is_const_v<T>)
  {
    return Span((const T *)_p, _w);
  }

  inline iterator begin() const noexcept { return _p; }

  inline iterator end() const noexcept { return _p + _w; }

  inline const_iterator cbegin() const noexcept { return _p; }

  inline const_iterator cend() const noexcept { return _p + _w; }

  inline reference front() const { return *_p; }

  inline reference back() const { return *(_p + _w - 1); }

  inline reference at(size_type pos) const {
    if (pos >= _w)
      throw std::out_of_range("index out of range");
    return *(_p + pos);
  }

  inline reference operator[](size_type idx) const { return *(_p + idx); }

  inline pointer data() const noexcept { return _p; }

  inline size_type size() const noexcept { return _w; }

  inline bool empty() const noexcept { return _w == 0; }

  inline Span first(size_type n) const { return Span(_p, n); }

  inline Span last(size_type n) const { return Span(_p + _w - n, n); }

  inline Span<reverse_iterator> rfirst(size_type n) const {
    return Span<reverse_iterator>(reverse_iterator(_p + n - 1), n);
  }

  inline Span<reverse_iterator> rlast(size_type n) const {
    return Span<reverse_iterator>(reverse_iterator(_p + _w - 1), n);
  }

  inline Span<reverse_iterator> rview() const {
    return Span<reverse_iterator>(reverse_iterator(_p + _w - 1), _w);
  }
};

template <typename T> struct Span<ReversePointer<T>> {
public:
  using iterator = ReversePointer<T>;
  using const_iterator = ReversePointer<const T>;
  using reverse_iterator = T *;
  using element_type = T;
  using value_type = std::remove_cv_t<T>;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using pointer = T *;
  using const_pointer = const T *;
  using reference = T &;

public:
  constexpr Span() = default;
  constexpr Span(ReversePointer<T> p_, size_type w_) : _p(p_), _w(w_) {}
  constexpr Span(const Span &obj) = default;
  constexpr Span(Span &&obj) = default;

protected:
  ReversePointer<T> _p;
  size_type _w;

public:
  inline Span &operator=(Span &&rhs) {
    _p = rhs._p;
    _w = rhs._w;
    return *this;
  }

  inline Span &operator=(const Span &rhs) {
    _p = rhs._p;
    _w = rhs._w;
    return *this;
  }

  inline Span &operator=(ReversePointer<T> rhs) {
    _p = rhs;
    return *this;
  }

  inline bool operator==(const T val) const {
    for (size_type i = 0; i < _w; ++i) {
      if (*(_p + i) != val)
        return false;
    }
    return true;
  }

  inline void write(Span<ReversePointer<value_type>> const &rhs) const {
    std::copy(rhs.cbegin(), rhs.cend(), begin());
  }

  inline void write(Span<ReversePointer<const value_type>> const &rhs) const {
    std::copy(rhs.cbegin(), rhs.cend(), begin());
  }

  inline bool operator==(Span<ReversePointer<value_type>> const &rhs) const {
    return std::equal(cbegin(), cend(), rhs.cbegin());
  }

  inline bool
  operator==(Span<ReversePointer<const value_type>> const &rhs) const {
    return std::equal(cbegin(), cend(), rhs.cbegin());
  }

  inline bool operator!=(Span<ReversePointer<value_type>> const &rhs) const {
    return !std::equal(cbegin(), cend(), rhs.cbegin());
  }

  inline bool
  operator!=(Span<ReversePointer<const value_type>> const &rhs) const {
    return !std::equal(cbegin(), cend(), rhs.cbegin());
  }

  inline operator Span<ReversePointer<const T>>() const
    requires(!std::is_const_v<T>)
  {
    return Span((ReversePointer<const T>)_p, _w);
  }

  inline iterator begin() const noexcept { return _p; }

  inline iterator end() const noexcept { return _p + _w; }

  inline const_iterator cbegin() const noexcept { return _p.asConst(); }

  inline const_iterator cend() const noexcept { return _p.asConst() + _w; }

  inline reference front() const { return *_p; }

  inline reference back() const { return *(_p + _w - 1); }

  inline reference at(size_type pos) const {
    if (pos >= _w)
      throw std::out_of_range("index out of range");
    return *(_p + pos);
  }

  inline reference operator[](size_type idx) const { return *(_p + idx); }

  inline pointer data() const noexcept { return _p; }

  inline size_type size() const noexcept { return _w; }

  inline bool empty() const noexcept { return _w == 0; }

  inline Span first(size_type n) const { return Span(_p, n); }

  inline Span last(size_type n) const { return Span(_p + _w - n, n); }

  inline Span<T *> rfirst(size_type n) const {
    return Span<T *>(static_cast<T *>(_p + n - 1), n);
  }

  inline Span<T *> rlast(size_type n) const {
    return Span<T *>(static_cast<T *>(_p + _w - 1), n);
  }

  inline Span<T *> rview() const {
    return Span<T *>(static_cast<T *>(_p + _w - 1), _w);
  }
};

// `T` can have any cv qualification, e.g. `const uint8_t` or `uint8_t`
template <typename T, typename Derived> struct StridedIteratorBase {
public:
  using span_type = Span<T *>;
  using iterator_category = std::random_access_iterator_tag;
  using size_type = std::size_t;
  using element_type = span_type;
  using value_type = span_type;
  using pointer_type = T *;
  using difference_type = std::ptrdiff_t;

public:
  StridedIteratorBase() : _data(nullptr), _w(0) {}

  StridedIteratorBase(T *p_, size_type w_) : _data(p_), _w(w_) {}

protected:
  T *_data;
  size_type _w;

public:
  operator Derived &() { return *static_cast<Derived *>(this); }

  // inline Derived& operator=(const Derived& rhs) { _data = rhs._data; _w =
  // rhs._w; return *this; }

  inline Derived operator+(difference_type rhs) const {
    Derived tmp(*static_cast<const Derived *>(this));
    tmp._data += rhs * _w;
    return tmp;
  }

  inline Derived operator-(difference_type rhs) const {
    Derived tmp(*static_cast<const Derived *>(this));
    tmp._data -= rhs * _w;
    return tmp;
  }

  inline Derived &operator=(T *rhs) {
    _data = rhs;
    return *static_cast<Derived *>(this);
  }

  inline Derived &operator+=(difference_type rhs) {
    _data += _w * rhs;
    return *static_cast<Derived *>(this);
  }

  inline Derived &operator-=(difference_type rhs) {
    _data -= _w * rhs;
    return *static_cast<Derived *>(this);
  }

  inline value_type operator*() const { return span_type(_data, _w); }

  inline value_type operator[](difference_type rhs) const {
    return span_type(_data + (_w * rhs), _w);
  }

  inline pointer_type operator->() const { return _data; }

  inline Derived &operator++() {
    _data += _w;
    return *static_cast<Derived *>(this);
  }

  inline Derived &operator--() {
    _data -= _w;
    return *static_cast<Derived *>(this);
  }

  inline Derived operator++(int) {
    Derived tmp(*static_cast<Derived *>(this));
    _data += _w;
    return tmp;
  }

  inline Derived operator--(int) {
    Derived tmp(*static_cast<Derived *>(this));
    _data -= _w;
    return tmp;
  }

  inline difference_type operator-(const Derived &rhs) const {
    return (_data - rhs._data) / _w;
  }

  inline bool operator==(const Derived &rhs) const {
    return _data == rhs._data;
  }

  inline bool operator!=(const Derived &rhs) const {
    return _data != rhs._data;
  }

  inline bool operator>(const Derived &rhs) const { return _data > rhs._data; }

  inline bool operator<(const Derived &rhs) const { return _data < rhs._data; }

  inline bool operator>=(const Derived &rhs) const {
    return _data >= rhs._data;
  }

  inline bool operator<=(const Derived &rhs) const {
    return _data <= rhs._data;
  }

  inline bool operator==(const T *rhs) const { return _data == rhs; }

  inline bool operator!=(const T *rhs) const { return _data != rhs; }
};

template <typename T>
struct StridedIterator : public StridedIteratorBase<T, StridedIterator<T>> {
  template <class... Args>
  StridedIterator(Args &&...args)
      : StridedIteratorBase<T, StridedIterator<T>>(
            std::forward<Args>(args)...) {}
};

template <template <class...> class IteratorBase = StridedIterator>
class Matrix {
public:
  using value_type = uint8_t;
  using size_type = std::size_t;
  using span_type = Span<value_type *>;
  using const_span_type = Span<const value_type *>;
  using difference_type = std::ptrdiff_t;
  using reference = value_type &;
  using const_reference = const value_type &;
  using pointer = value_type *;
  using const_pointer = const value_type *;
  using iterator = IteratorBase<value_type>;
  using const_iterator = IteratorBase<const value_type>;
  using sentinel_type = IteratorBase<value_type>;

public:
  Matrix(size_type nrows_, size_type ncols_, uint32_t alignment)
      : _rows(nrows_), _cols(ncols_), _size(_rows * _cols),
        _p(new value_type[_size + alignment]()) {
    _data = _p.get();
    while (reinterpret_cast<uintptr_t>(_data) % alignment)
      ++_data;
  }
  Matrix(size_type nrows_, size_type ncols_)
      : _rows(nrows_), _cols(ncols_), _size(_rows * _cols),
        _p(new value_type[_size]()) {
    _data = _p.get();
  }
  Matrix(std::initializer_list<std::initializer_list<value_type>> lists_) {
    _rows = lists_.size();
    if (_rows == 0) {
      _cols = 0;
    } else {
      // establish column size
      auto it = lists_.begin();
      _cols = it->size();
      while (it != lists_.end()) {
        _cols = std::min(_cols, it->size());
        ++it;
      }
      _size = _rows * _cols;
      // declare data
      _p = std::make_unique<value_type[]>(_size);
      _data = _p.get();
      // copy data
      value_type *o = _data;
      for (const auto &l_ : lists_) {
        o = std::copy_n(l_.begin(), _cols, o);
      }
    }
  }

  Matrix(const Matrix &obj) = delete;
  Matrix(Matrix &&obj) = default;

protected:
  size_type _rows, _cols, _size;
  std::unique_ptr<value_type[]> _p;
  value_type *_data;

public:
  inline void swap(Matrix &obj) {
    assert(_rows == obj._rows);
    assert(_cols == obj._cols);
    assert(_size == obj._size);
    _p.swap(obj._p);
    std::swap(_data, obj._data);
  }

  template <class... Args> inline iterator begin(Args &&...args) {
    return iterator(_data, _cols, std::forward<Args>(args)...);
  }
  template <class... Args> inline iterator end(Args &&...args) {
    return iterator(_data + _size, _cols, std::forward<Args>(args)...);
  }
  template <class... Args> inline const_iterator cbegin(Args &&...args) const {
    return const_iterator(_data, _cols, std::forward<Args>(args)...);
  }
  template <class... Args> inline const_iterator cend(Args &&...args) const {
    return const_iterator(_data + _size, _cols, std::forward<Args>(args)...);
  }
  template <class... Args> iterator at(size_type idx, Args &&...args) {
    assert(idx < _rows);
    return iterator(_data + (idx * _cols), _cols, std::forward<Args>(args)...);
  }
  template <class... Args>
  const_iterator at(size_type idx, Args &&...args) const {
    assert(idx < _rows);
    return const_iterator(_data + (idx * _cols), _cols,
                          std::forward<Args>(args)...);
  }
  inline pointer data() { return _data; }
  inline pointer data(size_type idx) { return _data + (idx * _cols); }
  inline const_pointer cdata() const { return _data; }
  inline const_pointer cdata(size_type idx) const {
    return _data + (idx * _cols);
  }
  inline size_type size() const noexcept { return _size; }
  inline size_type rows() const noexcept { return _rows; }
  inline size_type cols() const noexcept { return _cols; }
  inline span_type operator[](size_type idx) {
    return span_type(_data + (idx * _cols), _cols);
  }
  inline const_span_type operator[](size_type idx) const {
    return const_span_type(_data + (idx * _cols), _cols);
  }
};

} // namespace ndim

template <template <typename T> class Iter> class BaseKmerBuffer {
  using MyMatrix = ndim::Matrix<Iter>;

public:
  static constexpr uint8_t terminalEdge = 0b00000111;

public:
  using size_type = typename MyMatrix::size_type;
  using iterator = typename MyMatrix::iterator;
  using const_iterator = typename MyMatrix::const_iterator;
  using sentinel_type = typename MyMatrix::sentinel_type;

protected:
  BaseKmerBuffer(uint64_t num_records, uint32_t value_bytes, uint8_t _k,
                 uint8_t _k_eff, uint32_t key_padding)
      : _num_records(num_records), _record_capacity(num_records), k(_k),
        k_eff(_k_eff), _key_bytes(key_size(k_eff) + key_padding),
        _rec_size(record_size(_key_bytes, value_bytes)),
        _value_offset(_rec_size - value_bytes),
        _data(_num_records, _rec_size, ALIGNMENT) {}

  BaseKmerBuffer(uint32_t value_bytes, uint8_t _k, uint8_t _k_eff,
                 uint32_t key_padding,
                 std::initializer_list<std::initializer_list<uint8_t>> data)
      : _num_records(data.size()), _record_capacity(data.size()), k(_k),
        k_eff(_k_eff), _key_bytes(key_size(k_eff) + key_padding),
        _rec_size(record_size(_key_bytes, value_bytes)),
        _value_offset(_rec_size - value_bytes), _data(data) {}

public:
  BaseKmerBuffer(uint64_t num_records, uint32_t value_bytes, uint8_t _k,
                 uint8_t _k_eff)
      : _num_records(num_records), _record_capacity(num_records), k(_k),
        k_eff(_k_eff), _key_bytes(key_size(k_eff)),
        _rec_size(record_size(_key_bytes, value_bytes)),
        _value_offset(_rec_size - value_bytes),
        _data(_num_records, _rec_size, ALIGNMENT) {}

  BaseKmerBuffer(uint64_t num_records, uint32_t value_bytes, uint8_t k_eff)
      : BaseKmerBuffer(num_records, value_bytes, k_eff, k_eff) {}

  BaseKmerBuffer(uint32_t value_bytes, uint8_t k, uint8_t k_eff,
                 std::initializer_list<std::initializer_list<uint8_t>> data)
      : _num_records(data.size()), _record_capacity(data.size()), k(k),
        k_eff(k_eff), _key_bytes(key_size(k_eff)),
        _rec_size(record_size(_key_bytes, value_bytes)),
        _value_offset(_rec_size - value_bytes), _data(data) {}

  BaseKmerBuffer(const BaseKmerBuffer &obj) = delete;

  BaseKmerBuffer(BaseKmerBuffer &&obj) = default;

protected:
  uint64_t _num_records, _record_capacity;
  uint8_t k, k_eff;
  uint32_t _key_bytes, _rec_size, _value_offset;
  MyMatrix _data;

public:
  inline size_t size() const noexcept { return _num_records; }

  /**
   * Will not reallocate - can only resize up to the initially declared
   * capacity. Sets all values in this range to 0.
   */
  inline void resize(size_t new_size_) noexcept {
    assert(new_size_ <= _record_capacity);
    _num_records = new_size_;
    std::fill_n(_data.data(), new_size_ * _rec_size, 0);
  }

  inline size_t keyLeastSigByteOffset() const noexcept {
    return (sizeof(uint8_t) * 4) % (k_eff - 1);
  }

  inline uint8_t getK() const noexcept { return k; }
  inline uint8_t getEffK() const noexcept { return k_eff; }
  inline uint32_t recordBytes() const noexcept { return _rec_size; }
  inline uint32_t keyBytes() const noexcept { return _key_bytes; }
  inline uint32_t valueOffset() const noexcept { return _value_offset; }
  template <class... Args> inline uint8_t *data(Args &&...args) {
    return _data.data(std::forward<Args>(args)...);
  }
  template <class... Args> inline uint8_t const *cdata(Args &&...args) const {
    return _data.cdata(std::forward<Args>(args)...);
  }
  template <class... Args> inline auto operator[](Args &&...args) {
    return _data.operator[](std::forward<Args>(args)...);
  }
  inline const MyMatrix &memoryview() const { return _data; }
};

enum KmerDiffClass : uint8_t {
  IS_0 = 0b00u,  // corresponds to pattern `+  +  + ... +`
  IS_K = 0b01u,  // corresponds to pattern `-  +  + ... +`
  BW_0_K = 0b10u // corresponds to pattern `* ... - ... *`
};

template <typename It>
concept strides_bytes = requires(It __i, It __j, It &__k, std::size_t di) {
  typename It::difference_type;
  { __i + di } -> std::same_as<It>;
  { ++__i } -> std::same_as<It &>;
  { __i - __j } -> std::convertible_to<std::size_t>;
  { std::to_address(__k) } -> std::convertible_to<const uint8_t *>;
};

template <typename Op>
concept kmer_comparator =
    requires(Op op, const uint8_t *lhs, const uint8_t *rhs) {
      { op(lhs, rhs) } -> std::same_as<KmerDiffClass>;
    };

/**
 * In-place adjacent difference - note that `data` size is not
 * shrunk if over-allocated. First element defaults to `BW_0_K`.
 */
template <strides_bytes RandomIter>
void adjacentDifference(RandomIter begin, RandomIter end,
                        kmer_comparator auto &&diff,
                        sdsl::int_vector<2> &data) {
  // resize vector if necessary
  size_t numElements = end - begin;
  if (numElements == 0)
    return;
  else if (numElements > data.size())
    data.resize(numElements);

  // parallel fill
  data[0] = BW_0_K;
  static constexpr size_t blocksize = 8 * 512;
  LockedRegionManager mgr(numElements, blocksize);

  tbb::parallel_for(
      tbb::blocked_range<size_t>(1, numElements, KMER_OVERLAP_GRAINSIZE),
      [&](const tbb::blocked_range<size_t> &r) -> void {
        LockedRegionManager::Accessor locks(mgr);

        size_t i = r.begin() - 1;
        RandomIter it = begin + i;
        const uint8_t *prev = std::to_address(it);
        const uint8_t *curr;

        while (++i != r.end()) {
          curr = std::to_address(++it);
          KmerDiffClass dx = diff(curr, prev);
          locks.access(i);
          data[i] = dx;
          prev = curr;
        }
      });
}

/**
 * Out-of-place adjacent difference. Size of result reflects size of input
 * range. First element defaults to `BW_0_K`.
 */
template <strides_bytes RandomIter>
sdsl::int_vector<2> adjacentDifference(RandomIter begin, RandomIter end,
                                       kmer_comparator auto &&diff) {
  sdsl::int_vector<2> data;
  adjacentDifference(begin, end, diff, data);
  return data;
}
