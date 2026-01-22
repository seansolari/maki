
#pragma once
#include <cstdint>
#include <vector>
#include <string_view>
#include <syncstream>

#include <gtl/phmap.hpp>
#include <oneapi/tbb/tbb_allocator.h>
#include <oneapi/tbb/parallel_for.h>

/**
 * Map sets of values to integers. Groups of `Int`s form the keys, which are byte-packed for compression.
 * These tuples are mapped to integers of type `Val`.
 */
template <typename Int, typename Val>
class TupleMap
{
  using Key = std::vector<uint8_t, oneapi::tbb::tbb_allocator<uint8_t>>;

  struct ByteVectorHash
  {
    using is_transparent = void;
    size_t operator()(const Key &v) const noexcept
    {
      std::string_view sv(reinterpret_cast<const char *>(v.data()), v.size());
      return std::hash<std::string_view>{}(sv);
    }
    size_t operator()(std::string_view sv) const noexcept
    {
      return std::hash<std::string_view>{}(sv);
    }
  };

  struct ByteVectorEq
  {
    using is_transparent = void;
    bool operator()(const Key &a, const Key &b) const noexcept
    {
      return a == b;
    }
    bool operator()(const Key &a, std::string_view b) const noexcept
    {
      std::string_view sva(reinterpret_cast<const char *>(a.data()), a.size());
      return sva == b;
    }
  };

  using Map = gtl::parallel_flat_hash_map<
      Key,            // Key
      Val,            // Mapped value
      ByteVectorHash, // Hash functor
      ByteVectorEq,   // Equality comparator
      oneapi::tbb::tbb_allocator<std::pair<const Key, Val>>,
      /*N=*/7,   // number of submaps is 2^N => 128
      std::mutex // per-submap mutex for thread safety
      >;

  size_t _bypk; // bytes per key
  Map _data;

protected:
  // copy bytes in LSB-first order
  Key _compress(const Int *ptr_, std::size_t len_) const
  {
    Key out(len_ * _bypk);
    uint8_t *o = out.data();
    for (std::size_t i = 0; i < len_; ++i)
    {
      for (std::size_t b = 0; b < _bypk; ++b)
      {
        *o++ = (ptr_[i] >> (8 * b)) & 0xFFu;
      }
    }
    return out;
  }

  std::vector<Val> _decompress(const Key &bytes) const
  {
    const uint8_t *src = bytes.data();
    const size_t n = bytes.size() / _bypk;
    std::vector<Val> out(n, 0);
    Val *o = out.data();
    for (std::size_t b = 0; b < n; ++b)
    {
      o[b / _bypk] |= src[b] << (8 * (b % _bypk));
    }
    return out;
  }

public:
  TupleMap(size_t bitWidth_)
      : _bypk((bitWidth_ + 7) / 8), _data() {}

  // increment count for a key, initialising to 1
  void increment(const Int *ptr_, std::size_t len_)
  {
    auto key = _compress(ptr_, len_);
    std::string_view keyView(reinterpret_cast<const char *>(key.data()), key.size());
    _data.lazy_emplace_l(
        keyView,
        [](Map::value_type &kv)
        { ++kv.second; },
        [&](const Map::constructor &ctor)
        { ctor(std::move(key), /*initial*/(Val)1u); });
  }

  // Set meta-colour to value, constructing if it doesn't exist.
  void set(const Int *ptr_, std::size_t len_, Val z)
  {
    auto key = _compress(ptr_, len_);
    std::string_view keyView(reinterpret_cast<const char *>(key.data()), key.size());
    _data.lazy_emplace_l(
        keyView,
        [&](Map::value_type &kv)
        { kv.second = z; },
        [&](const Map::constructor &ctor)
        { ctor(std::move(key), /*initial*/z); });
  }

  // Only calls contructor if it DNE.
  template <class Ct>
  Val lazy_emplace(const Int *ptr_, std::size_t len_, Ct fn)
  {
    auto key = _compress(ptr_, len_);
    std::string_view keyView(reinterpret_cast<const char *>(key.data()), key.size());
    Val v;
    _data.lazy_emplace_l(
        keyView,
        [&](Map::value_type &kv)
        { v = kv.second; },
        [&](const Map::constructor &ctor)
        { v = fn(); ctor(std::move(key), /*initial*/v); });
    return v;
  }

  Map::const_iterator find(const Int *ptr_, std::size_t len_) const
  {
    auto key = _compress(ptr_, len_);
    std::string_view keyView(reinterpret_cast<const char *>(key.data()), key.size());
    return _data.find(key);
  }

  Map::const_iterator end() const
  {
    return _data.end();
  }

  void flush(std::ostream &bos) const
  {
    oneapi::tbb::parallel_for((size_t)0, _data.subcnt(), (size_t)1, [&](size_t submapIndex) -> void
                              {
                std::osyncstream os(bos);
                _data.with_submap(submapIndex, [&](const Map::EmbeddedSet& set)->void {
                    for (const auto &[key, count] : set) {
                        auto bigKey = _decompress(key);
                        auto it_ = bigKey.cbegin(), end_ = bigKey.cend();
                        os << *it_++;
                        while (it_ != end_)
                            os << ',' << *it_++;
                        os << '\t' << count << '\n';
                    }
                }); });
  }

  template <typename Fn>
  void forEach(Fn f_) const
  {
    for (const auto &[key, count] : _data)
    {
      auto bigKey = _decompress(key);
      f_(bigKey.cbegin(), bigKey.cend(), count);
    }
  }

  template <typename Fn>
  void pforEach(Fn f_) const
  {
    oneapi::tbb::parallel_for((size_t)0, _data.subcnt(), (size_t)1, [&](size_t submapIndex) -> void
                              { _data.with_submap(submapIndex, [&](const Map::EmbeddedSet &set) -> void
                                                  {
                    for (const auto &[key, count] : set) {
                        auto bigKey = _decompress(key);
                        auto it_ = bigKey.cbegin(), end_ = bigKey.cend();
                        f_(bigKey.cbegin(), bigKey.cend(), count);
                    } }); });
  }

  bool operator==(const TupleMap &rhs) const
  {
    return _data == rhs._data;
  }
};