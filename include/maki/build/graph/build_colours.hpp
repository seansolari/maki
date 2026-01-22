#pragma once
#include <cstdint>
#include <atomic>
#include <vector>
#include <string>
#include <mutex>

#include <oneapi/tbb/tbb_allocator.h>
#include <gtl/phmap.hpp>

#include "maki/core/utils/map.hpp"
#include "maki/build/utils/bits.hpp"

// ----------------------------------------------------------------
// Colours
// ----------------------------------------------------------------

// map seeds to colours
using ColourMap = gtl::parallel_flat_hash_map<
    std::string,
    uint64_t,
    gtl::priv::hash_default_hash<std::string>,
    gtl::priv::hash_default_eq<std::string>,
    oneapi::tbb::tbb_allocator<std::pair<const std::string, uint64_t>>,
    7UL,
    std::mutex>;

/**
 * Assign unique colours to seeds.
 */
struct Colours
{
  Colours() : ids(), _cid(0) {}
  Colours(const Colours&) =delete;
  Colours& operator=(const Colours&) =delete;

  ColourMap ids;

protected:
  std::atomic_uint64_t _cid;

public:
  uint64_t size() const noexcept
  {
    return _cid.load();
  }

  // Get colour for seed, assigning a new ID if it doesn't exist.
  uint64_t getOrAssign(std::string &&seed)
  {
    uint64_t id;
    ids.lazy_emplace_l(
        seed,
        [&](ColourMap::value_type &kv)
        { id = kv.second; },
        [&](const ColourMap::constructor &ctor)
        { id = ++_cid; ctor(std::move(seed), /*initial*/id); });
    return id;
  }
};

// Colour representation in a k-mer buffer
class BufferValue
{
  static constexpr uint64_t _msk = 0xFFu;

public:
  BufferValue() : _data(0ULL) {}
  BufferValue(uint64_t _edge) : _data(_edge) {}
  BufferValue(uint64_t _colour, uint64_t _edge) : _data((_colour << 3) | _edge) {}

  // read bytes from least- to most-significant
  explicit BufferValue(const uint8_t *src_, size_t n_) : _data(0)
  {
    for (size_t i = 0; i < n_; ++i)
    {
      uint64_t val = *(src_ + i);
      _data |= val << (8 * i);
    }
  }

public:
  bool operator==(const BufferValue &other) const noexcept { return _data == other._data; }
  bool operator!=(const BufferValue &other) const noexcept { return _data != other._data; }
  operator uint64_t() const { return _data; }
  uint64_t colour() const { return _data >> 3; }
  uint64_t edge() const { return _data & 0b111ULL; }
  uint64_t *data() { return &_data; }

public:
  // serialise bytes from least- to most-significant
  void flush(uint8_t *dest_, size_t n_) const
  {
    for (size_t i = 0; i < n_; ++i)
    {
      uint8_t v = static_cast<uint8_t>((_data >> (8 * i)) & _msk);
      *(dest_ + i) = v;
    }
  }

private:
  uint64_t _data;
};

// ----------------------------------------------------------------
// Colour Sets - Metacolours
// ----------------------------------------------------------------

// small buffer for arbitrary colour values
using ColourVector = std::vector<uint64_t, oneapi::tbb::tbb_allocator<uint64_t>>;

/**
 * Assign unique IDs to sets of colours (meta-colours), disjoint from the
 * original colour IDs.
 */
struct MetaColours
{
  MetaColours(ColourMap &&m_)
      : ids(std::move(m_)), _nid(ids.size()), r(ceil_log2(ids.size())) {}

protected:
  const ColourMap ids;
  std::atomic_uint64_t _nid;
  TupleMap<uint64_t, uint64_t> r;

public:
  // Get colour for seed, assigning a new ID if it doesn't exist.
  uint64_t getOrAssign(const ColourVector &v)
  {
    return r.lazy_emplace(v.data(), v.size(), [&]
                          { return ++_nid; });
  }
};
