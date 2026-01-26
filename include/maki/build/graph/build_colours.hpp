#pragma once
#include <atomic>
#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

#include <gtl/phmap.hpp>
#include <oneapi/tbb/tbb_allocator.h>

#include "maki/build/utils/bits.hpp"
#include "maki/core/graph/colours.hpp"
#include "maki/core/utils/map.hpp"

// ----------------------------------------------------------------
// Colours
// ----------------------------------------------------------------

// map seeds to colours
using ColourMap = gtl::parallel_flat_hash_map<
    std::string, colour_t, gtl::priv::hash_default_hash<std::string>,
    gtl::priv::hash_default_eq<std::string>,
    oneapi::tbb::tbb_allocator<std::pair<const std::string, colour_t>>, 7UL,
    std::mutex>;

/**
 * Assign unique colours to seeds.
 */
struct Colours {
  Colours() : ids(), _cid(1 /*colours start from 1*/) {}
  Colours(const Colours &) = delete;
  Colours &operator=(const Colours &) = delete;

  ColourMap ids;

protected:
  std::atomic<colour_t> _cid;

public:
  colour_t size() const noexcept { return _cid.load(); }

  // Get colour for seed, assigning a new ID if it doesn't exist.
  colour_t getOrAssign(std::string &&seed);
};

// Colour representation in a k-mer buffer
class BufferValue {
  static constexpr uint64_t _msk = 0xFFu;

public:
  BufferValue() : _data(0ULL) {}
  BufferValue(uint64_t _edge) : _data(_edge) {}
  BufferValue(uint64_t _colour, uint64_t _edge)
      : _data((_colour << 3) | _edge) {}

  // read bytes from least- to most-significant
  explicit BufferValue(const uint8_t *src_, size_t n_);

public:
  bool operator==(const BufferValue &other) const noexcept {
    return _data == other._data;
  }
  bool operator!=(const BufferValue &other) const noexcept {
    return _data != other._data;
  }
  operator uint64_t() const { return _data; }
  colour_t colour() const { return static_cast<colour_t>(_data >> 3); }
  uint64_t edge() const { return _data & 0b111ULL; }
  uint64_t *data() { return &_data; }

public:
  // serialise bytes from least- to most-significant
  void flush(uint8_t *dest_, size_t n_) const;

private:
  uint64_t _data;
};

// ----------------------------------------------------------------
// Colour Sets - Metacolours
// ----------------------------------------------------------------

// small buffer for arbitrary colour values
using ColourVector =
    std::vector<colour_t, oneapi::tbb::tbb_allocator<colour_t>>;

/**
 * Assign unique IDs to sets of colours (meta-colours), disjoint from the
 * original colour IDs.
 */
struct MetaColours {
  MetaColours(ColourMap &&m_, uint64_t numFeatures);

protected:
  ColourMap ids;                  // stores seed IDs
  uint64_t _mid;                  // marks colour IDs from features vs filters
  std::atomic_uint64_t _nid;      // number of nodes currently assigned
  TupleMap<colour_t, uint64_t> r; // map tuples of colours to meta colour IDs
  std::vector<uint64_t> _occs;    // occurrences of every colour

public:
  inline colour_t numColours() const { return ids.size(); }
  inline std::size_t maxColourWidth() const { return ceil_log2(numColours()); }

  // Do any colours come from features? Or are they all filters?
  bool assignable(const ColourVector &v) const;

  // Get colour for seed, assigning a new ID if it doesn't exist.
  uint64_t id(const ColourVector &v);

  // Insert a colour and increment occurrence counts
  uint64_t insert(const ColourVector &v);

  friend ColourRegistry toRegistry(MetaColours &&in_);
};

// ----------------------------------------------------------------
// Finalised Registry
// ----------------------------------------------------------------

ColourRegistry toRegistry(MetaColours &&in_);
