#pragma once
#include <cstdint>
#include <vector>
#include <string_view>
#include <mutex>
#include <oneapi/tbb/tbb_allocator.h>
#include <gtl/phmap.hpp>

// small buffer for arbitrary colour values
using ColourVector = std::vector<uint64_t, oneapi::tbb::tbb_allocator<uint64_t>>;

/**
 * Registry that maps a set of colours to a unique ID. Colour sets
 * are byte-packed for compression.
 */
class ColourSetRegistry
{
  using Key = std::vector<uint8_t, oneapi::tbb::tbb_allocator<uint8_t>>;

  struct ByteVectorHash
  {
    using is_transparent = void;
    size_t operator()(const Key &v) const noexcept;
    size_t operator()(std::string_view sv) const noexcept;
  };

  struct ByteVectorEq
  {
    using is_transparent = void;
    bool operator()(const Key &a, const Key &b) const noexcept;
    bool operator()(const Key &a, std::string_view b) const noexcept;
  };

  using Map = gtl::parallel_flat_hash_map<
      Key,            // Key
      uint64_t,       // Mapped value
      ByteVectorHash, // Hash functor
      ByteVectorEq,   // Equality comparator
      oneapi::tbb::tbb_allocator<std::pair<const Key, uint64_t>>,
      /*N=*/7,   // number of submaps is 2^N => 128
      std::mutex // per-submap mutex for thread safety
      >;

  size_t _bypk; // bytes per key
  Map _data;

protected:
  Key compressU64(const ColourVector &src) const;
  ColourVector decompressU64(const Key &bytes) const;

public:
  ColourSetRegistry(size_t bitWidth_);
  void increment(const ColourVector &v);
  void set(const ColourVector &v, uint32_t z);
  Map::const_iterator find(const ColourVector &v) const;
  Map::const_iterator end() const;
  void flush(std::ostream &) const;
  template <typename Fn>
  void forEach(Fn f_) const;
  template <typename Fn>
  void pforEach(size_t grainsize, Fn f_) const;
  bool operator==(const ColourSetRegistry &rhs) const;
};

struct FeatureIdAssigner
{
  FeatureIdAssigner() : _cid(0) {}

protected:
  std::atomic_uint64_t _cid;

public:
  uint64_t nextId() { return _cid++; }
  uint64_t current() const { return _cid.load(); }
};

class BufferValue
{
  static constexpr uint64_t _msk = 0xFFu;

public:
  BufferValue() : _data(0ULL) {}
  BufferValue(uint64_t _edge) : _data(_edge) {}
  BufferValue(uint64_t _colour, uint64_t _edge) : _data((_colour << 3) | _edge) {}
  // read bytes from least- to most-significant
  explicit BufferValue(const uint8_t *src_, size_t n_);

public:
  bool operator==(const BufferValue &other) const noexcept { return _data == other._data; }
  bool operator!=(const BufferValue &other) const noexcept { return _data != other._data; }
  operator uint64_t() const { return _data; }
  uint64_t colour() const { return _data >> 3; }
  uint64_t edge() const { return _data & 0b111ULL; }
  uint64_t *data() { return &_data; }

public:
  // serialise bytes from least- to most-significant
  void flush(uint8_t *dest_, size_t n_) const;

private:
  uint64_t _data;
};
