#pragma once
#include <atomic>
#include <cstdint>

namespace colour
{
  namespace encoding
  {
    struct FeatureIdAssigner
    {
      FeatureIdAssigner() : _cid(0) {}
    protected:
      std::atomic_uint64_t _cid;
    public:
      uint64_t nextId() { return _cid++; }
      uint64_t current() const { return _cid.load(); }
    };

  } // namespace encoding

} // namespace colour

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
