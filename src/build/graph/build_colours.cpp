
#include "maki/build/graph/build_colours.hpp"
#include "maki/build/utils/bits.hpp"

colour_t Colours::getOrAssign(std::string &&seed) {
  colour_t id;
  ids.lazy_emplace_l(
      seed, [&](ColourMap::value_type &kv) { id = kv.second; },
      [&](const ColourMap::constructor &ctor) {
        id = _cid++;
        ctor(std::move(seed), /*initial*/ id);
      });
  return id;
}

BufferValue::BufferValue(const uint8_t *src_, size_t n_) : _data(0) {
  for (size_t i = 0; i < n_; ++i) {
    uint64_t val = *(src_ + i);
    _data |= val << (8 * i);
  }
}

void BufferValue::flush(uint8_t *dest_, size_t n_) const {
  for (size_t i = 0; i < n_; ++i) {
    uint8_t v = static_cast<uint8_t>((_data >> (8 * i)) & _msk);
    *(dest_ + i) = v;
  }
}

MetaColours::MetaColours(ColourMap &&m_, colour_t numFeatures)
    : ids(std::move(m_)), _mid(numFeatures), _nid(ids.size() + 1),
      r(required_bits(ids.size())), _occs(ids.size() + 1, 0) {}

bool MetaColours::assignable(const ColourVector &v) const {
  for (const colour_t &c : v)
    if (c < _mid)
      return true;
  return false;
}

uint64_t MetaColours::insert(const ColourVector &v) {
  for (const colour_t &c : v) {
    assert(c < _occs.size());
    ++std::atomic_ref{_occs[c]};
  }
  return id(v);
}

uint64_t MetaColours::id(const ColourVector &v) {
  if (v.size() == 1) {
    return v.front();
  } else {
    return r.lazy_emplace(v.data(), v.size(), [&] { return _nid++; });
  }
}

ColourRegistry toRegistry(MetaColours &&in_) {
  ColourRegistry reg;
  std::size_t numSeeds = in_.numColours();

  // move seed names
  reg.seeds.resize(numSeeds);
  for (auto &[name, i] : in_.ids) {
    reg.seeds[i] = std::move(name);
  }

  // move seed occurrences
  reg.occs = std::move(in_._occs);
  assert(reg.occs.size() == reg.seeds.size());

  // move meta colour mapping
  reg.metas.resize(in_.r.size());
  in_.r.pforEach([&](std::vector<colour_t> &&key, uint64_t value) {
    reg.metas[value - numSeeds] = std::move(key);
  });

  return reg;
}
