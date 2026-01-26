
#include "maki/build/graph/build_colours.hpp"

uint64_t Colours::getOrAssign(std::string &&seed) {
  uint64_t id;
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

MetaColours::MetaColours(ColourMap &&m_, uint64_t numFeatures)
    : ids(std::move(m_)), _mid(numFeatures), _nid(ids.size()),
      r(ceil_log2(ids.size())), _occs(ids.size(), 0) {}

bool MetaColours::assignable(const ColourVector &v) const {
  for (const auto &c : v)
    if (c < _mid)
      return true;
  return false;
}

uint64_t MetaColours::insert(const ColourVector &v) {
  // increment occurrence counts
  for (const auto &c : v) {
    assert(c < _occs.size());
    ++std::atomic_ref{_occs[c]};
  }

  // assign ID
  if (v.size() == 1) {
    return v.front();
  } else {
    return id(v);
  }
}

uint64_t MetaColours::id(const ColourVector &v) {
  return r.lazy_emplace(v.data(), v.size(), [&] { return ++_nid; });
}

ColourRegistry toRegistry(MetaColours &&in_) {
  ColourRegistry reg;
  uint32_t numSeeds = in_.numColours();

  // move seed names
  reg.seeds.resize(numSeeds);
  for (auto &[name, i] : in_.ids) {
    reg.seeds[i] = std::move(name);
  }

  // move seed occurrences
  reg.occs = std::move(in_._occs);

  // move meta colour mapping
  reg.metas.resize(in_.r.size());
  in_.r.pforEach([&](std::vector<uint32_t> &&key, uint32_t value) {
    reg.metas[value - numSeeds] = std::move(key);
  });

  return reg;
}
