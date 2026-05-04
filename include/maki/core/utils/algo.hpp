#pragma once
#include <algorithm>
#include <sdsl/bit_vectors.hpp>

template <class O> void vector_append(O &lhs, const O &rhs) {
  lhs.reserve(lhs.size() + rhs.size());
  lhs.insert(lhs.end(), rhs.cbegin(), rhs.cend());
}

template <class O, typename V>
void vector_append_add(O &lhs, const O &rhs, V dv) {
  std::size_t i = lhs.size();
  vector_append(lhs, rhs);
  std::for_each(lhs.begin() + i, lhs.end(), [dv](auto &x) { x += dv; });
}
