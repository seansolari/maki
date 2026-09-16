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

template <typename F, typename Tuple, std::size_t... Is>
decltype(auto) call_tail_impl(F &&f, Tuple &&t, std::index_sequence<Is...>) {
  return std::forward<F>(f)(std::get<Is + 2>(std::forward<Tuple>(t))...);
}

template <typename F, typename Tuple>
decltype(auto) call_tail(F &&f, Tuple &&t) {
  constexpr std::size_t N = std::tuple_size_v<std::remove_reference_t<Tuple>>;

  static_assert(N >= 2);

  return call_tail_impl(std::forward<F>(f), std::forward<Tuple>(t),
                        std::make_index_sequence<N - 2>{});
}
