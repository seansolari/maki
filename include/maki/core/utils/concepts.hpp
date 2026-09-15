#pragma once

#include <concepts>


template <typename T, typename V>
concept span_like = requires(const T &obj, std::size_t i) {
  { obj.size() } -> std::same_as<std::size_t>;
  { obj[i] } -> std::convertible_to<V&>;
};
