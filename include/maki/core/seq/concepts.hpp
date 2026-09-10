
#pragma once
#include <concepts>
#include <cstdint>
#include <iterator>
#include <ranges>
#include <utility>

#include "seqan3/alphabet/nucleotide/dna4.hpp"

template <typename I>
concept random_dna4_iter =
    std::random_access_iterator<I> &&
    std::convertible_to<std::iter_value_t<I>, seqan3::dna4>;

template <typename R>
concept random_dna4_range =
    std::ranges::random_access_range<R> &&
    std::convertible_to<std::ranges::range_reference_t<R>, seqan3::dna4>;

template <typename T>
concept is_sequence_fragment_like = requires(const T &obj) {
  typename T::range_type;
  requires random_dna4_range<typename T::range_type>;

  { obj.data() } -> std::same_as<const typename T::range_type&>;
  { obj.id() } -> std::same_as<std::size_t>;
  { obj.endIsTerminal() } -> std::same_as<bool>;
  { obj.size() } -> std::same_as<std::size_t>;
};

template <random_dna4_range range_t> class SequenceFragment {
public:
  using range_type = range_t;

  explicit SequenceFragment(range_t &&rng, uint64_t id, bool terminal)
      : _rng(std::move(rng)), _fid(id), _terminal(terminal) {}

  range_t const &data() const noexcept { return _rng; }
  uint64_t id() const noexcept { return _fid; }
  bool endIsTerminal() const noexcept { return _terminal; }
  std::size_t size() const noexcept { return std::ranges::size(_rng); }

protected:
  range_t _rng;
  uint64_t _fid;
  bool _terminal;
};

template <typename R>
concept sequence_fragment_input_range =
    std::ranges::input_range<R> &&
    is_sequence_fragment_like<std::ranges::range_value_t<R>>;

template <typename T>
concept countable_sequence_holder = requires(T obj, std::size_t k_) {
  { obj.numTerminals(k_) } -> std::same_as<std::size_t>;
  { obj.numKmers(k_) } -> std::same_as<std::size_t>;
};

template <typename T>
concept sequence_container_like =
    countable_sequence_holder<T> && requires(T obj, std::size_t k_) {
      { obj.terminals() } -> sequence_fragment_input_range;
      { obj.fragments(k_) } -> sequence_fragment_input_range;
    };

template <typename T>
concept sequence_like = countable_sequence_holder<T> && requires (T obj, std::size_t k_) {
  { obj.terminals() } -> is_sequence_fragment_like;
  { obj.fragments(k_) } -> is_sequence_fragment_like;
};
