
#pragma once
#include <concepts>
#include <cstdint>
#include <iterator>
#include <ranges>
#include <type_traits>
#include <utility>

#include "seqan3/alphabet/nucleotide/dna4.hpp"
#include "seqan3/alphabet/nucleotide/dna5.hpp"

template <typename I>
concept random_dna4_iter =
    std::random_access_iterator<I> &&
    std::convertible_to<std::iter_value_t<I>, seqan3::dna4>;

template <typename I>
concept random_dna5_iter =
    std::random_access_iterator<I> &&
    std::convertible_to<std::iter_value_t<I>, seqan3::dna5>;

template <typename R>
concept random_dna4_range =
    std::ranges::random_access_range<R> &&
    std::convertible_to<std::ranges::range_reference_t<R>, seqan3::dna4>;

template <typename R>
concept random_dna5_range =
    std::ranges::random_access_range<R> &&
    std::convertible_to<std::ranges::range_reference_t<R>, seqan3::dna5>;

template <random_dna4_range range_t> class SequenceFragment {
public:
  using range_type = range_t;

  explicit SequenceFragment(range_t &&rng, uint64_t id, bool terminal)
      : _rng(std::move(rng)), _fid(id), _terminal(terminal) {}

  inline constexpr range_t &data() noexcept { return _rng; }
  inline constexpr range_t const &data() const noexcept { return _rng; }
  inline constexpr auto begin() noexcept { return std::ranges::begin(_rng); }
  inline constexpr auto begin() const noexcept { return std::ranges::begin(_rng); }
  inline constexpr auto cbegin() const noexcept { return std::ranges::cbegin(_rng); }
  inline constexpr auto end() noexcept { return std::ranges::end(_rng); }
  inline constexpr auto end() const noexcept { return std::ranges::end(_rng); }
  inline constexpr auto cend() const noexcept { return std::ranges::end(_rng); }
  inline constexpr uint64_t id() const noexcept { return _fid; }
  inline constexpr bool endIsTerminal() const noexcept { return _terminal; }
  inline constexpr std::size_t size() const noexcept { return std::ranges::size(_rng); }

protected:
  range_t _rng;
  uint64_t _fid;
  bool _terminal;
};

template <typename T> struct is_sequence_fragment_t : std::false_type {};

template <random_dna4_range R>
struct is_sequence_fragment_t<SequenceFragment<R>> : std::true_type {};

template <typename T>
inline constexpr bool is_sequence_fragment_v = is_sequence_fragment_t<T>::value;

template <typename T>
concept is_sequence_fragment = is_sequence_fragment_v<T>;

template <typename R>
concept sequence_fragment_input_range =
    std::ranges::input_range<R> &&
    is_sequence_fragment<std::ranges::range_value_t<R>>;

template <typename T>
concept countable_sequence_holder = requires(T obj, std::size_t k_) {
  { obj.numTerminals(k_) } -> std::same_as<std::size_t>;
  { obj.numKmers(k_) } -> std::same_as<std::size_t>;
};

template <typename T>
concept sequence_like =
    countable_sequence_holder<T> && requires(T obj, std::size_t k_) {
      { obj.terminals() } -> is_sequence_fragment;
      { obj.fragments(k_) } -> is_sequence_fragment;
    };

template <typename T>
concept sequence_container_like =
    countable_sequence_holder<T> && requires(T obj, std::size_t k_) {
      { obj.terminals() } -> sequence_fragment_input_range;
      { obj.fragments(k_) } -> sequence_fragment_input_range;
    };

template <typename T>
concept stranded_sequence_container_like =
    countable_sequence_holder<T> && requires(T obj, std::size_t k_) {
      { obj.forwardTerminals() } -> sequence_fragment_input_range;
      { obj.reverseTerminals() } -> sequence_fragment_input_range;
      { obj.forwardFragments(k_) } -> sequence_fragment_input_range;
      { obj.reverseFragments(k_) } -> sequence_fragment_input_range;
    };

template <typename T>
concept sequence_fragment_container =
    sequence_like<T> || sequence_container_like<T> ||
    stranded_sequence_container_like<T>;

template <typename T>
concept container_span = requires(const T &obj, std::size_t i) {
  { obj.size() } -> std::same_as<std::size_t>;
  { obj[i] } -> sequence_fragment_container;
};
