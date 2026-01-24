
#pragma once
#include <cstdint>
#include <memory>
#include <optional>
#include <ranges>
#include <type_traits>
#include <utility>

#include "seq_io.hpp"

/**
 * Models any input range over type `T`.
 */
template <class T> class poly_input_range {
  /**
   * Interface required for a type-erased input range.
   */
  struct concept_t {
    virtual ~concept_t() = default;
    virtual std::optional<T> next() = 0; // pull-based cursor
  };

  /**
   * Implements `concept_t` for any `std::ranges::input_range` whose
   * range type can be converted to `T`.
   */
  template <std::ranges::input_range R>
    requires std::convertible_to<std::ranges::range_value_t<R>, T>
  struct model final : concept_t {
    std::ranges::iterator_t<R> it_;
    std::ranges::sentinel_t<R> end_;

    explicit model(R r)
        : it_(std::ranges::begin(r)), end_(std::ranges::end(r)) {}

    std::optional<T> next() override {
      if (it_ == end_)
        return std::nullopt;
      if constexpr (std::is_reference_v<std::ranges::range_reference_t<R>>) {
        return *it_++;
      } else {
        T v = *it_;
        ++it_;
        return v;
      }
    }
  };

  struct model_single_owned final : concept_t {
    T value_;
    bool emitted_ = false;
    explicit model_single_owned(T v) : value_(std::move(v)) {}

    std::optional<T> next() override {
      if (emitted_)
        return std::nullopt;
      emitted_ = true;
      if constexpr (std::is_move_constructible_v<T>) {
        return std::move(value_);
      } else {
        return value_;
      }
    }
  };

  struct model_single_ref final : concept_t {
    const T *ptr_ = nullptr;
    bool emitted_ = false;
    explicit model_single_ref(const T *p) : ptr_(p) {}

    std::optional<T> next() override {
      if (emitted_ || !ptr_)
        return std::nullopt;
      emitted_ = true;
      return *ptr_; // note: returns by value (input range)
    }
  };

  std::unique_ptr<concept_t> self_;

  /**
   * Iterator over input range that pulls items via next().
   */
  class iter {
    concept_t *p_ = nullptr;
    std::optional<T> cur_;

  public:
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::input_iterator_tag;

    iter() = default;
    explicit iter(concept_t *p, bool at_end) : p_(p) {
      if (p_ && !at_end)
        cur_ = p_->next();
    }
    T operator*() const { return *cur_; }
    iter &operator++() {
      cur_ = p_->next();
      return *this;
    }
    void operator++(int) { ++(*this); }
    friend bool operator==(const iter &a, const iter &b) {
      // end == end or both cur_ empty means end
      if (a.p_ != b.p_)
        return false; // different ranges
      return (!a.cur_.has_value()) == (!b.cur_.has_value());
    }
  };

public:
  poly_input_range() = default;

  template <std::ranges::input_range R>
    requires std::convertible_to<std::ranges::range_value_t<R>, T>
  poly_input_range(R &&r) : self_(std::make_unique<model<R>>(std::move(r))) {}

  explicit poly_input_range(T &&v)
      : self_(std::make_unique<model_single_owned>(std::move(v))) {}

  explicit poly_input_range(const T &v)
      : self_(std::make_unique<model_single_ref>(&v)) {}

  // range interface
  iter begin() {
    if (!self_)
      return iter{};
    else
      return iter(self_.get(), false);
  }
  iter end() { return iter(self_.get(), true); }
};

class SequenceFragment {
public:
  SequenceFragment(Dna4SequenceConstIter begin, Dna4SequenceConstIter end,
                   uint64_t id, bool terminal)
      : _begin(begin), _end(end), _fid(id), _terminal(terminal) {}

  Dna4SequenceConstIter begin() const noexcept { return _begin; }
  Dna4SequenceConstIter end() const noexcept { return _end; }
  uint64_t id() const noexcept { return _fid; }
  bool endIsTerminal() const noexcept { return _terminal; }
  std::size_t size() const noexcept { return _end - _begin; }

protected:
  Dna4SequenceConstIter _begin;
  Dna4SequenceConstIter _end;
  uint64_t _fid;
  bool _terminal;
};

struct SequenceContainer {
  virtual ~SequenceContainer() = default;

  // ------------------------------------------------------------
  // Terminals
  // ------------------------------------------------------------

  /**
   * Count the number of starting terminals in the container.
   */
  virtual std::size_t numTerminals(std::size_t k_) const = 0;

  /**
   * Return terminal sequence fragments.
   */
  virtual poly_input_range<SequenceFragment> terminals() const = 0;

  // ------------------------------------------------------------
  // Terminals
  // ------------------------------------------------------------

  /**
   * Count the number of k-mers in the container, including end
   * terminal edges.
   */
  virtual std::size_t numKmers(std::size_t k_) const = 0;

  /**
   * View over fragments in the container.
   */
  virtual poly_input_range<SequenceFragment>
  fragments(std::size_t k_) const = 0;
};
