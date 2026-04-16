
#pragma once
#include <cstdint>
#include <iterator>
#include <memory>
#include <ranges>
#include <utility>

#include "seq_io.hpp"

template <typename T> class poly_input_range {
public:
  using value_type = T;

  poly_input_range() = delete;

  // Construct from any range/view
  template <std::ranges::input_range R>
    requires std::same_as<std::ranges::range_value_t<R>, T>
  poly_input_range(R &&r)
      : self_(std::make_unique<range_model<R>>(std::forward<R>(r))) {}

  // Construct from a single value (owned)
  poly_input_range(T value)
      : self_(std::make_unique<single_value_model>(std::move(value))) {}

  // Construct from a single reference (non-owning)
  poly_input_range(std::reference_wrapper<T> ref)
      : self_(std::make_unique<single_ref_model>(ref.get())) {}

  // empty factory
  static poly_input_range empty_range() {
    return poly_input_range(std::ranges::empty_view<T>{});
  }

  struct iterator {
    using iterator_category = std::input_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;

    struct iter_concept {
      virtual ~iter_concept() = default;
      virtual T deref() const = 0;
      virtual void inc() = 0;
      virtual bool is_end() const = 0;
    };

  private:
    std::unique_ptr<iter_concept> iter_;

  public:
    iterator() = default;
    explicit iterator(std::unique_ptr<iter_concept> p) : iter_(std::move(p)) {}

    T operator*() const { return iter_->deref(); }
    iterator &operator++() {
      iter_->inc();
      return *this;
    }

    bool operator==(std::default_sentinel_t) const { return iter_->is_end(); }
  };

public:
  iterator begin() { return iterator(self_->begin()); }

  std::default_sentinel_t end() { return {}; }

private:
  struct range_concept {
    virtual ~range_concept() = default;
    virtual std::unique_ptr<typename iterator::iter_concept> begin() = 0;
  };

  template <typename R> struct range_model : range_concept {
    R range;

    explicit range_model(R &&r) : range(std::move(r)) {}

    struct iter_model : iterator::iter_concept {
      std::ranges::iterator_t<R> it;
      std::ranges::sentinel_t<R> end;

      iter_model(R &r) : it(std::ranges::begin(r)), end(std::ranges::end(r)) {}

      T deref() const override { return *it; }

      void inc() override { ++it; }

      bool is_end() const override { return it == end; }
    };

    std::unique_ptr<typename iterator::iter_concept> begin() override {
      return std::make_unique<iter_model>(range);
    }
  };

  struct single_value_model : range_concept {
    T value;

    explicit single_value_model(T v) : value(std::move(v)) {}

    struct iter_model : iterator::iter_concept {
      const T *ptr;
      bool done = false;

      explicit iter_model(const T &v) : ptr(&v) {}

      T deref() const override { return *ptr; }

      void inc() override { done = true; }

      bool is_end() const override { return done; }
    };

    std::unique_ptr<typename iterator::iter_concept> begin() override {
      return std::make_unique<iter_model>(value);
    }
  };

  struct single_ref_model : range_concept {
    T *ptr;

    explicit single_ref_model(T *p) : ptr(p) {}

    struct iter_model : iterator::iter_concept {
      T *ptr;
      bool done = false;

      explicit iter_model(T *p) : ptr(p) {}

      T deref() const override { return *ptr; }

      void inc() override { done = true; }

      bool is_end() const override { return done; }
    };

    std::unique_ptr<typename iterator::iter_concept> begin() override {
      return std::make_unique<iter_model>(ptr);
    }
  };

  std::unique_ptr<range_concept> self_;
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
