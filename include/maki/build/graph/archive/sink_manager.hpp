
#pragma once

#include <cassert>
#include <concepts>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <mutex>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <oneapi/tbb/parallel_pipeline.h>

#if defined(__linux__) || defined(__APPLE__) || defined(__unix__)
#define MSINK_POSIX 1
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#if defined(__linux__)
#include <sys/mman.h>
#endif
#else
#define MSINK_POSIX 0
#endif

// ============================================================================
// Sink concepts and Sink-Bundle association
// ============================================================================

template <class Sink> struct sink_payload {};

// Helper alias to build a tuple of payloads for a sink pack
template <class... Sinks>
using BundlePayloadsFor = std::tuple<typename sink_payload<Sinks>::type...>;

template <class S, class P>
concept SinkAccepts = requires(S s, P &p) {
  { s.write(p) } -> std::same_as<std::size_t>;
  { s.finalize() } -> std::same_as<void>;
};

template <class... Sinks>
concept SinkSet =
    (SinkAccepts<Sinks, typename sink_payload<Sinks>::type> && ...);

// ============================================================================
// Sink and bundle managers
// ============================================================================

namespace detail {

template <class Tuple, class F, size_t... Is>
constexpr void tuple_for_each_impl(Tuple &&t, F &&f,
                                   std::index_sequence<Is...>) {
  (f(std::get<Is>(std::forward<Tuple>(t))), ...);
}

template <class Tuple, class F>
constexpr void tuple_for_each(Tuple &&t, F &&f) {
  constexpr size_t N = std::tuple_size_v<std::remove_reference_t<Tuple>>;
  tuple_for_each_impl(std::forward<Tuple>(t), std::forward<F>(f),
                      std::make_index_sequence<N>{});
}

template <class TupleA, class TupleB, class F, size_t... Is>
constexpr void tuple_for_each_pair_impl(TupleA &&a, TupleB &&b, F &&f,
                                        std::index_sequence<Is...>) {
  (f(std::get<Is>(std::forward<TupleA>(a)),
     std::get<Is>(std::forward<TupleB>(b))),
   ...);
}

template <class TupleA, class TupleB, class F>
constexpr void tuple_for_each_pair(TupleA &&a, TupleB &&b, F &&f) {
  static_assert(std::tuple_size_v<std::remove_reference_t<TupleA>> ==
                    std::tuple_size_v<std::remove_reference_t<TupleB>>,
                "tuple_for_each_pair: tuple sizes must match");
  constexpr size_t N = std::tuple_size_v<std::remove_reference_t<TupleA>>;
  tuple_for_each_pair_impl(std::forward<TupleA>(a), std::forward<TupleB>(b),
                           std::forward<F>(f), std::make_index_sequence<N>{});
}

} // namespace detail

/** Un-throttled container generation. If `obtain` is called when
 * no buffers are available, a new one is created.
 */
template <typename Container> struct Factory {
  Factory() =default;

  void release(std::unique_ptr<Container> &p) {
    std::lock_guard lock(mtx_);
    free_.push_back(std::move(p));
  }

  template <class... Args>
  std::unique_ptr<Container> obtain(Args&&... args) {
    {
      std::lock_guard lock(mtx_);
      if (!free_.empty()) {
        auto p = std::move(free_.back());
        free_.pop_back();
        return p;
      }
    }
    return std::make_unique<Container>(std::forward<Args>(args)...);
  }

protected:
  std::mutex mtx_;
  std::vector<std::unique_ptr<Container>> free_;
};

/**
 * Throttling stack that will block if `pop` is called when no data
 * is available.
 */
template <typename Container> class Stack {
public:
  // reserve (but don't allocate) space for this many items
  Stack(std::size_t reserve_) { free_.reserve(reserve_); }

  // allocate this number of items, in-place constructed
  template <class... Args>
  Stack(std::in_place_t, std::size_t size_, Args &...args) : free_(size_) {
    for (auto &p : free_) {
      p = std::make_unique<Container>(args...);
    }
  }

  std::unique_ptr<Container> pop() {
    std::unique_lock<std::mutex> lock(m_);
    cv_.wait(lock, [&] { return !free_.empty(); });
    auto p = std::move(free_.back());
    free_.pop_back();
    return p;
  }

  void push(std::unique_ptr<Container> p) {
    {
      std::lock_guard<std::mutex> lock(m_);
      free_.push_back(std::move(p));
    }
    cv_.notify_one();
  }

protected:
  std::mutex m_;
  std::condition_variable cv_;
  std::vector<std::unique_ptr<Container>> free_;
};

/**
 * Bundle type for a collection of payloads, associated with a chunk id.
 * Payload types are specified by the Sink they are associated with,
 * which is specified by overloading`sink_payload` and associating the
 * Chunk/Payload as the `type` member.
 */
template <class... Sinks> struct ChunkBundleT {
  static_assert(
      SinkSet<Sinks...>,
      "All sinks must satisfy SinkAccepts concept with their payload type");

  using Payloads = BundlePayloadsFor<Sinks...>;

  uint64_t id = 0;
  Payloads payloads;
};

/**
 * Payload bundles that threads request from. Prevents re-allocation for each
 * workload.
 */
template <class... Sinks> class BundlePool {
public:
  using Bundle = ChunkBundleT<Sinks...>;
  static constexpr size_t N = sizeof...(Sinks);

  // Reserve capacities per sink payload.bytes (in bytes)
  explicit BundlePool(size_t pool_size, size_t per_sink_reserve = 0)
      : stack_(pool_size) {
    for (size_t i = 0; i < pool_size; ++i) {
      auto b = std::make_unique<Bundle>();
      // Reserve bytes for each payload if it is PackedPayload (or has .bytes
      // vector)
      reserve_payloads_(b->payloads, per_sink_reserve);
      stack_.push(std::move(b));
    }
  }

  std::unique_ptr<Bundle> acquire() { return stack_.pop(); }

  void release(std::unique_ptr<Bundle> b) {
    clear_payloads_(b->payloads);
    b->id = 0;
    stack_.push(std::move(b));
  }

private:
  template <class PayloadTuple>
  static void reserve_payloads_(PayloadTuple &t, size_t reserve) {
    detail::tuple_for_each(t, [reserve](auto &payload) {
      if constexpr (requires { payload.reserve((std::size_t)0); }) {
        if (reserve)
          payload.reserve(reserve);
      }
    });
  }

  template <class PayloadTuple> static void clear_payloads_(PayloadTuple &t) {
    detail::tuple_for_each(t, [](auto &payload) {
      if constexpr (requires { payload.clear(); }) {
        payload.clear(); // keep capacity
      }
    });
  }

  Stack<Bundle> stack_;
};

/**
 * Manages distribution of chunks to sinks.
 */
template <class... Sinks> class MultiSink {
public:
  using Bundle = ChunkBundleT<Sinks...>;
  MultiSink(Sinks&& ...sinks) : sinks_(std::forward<Sinks>(sinks)...) {}

  MultiSink(const MultiSink&) =delete;
  MultiSink& operator=(const MultiSink&) =delete;

  void write_bundle(Bundle &b) {
    detail::tuple_for_each_pair(b.payloads, sinks_,
                                [](auto &payload, auto &sink) {
                                  // Skip empty payloads if they expose .empty()
                                  if constexpr (requires { payload.empty(); }) {
                                    if (payload.empty())
                                      return;
                                  }
                                  sink.write(payload);
                                });
  }

  void finalize() {
    detail::tuple_for_each(sinks_, [](auto &sink) { sink.finalize(); });
  }

  auto &sinks() { return sinks_; }
  const auto &sinks() const { return sinks_; }

private:
  std::tuple<Sinks...> sinks_;
};

/**
 * These class manages writing of payloads to sinks. It waits until contiguous
 * payloads are returned, holding non-contiguous payloads in a queue, and then
 * sinking as many contiguous payloads as possible.
 */

struct NextChunk {
  NextChunk(size_t chunks)
      : chunk(std::make_shared<uint64_t>(0u)), total(chunks) {}

  uint64_t operator()(oneapi::tbb::flow_control &fc) const {
    if (*chunk == total) {
      fc.stop();
      return 0;
    } else {
      return (*chunk)++;
    }
  }

protected:
  std::shared_ptr<uint64_t> chunk;
  const uint64_t total;
};

template <class... Sinks> struct FlushBundle {
  using Bundle = ChunkBundleT<Sinks...>;
  using Cache = std::unordered_map<uint64_t, std::unique_ptr<Bundle>>;

  FlushBundle(const std::shared_ptr<BundlePool<Sinks...>> &pool,
              MultiSink<Sinks...> *multi)
      : next_(std::make_shared<uint64_t>(0)), pool_(pool), multi_(multi),
        pending_(std::make_shared<Cache>()) {}

  void operator()(std::unique_ptr<Bundle> b) const {
    pending_->emplace(b->id, std::move(b));
    flush_();
  }

protected:
  void flush_() const {
    while (!pending_->empty()) {
      auto it = pending_->find(*next_);
      if (it == pending_->end())
        break;
      
      multi_->write_bundle(*it->second);
      pool_->release(std::move(it->second));
      pending_->erase(it);

      ++(*next_);
    }
  }

  std::shared_ptr<uint64_t> next_;
  std::shared_ptr<BundlePool<Sinks...>> pool_;
  MultiSink<Sinks...> *multi_;
  std::shared_ptr<Cache> pending_;
};

template <class T, class... Sinks>
concept BundleProducer = requires(const T fn, T gn, uint64_t i, std::shared_ptr<BundlePool<Sinks...>> &p) {
  { fn(i) } -> std::same_as<std::unique_ptr<ChunkBundleT<Sinks...>>>;
  { gn.setPool(p) } -> std::same_as<void>;
};

template <class T, class... Sinks>
  requires BundleProducer<T, Sinks...>
void ProcessChunks(T op, MultiSink<Sinks...> &sinks, std::size_t pool_size,
                   size_t per_sink_reserve, std::uint64_t expected_chunks) {
  using Bundle = ChunkBundleT<Sinks...>;
  auto pool =
      std::make_shared<BundlePool<Sinks...>>(pool_size, per_sink_reserve);
  op.setPool(pool);
  oneapi::tbb::parallel_pipeline(
      pool_size,
      oneapi::tbb::make_filter<void, uint64_t>(
          oneapi::tbb::filter_mode::serial_in_order, NextChunk{expected_chunks}) &
          oneapi::tbb::make_filter<uint64_t, std::unique_ptr<Bundle>>(
              oneapi::tbb::filter_mode::parallel, std::move(op)) &
          oneapi::tbb::make_filter<std::unique_ptr<Bundle>, void>(
              oneapi::tbb::filter_mode::serial_in_order, FlushBundle{pool, &sinks}));
}
