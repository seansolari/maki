
#pragma once

#include <array>
#include <cassert>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <exception>
#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

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

template <class Sink>
struct sink_payload
{
};

// Helper alias to build a tuple of payloads for a sink pack
template <class... Sinks>
using BundlePayloadsFor = std::tuple<typename sink_payload<Sinks>::type...>;

template <class S, class P>
concept SinkAccepts = requires(S s, const P &p) {
  { s.write(p) } -> std::same_as<void>;
  { s.finalize() } -> std::same_as<void>;
};

template <class... Sinks>
concept SinkSet = (SinkAccepts<Sinks, typename sink_payload<Sinks>::type> && ...);

// ============================================================================
// Sink and bundle managers
// ============================================================================

namespace detail
{

  template <class Tuple, class F, size_t... Is>
  constexpr void tuple_for_each_impl(Tuple &&t, F &&f, std::index_sequence<Is...>)
  {
    (f(std::get<Is>(std::forward<Tuple>(t))), ...);
  }

  template <class Tuple, class F>
  constexpr void tuple_for_each(Tuple &&t, F &&f)
  {
    constexpr size_t N = std::tuple_size_v<std::remove_reference_t<Tuple>>;
    tuple_for_each_impl(std::forward<Tuple>(t), std::forward<F>(f), std::make_index_sequence<N>{});
  }

  template <class TupleA, class TupleB, class F, size_t... Is>
  constexpr void tuple_for_each_pair_impl(TupleA &&a, TupleB &&b, F &&f, std::index_sequence<Is...>)
  {
    (f(std::get<Is>(std::forward<TupleA>(a)), std::get<Is>(std::forward<TupleB>(b))), ...);
  }

  template <class TupleA, class TupleB, class F>
  constexpr void tuple_for_each_pair(TupleA &&a, TupleB &&b, F &&f)
  {
    static_assert(std::tuple_size_v<std::remove_reference_t<TupleA>> ==
                      std::tuple_size_v<std::remove_reference_t<TupleB>>,
                  "tuple_for_each_pair: tuple sizes must match");
    constexpr size_t N = std::tuple_size_v<std::remove_reference_t<TupleA>>;
    tuple_for_each_pair_impl(std::forward<TupleA>(a), std::forward<TupleB>(b), std::forward<F>(f),
                             std::make_index_sequence<N>{});
  }

} // namespace detail

/**
 * Bundle type for a collection of chunks. Chunk types are specified by the
 * Sink they are associated with, which is specified by overloading `sink_payload`
 * and associating the Chunk/Payload as the `type` member.
 */
template <class... Sinks>
struct ChunkBundleT
{
  static_assert(SinkSet<Sinks...>, "All sinks must satisfy SinkAccepts concept with their payload type");

  using Payloads = BundlePayloadsFor<Sinks...>;

  uint64_t id = 0;
  Payloads payloads;
};

/**
 * Payload bundles that threads request from. Prevents re-allocation for each workload.
 */
template <class... Sinks>
class BundlePool
{
public:
  using Bundle = ChunkBundleT<Sinks...>;
  static constexpr size_t N = sizeof...(Sinks);

  // Reserve capacities per sink payload.bytes (in bytes)
  explicit BundlePool(size_t pool_size, const std::array<size_t, N> &per_sink_reserve = {})
  {
    free_.reserve(pool_size);
    for (size_t i = 0; i < pool_size; ++i)
    {
      auto b = std::make_unique<Bundle>();
      // Reserve bytes for each payload if it is PackedPayload (or has .bytes vector)
      reserve_payloads_(b->payloads, per_sink_reserve);
      free_.push_back(std::move(b));
    }
  }

  std::unique_ptr<Bundle> acquire()
  {
    std::unique_lock<std::mutex> lock(m_);
    cv_.wait(lock, [&]
             { return !free_.empty(); });
    auto p = std::move(free_.back());
    free_.pop_back();
    return p;
  }

  void release(std::unique_ptr<Bundle> b)
  {
    clear_payloads_(b->payloads);
    b->id = 0;
    {
      std::lock_guard<std::mutex> lock(m_);
      free_.push_back(std::move(b));
    }
    cv_.notify_one();
  }

private:
  template <class PayloadTuple>
  static void reserve_payloads_(PayloadTuple &t, const std::array<size_t, N> &reserve)
  {
    detail::tuple_for_each_pair(t, reserve, {
      // Default payload is PackedPayload -> has .bytes
      if constexpr (requires { payload.bytes.reserve(0); })
      {
        if (cap)
          payload.bytes.reserve(cap);
      }
    });
  }

  template <class PayloadTuple>
  static void clear_payloads_(PayloadTuple &t)
  {
    detail::tuple_for_each(t, {
      if constexpr (requires { payload.bytes.clear(); })
      {
        payload.bytes.clear(); // keep capacity
      }
      if constexpr (requires { payload.elem_count = 0; })
      {
        payload.elem_count = 0;
      }
      if constexpr (requires { payload.bit_width = 0; })
      {
        payload.bit_width = 0;
      }
    });
  }

  std::mutex m_;
  std::condition_variable cv_;
  std::vector<std::unique_ptr<Bundle>> free_;
};

/**
 * Holds payload bundles that can't be written to a sink yet.
 */
template <class Bundle>
class BundleQueue
{
public:
  void submit(std::unique_ptr<Bundle> b)
  {
    {
      std::lock_guard<std::mutex> lock(m_);
      q_.push_back(std::move(b));
    }
    cv_.notify_one();
  }

  bool take(std::unique_ptr<Bundle> &out)
  {
    std::unique_lock<std::mutex> lock(m_);
    cv_.wait(lock, [&]
             { return !q_.empty() || closed_; });
    if (q_.empty())
      return false;
    out = std::move(q_.front());
    q_.pop_front();
    return true;
  }

  void close()
  {
    std::lock_guard<std::mutex> lock(m_);
    closed_ = true;
    cv_.notify_all();
  }

private:
  std::mutex m_;
  std::condition_variable cv_;
  std::deque<std::unique_ptr<Bundle>> q_;
  bool closed_ = false;
};

/**
 * Manages distribution of chunks to sinks.
 */
template <class... Sinks>
class MultiSink
{
public:
  using Bundle = ChunkBundleT<Sinks...>;
  MultiSink(std::tuple<Sinks...> sinks) : sinks_(std::move(sinks)) {}

  void write_bundle(const Bundle &b)
  {
    detail::tuple_for_each_pair(b.payloads, sinks_, {
      // Skip empty payloads if they expose .empty()
      if constexpr (requires { payload.empty(); })
      {
        if (payload.empty())
          return;
      }
      sink.write(payload);
    });
  }

  void finalize()
  {
    detail::tuple_for_each(sinks_, { sink.finalize(); });
  }

  auto &sinks() { return sinks_; }
  const auto &sinks() const { return sinks_; }

private:
  std::tuple<Sinks...> sinks_;
};

/**
 * This class manages writing of payloads to sinks. It waits until contiguous payloads are returned,
 * holding non-contiguous payloads in a queue, and then sinking as many contiguous payloads
 * as possible.
 */
template <class... Sinks>
class OrderedBundleWriter
{
public:
  using Bundle = ChunkBundleT<Sinks...>;

  OrderedBundleWriter(BundleQueue<Bundle> &q, BundlePool<Sinks...> &pool, MultiSink<Sinks...> &multi, uint64_t expected_bundles)
      : q_(q), pool_(pool), multi_(multi), expected_(expected_bundles)
  {
    assert(expected_ > 0);
  }

  void run()
  {
    std::unique_ptr<Bundle> b;
    while (q_.take(b))
    {
      pending_.emplace(b->id, std::move(b));
      flush_();
      if (done_)
        break;
    }
    flush_();
    multi_.finalize();
  }

private:
  void flush_()
  {
    for (;;)
    {
      auto it = pending_.find(next_);
      if (it == pending_.end())
        break;

      multi_.write_bundle(*it->second);
      pool_.release(std::move(it->second));
      pending_.erase(it);

      ++next_;
      if (expected_ && next_ == expected_)
      {
        done_ = true;
        break;
      }
    }
  }

  BundleQueue<Bundle> &q_;
  BundlePool<Sinks...> &pool_;
  MultiSink<Sinks...> &multi_;
  const uint64_t expected_;
  uint64_t next_ = 0;
  bool done_ = false;
  std::unordered_map<uint64_t, std::unique_ptr<Bundle>> pending_;
};

template <class... Sinks>
class PipelineFacade
{
public:
  using Bundle = ChunkBundleT<Sinks...>;
  using Pool = BundlePool<Sinks...>;
  using Queue = BundleQueue<Bundle>;
  using Multi = MultiSink<Sinks...>;
  using Writer = OrderedBundleWriter<Sinks...>;

  PipelineFacade(std::tuple<Sinks...> sinks,
                 std::size_t pool_size,
                 const std::array<std::size_t, sizeof...(Sinks)> &per_sink_reserve,
                 std::uint64_t expected_bundles = 0)
      : pool_(pool_size, per_sink_reserve),
        multi_(std::move(sinks)),
        writer_(queue_, pool_, multi_, expected_bundles)
  {
  }

  Pool &pool() { return pool_; }
  Queue &queue() { return queue_; }
  Multi &multi() { return multi_; }
  Writer &writer() { return writer_; }

private:
  Pool pool_;
  Queue queue_;
  Multi multi_;
  Writer writer_;
};
