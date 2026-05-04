
#pragma once
#include "maki/build/graph/archive/sink_manager.hpp"
#include "maki/core/graph/archive/counts.hpp"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/util.hpp>

struct CountBufferSink {
  using Payload = CountBuffer;
  std::size_t write(const CountBuffer &);
  void finalize() {} // do nothing
  CountBuffer data;
};

template <> struct sink_payload<CountBufferSink> {
  using type = CountBuffer;
};
