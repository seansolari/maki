
#pragma once

// -----------------------------------------------------------------------------
// Generalized SDSL int_vector sinks:
//  1) SdslIntVectorAppendOnDiskSink  : stream-append many sdsl::int_vector<>
//     into ONE SDSL-compatible on-disk int_vector using int_vector_buffer<0>.
//     The file has a single global width; this sink enforces/chooses it.
//  2) SdslIntVectorAppendInMemorySink: append many sdsl::int_vector<> payloads
//     into one in-memory sdsl::int_vector<>, widening automatically if needed.
// -----------------------------------------------------------------------------
// Notes:
//  - Payload type for both sinks: sdsl::int_vector<> (runtime width).
//  - On-disk sink: global width is fixed for the whole file; choose at ctor
//    or derived from the first payload. Later payloads must have width <=
//    global width (values must fit).
//  - In-memory sink: if a later payload needs more bits, we transparently
//    widen by repacking into a new vector and swapping.
//  - Both sinks expose write(const sdsl::int_vector<>&) and finalize().
//  - sink_payload<> specializations are provided.
// -----------------------------------------------------------------------------

#include <cstdint>
#include <string>

#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include "sink_manager.hpp"

// ---------------------- Utility: fast append helpers -------------------------
//
// We provide small helpers that do the right thing for arbitrary width.
// For simplicity and portability, we use element-wise reads/writes. If you
// need maximum throughput and you know widths are equal, you can implement
// a specialized block-copy path later.

namespace detail {

template <uint8_t mwidth>
inline void append_values_to_buffer(sdsl::int_vector_buffer<mwidth> &buf,
                                    const sdsl::int_vector<mwidth> &src) {
  const std::size_t n = src.size();
  for (std::size_t i = 0; i < n; ++i) {
    buf.push_back(src[i]); // width enforced by buffer configuration
  }
}

// assumes `dst` has been extended to be able to append `src`
template <uint8_t mwidth>
inline void append_values_to_int_vector(sdsl::int_vector<mwidth> &dst,
                                        std::size_t dst_offset,
                                        const sdsl::int_vector<mwidth> &src) {
  const std::size_t n = src.size();
  for (std::size_t i = 0; i < n; ++i) {
    dst[dst_offset + i] = src[i];
  }
}

} // namespace detail

// -----------------------------------------------------------------------------
// 1) On-disk sink: append many payloads into one int_vector file (width fixed)
// -----------------------------------------------------------------------------
template <uint8_t mwidth> class SdslIntVectorOnDiskSink {
  static_assert(mwidth > 0);

public:
  using Payload = sdsl::int_vector<mwidth>;

  // Enforce fixed width for the entire file.
  explicit SdslIntVectorOnDiskSink(const std::string &path,
                                   std::uint32_t mode = 0644,
                                   std::size_t buffer_bytes = (8ull << 20))
      : path_(path), mode_(mode),
        buf_(path, std::ios::out,
             buffer_bytes) // SDSL manages file stream/header
  {}

  // Append one payload; returns number of elements appended.
  std::uint64_t write(const Payload &iv) {
    const std::uint64_t n = iv.size();
    if (n == 0)
      return 0;
    detail::append_values_to_buffer(buf_, iv);
    return n;
  }

  void finalize() {
    buf_.close(); // flush & finalize header
  }

  std::uint64_t size() const noexcept { return buf_.size(); }
  const std::string &path() const noexcept { return path_; }

private:
  std::string path_;
  std::uint32_t mode_ = 0644;
  sdsl::int_vector_buffer<mwidth> buf_;
};

// -----------------------------------------------------------------------------
// In-memory sink: append payloads into a single sdsl::int_vector<>
// -----------------------------------------------------------------------------
template <uint8_t mwidth> class SdslIntVectorInMemorySink {
  static_assert(mwidth > 0);

public:
  using Payload = sdsl::int_vector<mwidth>;

  std::uint64_t write(const Payload &iv) {
    const std::uint64_t n = iv.size();
    if (n == 0)
      return 0;
    const std::uint64_t old_n = acc_.size();
    acc_.resize(old_n + n);
    detail::append_values_to_int_vector(acc_, old_n, iv);
    return size();
  }

  void finalize() {
    // Nothing to do; vector is already ready for use.
  }

  const sdsl::int_vector<mwidth> &data() const noexcept { return acc_; }
  sdsl::int_vector<mwidth> &data() noexcept { return acc_; }

  std::uint64_t size() const noexcept { return acc_.size(); }
  static constexpr std::uint8_t width() noexcept { return mwidth; }

private:
  sdsl::int_vector<mwidth> acc_; // accumulated data
};

// -----------------------------------------------------------------------------
// sink_payload<> specializations
// -----------------------------------------------------------------------------

template <uint8_t w> struct sink_payload<SdslIntVectorOnDiskSink<w>> {
  using type = sdsl::int_vector<w>;
};

template <uint8_t w> struct sink_payload<SdslIntVectorInMemorySink<w>> {
  using type = sdsl::int_vector<w>;
};
