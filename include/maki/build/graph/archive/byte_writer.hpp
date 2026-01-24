
#pragma once

// -----------------------------------------------------------------------------
// ByteArraySink: append std::vector<uint8_t> payloads into a massive
// on-disk byte array (sequential append, large-file optimized).
// -----------------------------------------------------------------------------
// - Policies: Alignment, Preallocation, Fadvise (sequential).
// - Integrates with compile-time multi-sink pipeline via sink_payload<>.
// -----------------------------------------------------------------------------

#include <cstdint>
#include <string>
#include <vector>

#if defined(__linux__) || defined(__APPLE__) || defined(__unix__)
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#define MSINK_POSIX 1
#else
#define MSINK_POSIX 0
#endif

#include "sink_manager.hpp"
#include "sink_policies.hpp"

// ---------------------------------------------------------------------------
// Sink: ByteArraySink
// ---------------------------------------------------------------------------
//
// Appends std::vector<uint8_t> payloads back-to-back into a single file,
// with optional alignment (e.g., 4KiB or 2MiB), preallocation, and fadvise.
//
// The output behaves like a gigantic std::vector<uint8_t> on disk: trivial to
// mmap later and do random access.
//
// Integrates with the compile-time pipeline via sink_payload specialization
// mapping this sink to std::vector<uint8_t>.
template <class AlignmentPolicy = NoAlignment,
          class PreallocPolicy = NoPreallocate,
          class FadvisePolicy = FadviseSequential>
class ByteArraySink {
public:
  using Payload = std::vector<std::uint8_t>;

  explicit ByteArraySink(const std::string &path, std::uint32_t mode = 0644)
      : path_(path), mode_(mode) {
    fd_ = io::open_writable_posix(path_, mode_);
    PreallocPolicy::on_open(fd_);
    FadvisePolicy::on_open(fd_);
  }

  ~ByteArraySink() {
    if (fd_ >= 0)
      io::close_posix(fd_);
  }

  void write(const Payload &bytes) {
    // Optional chunk alignment
    AlignmentPolicy::before_write(fd_, offset_, [&](uint64_t pad){ write_zeros_(pad); });

    if (!bytes.empty()) {
      io::write_all(fd_, bytes.data(), bytes.size());
      offset_ += bytes.size();
    }
  }

  void finalize() {
#if MSINK_POSIX
    ::fdatasync(fd_);
#endif
  }

  std::uint64_t size_bytes() const noexcept { return offset_; }

private:
  void write_zeros_(std::uint64_t count) {
    static constexpr std::size_t ZB = 256 * 1024;
    static const std::uint8_t Z[ZB] = {0};
    std::uint64_t left = count;
    while (left) {
      const std::size_t n = left > ZB ? ZB : static_cast<std::size_t>(left);
      io::write_all(fd_, Z, n);
      offset_ += n;
      left -= n;
    }
  }

  std::string path_;
  std::uint32_t mode_ = 0644;
  int fd_ = -1;
  std::uint64_t offset_ = 0;
};

// ---------------------------------------------------------------------------
// Multi-sink integration: sink_payload specialization
// ---------------------------------------------------------------------------

template <class A, class P, class F>
struct sink_payload<ByteArraySink<A, P, F>> {
  using type = std::vector<std::uint8_t>;
};
