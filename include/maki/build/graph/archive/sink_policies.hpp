
#pragma once
#include <cstdint>
#include <string>
#include <system_error>

#if defined(__linux__) || defined(__APPLE__) || defined(__unix__)
#define MSINK_POSIX 1
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#else
#define MSINK_POSIX 0
#endif

#include "maki/core/utils/crc32.hpp"

// --------------------------- Alignment Policy -------------------------------

struct NoAlignment {
  static constexpr std::uint32_t value = 0;
  template <class WriteZerosFn>
  static void before_write(int /*fd*/, std::uint64_t /*offset*/,
                           WriteZerosFn && /*write_zeros*/) {}
};

template <std::uint32_t AlignBytes> struct AlignTo {
  static_assert((AlignBytes & (AlignBytes - 1)) == 0,
                "AlignTo<N>: N must be power of two");
  static constexpr std::uint32_t value = AlignBytes;
  template <class WriteZerosFn>
  static void before_write(int /*fd*/, std::uint64_t current_offset,
                           WriteZerosFn &&write_zeros) {
    if constexpr (AlignBytes > 0) {
      const std::uint64_t mis = current_offset & (AlignBytes - 1);
      if (mis)
        write_zeros(AlignBytes - mis);
    }
  }
};

// --------------------------- Preallocation Policy ---------------------------

struct NoPreallocate {
  static void on_open(int /*fd*/) {}
};

template <std::uint64_t Bytes> struct Preallocate {
  static void on_open(int fd) {
#if MSINK_POSIX
    (void)::posix_fallocate(fd, 0, static_cast<off_t>(Bytes)); // best-effort
#else
    (void)fd;
#endif
  }
};

template <> struct Preallocate<0> {
  std::uint64_t bytes = 0;

  void on_open(int fd) {
#if MSINK_POSIX
    (void)::posix_fallocate(fd, 0, static_cast<off_t>(bytes)); // best-effort
#else
    (void)fd;
#endif
  }
};

// --------------------------- Fadvise Policy ---------------------------------

struct NoFadvise {
  static void on_open(int /*fd*/) {}
};

struct FadviseSequential {
  static void on_open(int fd) {
#if MSINK_POSIX
    (void)::posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#else
    (void)fd;
#endif
  }
};

// --------------------------- Fadvise Policy ---------------------------------

struct CRC32 {
  static std::uint32_t compute(const std::uint8_t *data, std::size_t len) {
    return crc32_compute(data, len);
  }
};

struct NoChecksum {
  static std::uint32_t compute(const std::uint8_t * /*data*/,
                               std::size_t /*len*/) {
    return 0;
  }
};

// --------------------------- POSIX I/O Helpers ------------------------------

namespace io {

inline int open_writable_posix(const std::string &path, std::uint32_t mode) {
#if MSINK_POSIX
  int fd = ::open(path.c_str(), O_CREAT | O_TRUNC | O_WRONLY, mode);
  if (fd < 0)
    throw std::system_error(errno, std::generic_category(), "open");
  return fd;
#else
  (void)path;
  (void)mode;
  throw std::runtime_error("POSIX not available");
#endif
}

inline void close_posix(int fd) {
#if MSINK_POSIX
  ::close(fd);
#else
  (void)fd;
#endif
}

inline void write_all(int fd, const void *buf, std::size_t len) {
#if MSINK_POSIX
  const std::uint8_t *p = static_cast<const std::uint8_t *>(buf);
  std::size_t left = len;
  while (left > 0) {
    ssize_t w = ::write(fd, p, left);
    if (w < 0)
      throw std::system_error(errno, std::generic_category(), "write");
    p += static_cast<std::size_t>(w);
    left -= static_cast<std::size_t>(w);
  }
#else
  (void)fd;
  (void)buf;
  (void)len;
  throw std::runtime_error("POSIX not available");
#endif
}

inline void fdatasync_posix(int fd) {
#if MSINK_POSIX
  ::fdatasync(fd);
#else
  (void)fd;
#endif
}

} // namespace io
