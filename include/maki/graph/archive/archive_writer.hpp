
#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <stdexcept>
#include <system_error>
#include <type_traits>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "archive_format.hpp"
#include "crc32.hpp"
#include "sink_manager.hpp"

// ============================ Policies ======================================

struct NoAlignment
{
  static constexpr std::uint32_t value = 0;
  template <class WriteZerosFn>
  static void before_write(int /*fd*/, std::uint64_t /*offset*/, WriteZerosFn && /*wz*/) {}
};

template <std::uint32_t AlignBytes>
struct AlignTo
{
  static_assert((AlignBytes & (AlignBytes - 1)) == 0, "AlignTo<N>: N must be power of two");
  static constexpr std::uint32_t value = AlignBytes;
  template <class WriteZerosFn>
  static void before_write(int /*fd*/, std::uint64_t current_offset, WriteZerosFn &&write_zeros)
  {
    if constexpr (AlignBytes)
    {
      const std::uint64_t mis = current_offset & (AlignBytes - 1);
      if (mis)
        write_zeros(AlignBytes - mis);
    }
  }
};

struct NoPreallocate
{
  static void on_open(int /*fd*/) {}
};

template <std::uint64_t Bytes>
struct Preallocate
{
  static void on_open(int fd)
  {
    // Best-effort; ignore ENOSYS/EOPNOTSUPP
    (void)::posix_fallocate(fd, 0, static_cast<off_t>(Bytes));
  }
};

template <>
struct Preallocate<0>
{
  std::uint64_t bytes;

  void on_open(int fd)
  {
    // Best-effort; ignore ENOSYS/EOPNOTSUPP
    (void)::posix_fallocate(fd, 0, static_cast<off_t>(bytes));
  }
};

struct NoFadvise
{
  static void on_open(int /*fd*/) {}
};

struct FadviseSequential
{
  static void on_open(int fd)
  {
    (void)::posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
  }
};

// CRC32 policy
struct CRC32
{
  static std::uint32_t compute(const std::uint8_t *data, std::size_t len)
  {
    return crc32_compute(data, len);
  }
};

struct NoChecksum
{
  static std::uint32_t compute(const std::uint8_t * /*data*/, std::size_t /*len*/)
  {
    return 0;
  }
};

// ========================== POSIX helpers ===================================

namespace detail
{

  inline int open_writable_posix(const std::string &path, std::uint32_t mode)
  {
    int fd = ::open(path.c_str(), O_CREAT | O_TRUNC | O_WRONLY, mode);
    if (fd < 0)
      throw std::system_error(errno, std::generic_category(), "open");
    return fd;
  }

  inline void close_posix(int fd)
  {
    ::close(fd);
  }

  inline void write_all(int fd, const void *buf, std::size_t len)
  {
    const std::uint8_t *p = static_cast<const std::uint8_t *>(buf);
    std::size_t left = len;
    while (left > 0)
    {
      ssize_t w = ::write(fd, p, left);
      if (w < 0)
        throw std::system_error(errno, std::generic_category(), "write");
      p += static_cast<std::size_t>(w);
      left -= static_cast<std::size_t>(w);
    }
  }

} // namespace detail


/**
 * Compressed data to be inserted into the archive.
 */
class ArchivePayload
{
  std::vector<uint8_t> bytes;
  uint64_t elem_count;
  uint64_t bit_width;
};

// ========================= ArchiveWriter ====================================
//
// Template parameters (policies):
//  - AlignmentPolicy: AlignTo<4096>, AlignTo<2*1024*1024>, or NoAlignment
//  - PreallocPolicy:  Preallocate<bytes> or NoPreallocate
//  - FadvisePolicy:   FadviseSequential or NoFadvise
//  - ChecksumPolicy:  CRC32 or NoChecksum
//
// Optimized for large sequential writes (hundreds of MB to multi-GB).
//
template <
    class AlignmentPolicy = AlignTo<4096>,
    class PreallocPolicy = NoPreallocate,
    class FadvisePolicy = FadviseSequential,
    class ChecksumPolicy = CRC32>
class ArchiveWriter
{
public:
  explicit ArchiveWriter(const std::string &path, std::uint32_t mode = 0644,
                         std::uint32_t version = VERSION, std::uint32_t magic = MAGIC_CHAR,
                         PreallocPolicy prealloc = PreallocPolicy{})
      : version_(version), magic_(magic), prealloc_(prealloc)
  {
    fd_ = detail::open_writable_posix(path, mode);
    // File-level policies
    prealloc_.on_open(fd_);
    FadvisePolicy::on_open(fd_);
    // Reserve small TOC by default; user can call reserve_chunks()
    toc_.reserve(1024);
  }

  ~ArchiveWriter()
  {
    if (fd_ >= 0)
      detail::close_posix(fd_);
  }

  // Optional: reserve TOC capacity upfront for scalability (millions of chunks)
  void reserve_chunks(std::size_t n) { toc_.reserve(n); }

  // Append one RAW_PACKED chunk (packed bytes only).
  // Returns chunk id (0-based).
  std::size_t append_raw_packed(const std::uint8_t *data,
                                std::size_t len,
                                std::uint64_t elem_count,
                                std::uint8_t bit_width)
  {
    if (len && !data)
      throw std::invalid_argument("append_raw_packed: null data with nonzero len");

    // Align (if requested) before writing the chunk
    AlignmentPolicy::before_write(fd_, offset_, &{ write_zeros_(pad); });

    const std::uint64_t chunk_off = offset_;
    if (len)
    {
      detail::write_all(fd_, data, len);
      offset_ += len;
    }

    // Compute checksum via policy (CRC32 or NoChecksum)
    const std::uint32_t csum = ChecksumPolicy::compute(data, len);

    ChunkMeta meta{};
    meta.start_index = elems_prefix_;
    meta.elem_count = elem_count;
    meta.file_offset = chunk_off;
    meta.byte_len = len;
    meta.bit_width = bit_width;
    meta.encoding = static_cast<std::uint8_t>(Encoding::RAW_PACKED);
    meta.alignment = static_cast<std::uint16_t>(AlignmentPolicy::value <= 0xFFFF ? AlignmentPolicy::value : 0);
    meta.crc32 = csum;

    toc_.push_back(meta);
    elems_prefix_ += elem_count;

    return toc_.size() - 1;
  }

  inline std::size_t write(const ArchivePayload &pld)
  {
    return append_raw_packed(pld.bytes.data(), pld.bytes.size(), pld.elem_count, pld.bit_width);
  }

  // Finalize: write TOC and footer
  void finalize()
  {
    if (finalized_)
      return;

    const std::uint64_t toc_off = offset_;
    if (!toc_.empty())
    {
      // Write TOC as contiguous array
      detail::write_all(fd_, toc_.data(), toc_.size() * sizeof(ChunkMeta));
      offset_ += toc_.size() * sizeof(ChunkMeta);
    }

    Footer f{};
    f.toc_offset = toc_off;
    f.chunk_count = toc_.size();
    f.version = version_;
    f.magic = magic_;

    detail::write_all(fd_, &f, sizeof(Footer));
    offset_ += sizeof(Footer);

    // Flush metadata for durability
    ::fdatasync(fd_);
    finalized_ = true;
  }

private:
  void write_zeros_(std::uint64_t bytes)
  {
    static constexpr std::size_t ZB = 256 * 1024; // large zero buffer for throughput
    static const std::uint8_t Z[ZB] = {0};

    std::uint64_t left = bytes;
    while (left)
    {
      const std::size_t n = left > ZB ? ZB : static_cast<std::size_t>(left);
      detail::write_all(fd_, Z, n);
      offset_ += n;
      left -= n;
    }
  }

  int fd_ = -1;
  std::uint64_t offset_ = 0;
  std::uint64_t elems_prefix_ = 0;
  bool finalized_ = false;

  std::uint32_t version_;
  std::uint32_t magic_;
  PreallocPolicy prealloc_;

  std::vector<ChunkMeta> toc_;
};

template <class A, class P, class F, class C>
struct sink_payload<ArchiveWriter<A,P,F,C>> { using type = ArchivePayload; };
