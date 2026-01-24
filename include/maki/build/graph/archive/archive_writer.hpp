
#pragma once
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "maki/core/graph/archive_format.hpp"
#include "sink_manager.hpp"
#include "sink_policies.hpp"

struct ArchivePayload {
  std::vector<uint8_t> bytes;
  std::uint64_t elem_count;
  std::uint8_t bit_width;
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
template <class AlignmentPolicy = AlignTo<4096>,
          class PreallocPolicy = NoPreallocate,
          class FadvisePolicy = FadviseSequential, class ChecksumPolicy = CRC32>
class ArchiveWriter {
public:
  explicit ArchiveWriter(const std::string &path, std::uint32_t mode = 0644,
                         std::uint32_t version = VERSION,
                         std::uint32_t magic = MAGIC_CHAR,
                         PreallocPolicy prealloc = PreallocPolicy{})
      : version_(version), magic_(magic), prealloc_(prealloc) {
    fd_ = io::open_writable_posix(path, mode);
    // File-level policies
    prealloc_.on_open(fd_);
    FadvisePolicy::on_open(fd_);
    // Reserve small TOC by default; user can call reserve_chunks()
    toc_.reserve(1024);
  }

  ~ArchiveWriter() {
    if (fd_ >= 0)
      io::close_posix(fd_);
  }

  // Optional: reserve TOC capacity upfront for scalability (millions of chunks)
  void reserve_chunks(std::size_t n) { toc_.reserve(n); }

  // Append one RAW_PACKED chunk (packed bytes only).
  // Returns chunk id (0-based).
  std::size_t append_raw_packed(const std::uint8_t *data, std::size_t len,
                                std::uint64_t elem_count,
                                std::uint8_t bit_width) {
    if (len && !data)
      throw std::invalid_argument(
          "append_raw_packed: null data with nonzero len");

    // Align (if requested) before writing the chunk
    AlignmentPolicy::before_write(
        fd_, offset_, [&](std::uint64_t pad) { write_zeros_(pad); });

    const std::uint64_t chunk_off = offset_;
    if (len) {
      io::write_all(fd_, data, len);
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
    meta.alignment = static_cast<std::uint16_t>(
        AlignmentPolicy::value <= 0xFFFF ? AlignmentPolicy::value : 0);
    meta.crc32 = csum;

    toc_.push_back(meta);
    elems_prefix_ += elem_count;

    return toc_.size() - 1;
  }

  inline std::size_t write(const ArchivePayload &pld) {
    return append_raw_packed(pld.bytes.data(), pld.bytes.size(), pld.elem_count,
                             pld.bit_width);
  }

  // Finalize: write TOC and footer
  void finalize() {
    if (finalized_)
      return;

    const std::uint64_t toc_off = offset_;
    if (!toc_.empty()) {
      // Write TOC as contiguous array
      io::write_all(fd_, toc_.data(), toc_.size() * sizeof(ChunkMeta));
      offset_ += toc_.size() * sizeof(ChunkMeta);
    }

    Footer f{};
    f.toc_offset = toc_off;
    f.chunk_count = toc_.size();
    f.version = version_;
    f.magic = magic_;

    io::write_all(fd_, &f, sizeof(Footer));
    offset_ += sizeof(Footer);

    // Flush metadata for durability
    ::fdatasync(fd_);
    finalized_ = true;
  }

private:
  void write_zeros_(std::uint64_t bytes) {
    static constexpr std::size_t ZB =
        256 * 1024; // large zero buffer for throughput
    static const std::uint8_t Z[ZB] = {0};

    std::uint64_t left = bytes;
    while (left) {
      const std::size_t n = left > ZB ? ZB : static_cast<std::size_t>(left);
      io::write_all(fd_, Z, n);
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
struct sink_payload<ArchiveWriter<A, P, F, C>> {
  using type = ArchivePayload;
};
