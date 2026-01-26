
#pragma once
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <future>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "archive_format.hpp"
#include "array.hpp"
#include "maki/core/utils/crc32.hpp"

// ------------------------------ Reader --------------------------------------

class ArchiveReader {
public:
  explicit ArchiveReader(const std::string &path) {
    // Open
    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0)
      throw std::system_error(errno, std::generic_category(), "open");

    // Size
    struct stat st{};
    if (::fstat(fd_, &st) != 0) {
      int e = errno;
      ::close(fd_);
      throw std::system_error(e, std::generic_category(), "fstat");
    }
    file_size_ = static_cast<std::size_t>(st.st_size);
    if (file_size_ < sizeof(Footer))
      throw std::runtime_error("Archive too small (no footer)");

    // Map
    map_ = static_cast<const std::uint8_t *>(
        ::mmap(nullptr, file_size_, PROT_READ, MAP_SHARED, fd_, 0));
    if (map_ == MAP_FAILED) {
      int e = errno;
      ::close(fd_);
      throw std::system_error(e, std::generic_category(), "mmap");
    }

    // Footer
    const std::size_t footer_off = file_size_ - sizeof(Footer);
    footer_ = reinterpret_cast<const Footer *>(map_ + footer_off);
    if (footer_->magic != MAGIC_CHAR || footer_->version != VERSION)
      throw std::runtime_error("Invalid archive (magic/version mismatch)");

    // TOC bounds
    const std::size_t toc_off = static_cast<std::size_t>(footer_->toc_offset);
    const std::size_t toc_len =
        static_cast<std::size_t>(footer_->chunk_count) * sizeof(ChunkMeta);
    if (toc_off > footer_off || toc_off + toc_len > footer_off)
      throw std::runtime_error("Invalid TOC bounds");

    toc_ = reinterpret_cast<const ChunkMeta *>(map_ + toc_off);
    chunk_count_ = static_cast<std::size_t>(footer_->chunk_count);

    // Total elements
    total_elems_ = 0;
    if (chunk_count_) {
      const auto &last = toc_[chunk_count_ - 1];
      total_elems_ = last.start_index + last.elem_count;
    }

    // Access pattern hint: random across chunks by default
    ::posix_madvise(const_cast<std::uint8_t *>(map_), file_size_,
                    POSIX_MADV_RANDOM);
  }

  ~ArchiveReader() {
    if (map_ && map_ != MAP_FAILED)
      ::munmap(const_cast<std::uint8_t *>(map_), file_size_);
    if (fd_ >= 0)
      ::close(fd_);
  }

  // --------- Basic metadata ---------
  std::size_t chunk_count() const noexcept { return chunk_count_; }
  std::uint64_t total_elems() const noexcept { return total_elems_; }

  const ChunkMeta &meta(std::size_t i) const {
    if (i >= chunk_count_)
      throw std::out_of_range("meta: index OOB");
    return toc_[i];
  }

  // Return (ptr, len) for chunk i’s data region
  std::pair<const std::uint8_t *, std::size_t>
  chunk_region(std::size_t i) const {
    const auto &m = meta(i);
    const std::size_t off = static_cast<std::size_t>(m.file_offset);
    const std::size_t len = static_cast<std::size_t>(m.byte_len);
    const std::size_t max_end =
        file_size_ - sizeof(Footer); // data must end before footer
    if (off > max_end || off + len > max_end)
      throw std::runtime_error("chunk_region: invalid bounds");
    return {map_ + off, len};
  }

  // O(log N) global index -> chunk id (or -1 if not found)
  long find_chunk_by_global_index(std::uint64_t idx) const {
    if (!chunk_count_)
      return -1;
    std::size_t lo = 0, hi = chunk_count_;
    while (lo < hi) {
      std::size_t mid = (lo + hi) / 2;
      if (toc_[mid].start_index <= idx)
        lo = mid + 1;
      else
        hi = mid;
    }
    if (lo == 0)
      return -1;
    std::size_t k = lo - 1;
    const auto &m = toc_[k];
    if (idx < m.start_index + m.elem_count)
      return static_cast<long>(k);
    return -1;
  }

  // --------- Read-only view ---------

  // Obtain array view on chunk i
  RawPackedArray view(std::uint64_t i) const {
    const auto &m = meta(i);
    const std::size_t off = static_cast<std::size_t>(m.file_offset);
    const std::size_t len = static_cast<std::size_t>(m.byte_len);
    const std::size_t max_end =
        file_size_ - sizeof(Footer); // data must end before footer
    if (off > max_end || off + len > max_end)
      throw std::runtime_error("chunk_region: invalid bounds");
    return RawPackedArray{map_ + off, len, m.elem_count, m.bit_width};
  }

  // --------- Checksum validation APIs ---------

  struct ValidationResult {
    // indices of chunks whose CRC did not match
    std::vector<std::size_t> mismatches;
    // total chunks validated
    std::size_t validated = 0;
  };

  // Validate a single chunk’s CRC32; returns true if OK.
  bool validate_chunk(std::size_t i) const {
    const auto &m = meta(i);
    auto [ptr, len] = chunk_region(i);
    const std::uint32_t crc = crc32_compute(ptr, len);
    return crc == m.crc32;
  }

  // Validate [first, last) chunk indices (single-thread)
  ValidationResult validate_range(std::size_t first, std::size_t last) const {
    if (first > last || last > chunk_count_)
      throw std::out_of_range("validate_range OOB");
    ValidationResult r;
    r.mismatches.reserve(last - first);
    for (std::size_t i = first; i < last; ++i) {
      if (!validate_chunk(i))
        r.mismatches.push_back(i);
    }
    r.validated = last - first;
    return r;
  }

  // Validate all chunks (single-thread)
  ValidationResult validate_all() const {
    return validate_range(0, chunk_count_);
  }

  // Validate all chunks in parallel. The file is read-only mmapped,
  // so this is safe across threads. Choose threads to match I/O and CPU.
  ValidationResult validate_all_parallel(
      unsigned threads = std::thread::hardware_concurrency()) const {
    if (chunk_count_ == 0)
      return {};
    if (threads == 0)
      threads = 1;

    const std::size_t N = chunk_count_;
    const std::size_t per = (N + threads - 1) / threads;

    std::vector<std::future<ValidationResult>> futs;
    futs.reserve(threads);

    for (unsigned t = 0; t < threads; ++t) {
      const std::size_t start = t * per;
      if (start >= N)
        break;
      const std::size_t end = std::min(N, start + per);

      futs.emplace_back(std::async(std::launch::async, [this, start, end] {
        // Hint OS that this subrange may be scanned sequentially
        try {
          // Calling madvise per chunk is expensive; do nothing here by default.
          // If you have very large contiguous ranges, you can add:
          // auto [ptr, len] = chunk_region(start);
          // ::posix_madvise(const_cast<std::uint8_t*>(ptr), (end-start ?
          // (chunk_region(end-1).first+chunk_region(end-1).second - ptr) : 0),
          // POSIX_MADV_SEQUENTIAL);
        } catch (...) {
          // ignore hints errors
        }
        return this->validate_range(start, end);
      }));
    }

    ValidationResult all;
    for (auto &f : futs) {
      auto part = f.get();
      all.validated += part.validated;
      // Concatenate mismatches
      all.mismatches.insert(all.mismatches.end(), part.mismatches.begin(),
                            part.mismatches.end());
    }
    // Keep mismatches sorted
    std::sort(all.mismatches.begin(), all.mismatches.end());
    return all;
  }

  // Recompute and return CRC32 for a chunk (helper)
  std::uint32_t recompute_crc(std::size_t i) const {
    auto [ptr, len] = chunk_region(i);
    return crc32_compute(ptr, len);
  }

private:
  int fd_ = -1;
  std::size_t file_size_ = 0;
  const std::uint8_t *map_ = nullptr;

  const Footer *footer_ = nullptr;
  const ChunkMeta *toc_ = nullptr;
  std::size_t chunk_count_ = 0;
  std::uint64_t total_elems_ = 0;
};
