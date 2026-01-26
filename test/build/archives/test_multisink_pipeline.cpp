
#include "maki/build/graph/archive/byte_writer.hpp"
#include "maki/core/graph/archive/archive_reader.hpp"
#include "test_common.hpp"
#include <chrono>
#include <gtest/gtest.h>
#include <random>
#include <thread>

#include "maki/build/graph/archive/archive_writer.hpp"
#include "maki/build/graph/archive/sink_manager.hpp"

// Type aliases for sinks in this test

using ArchiveSinkT =
    ArchiveWriter<AlignTo<4096>, NoPreallocate, FadviseSequential, CRC32>;
using ByteSinkT = ByteArraySink<>;

using Sinks = std::tuple<ArchiveSinkT, ByteSinkT>;
using Bundle = ChunkBundleT<ArchiveSinkT, ByteSinkT>;
using Pool = BundlePool<ArchiveSinkT, ByteSinkT>;
using Multi = MultiSink<ArchiveSinkT, ByteSinkT>;

struct FakeJob {
  void setPool(std::shared_ptr<Pool> &p_) { p = p_; }

  std::unique_ptr<Bundle> operator()(uint64_t i) const {
    std::mt19937 rng(1234 + i);

    auto bnd = p->acquire();
    bnd->id = i;

    // Fill archive payload
    auto &[pa, pb] = bnd->payloads;

    // generate some values with width 12
    pa.raw = GenValues(1000 + (i * 10), (1u << 12) - 1);
    // Fill byte sink payload
    pb.assign(1024 + (i % 5) * 100, (uint8_t)(i & 0xFF));

    // Random jitter to force out-of-order arrivals
    std::this_thread::sleep_for(std::chrono::milliseconds(rng() % 100));

    return bnd;
  }

  std::shared_ptr<Pool> p;
};

class MultiSinkPipeline : public testing::Test {
protected:
  MultiSinkPipeline()
      : archive_path(MakeTempPath(".char")), bytes_path(MakeTempPath(".u8")) {}
  ~MultiSinkPipeline() {
    std::filesystem::remove(archive_path);
    std::filesystem::remove(bytes_path);
  }
  std::string archive_path, bytes_path;
};

TEST_F(MultiSinkPipeline, ValidSinks) {
  Multi sinks{ArchiveSinkT{archive_path}, ByteSinkT{bytes_path}};
  sinks.finalize();
}

TEST_F(MultiSinkPipeline, OutOfOrderSubmissionsAreOrderedAndFannedOut) {
  std::size_t N = 24;
  // create archive
  {
    Multi sinks{ArchiveSinkT{archive_path}, ByteSinkT{bytes_path}};
    ProcessChunks(FakeJob(), sinks, 4, 1 << 20, N);
    sinks.finalize();
  }

  // Validate archive
  {
    ArchiveReader r(archive_path);
    ASSERT_EQ(r.chunk_count(), N);
    auto chk = r.validate_all_parallel();
    EXPECT_TRUE(chk.mismatches.empty());
    // Ensure TOC is in increasing order by file_offset and start_index
    uint64_t prev_off = 0, prev_idx = 0;
    for (size_t i = 0; i < r.chunk_count(); ++i) {
      const auto &m = r.meta(i);
      if (i > 0) {
        EXPECT_LE(prev_off, m.file_offset);
        EXPECT_LE(prev_idx, m.start_index);
      }
      prev_off = m.file_offset;
      prev_idx = m.start_index;
    }
  }

  // Validate bytes file is pure concatenation of the submitted byte payloads in
  // order
  {
    auto raw = ReadFile(bytes_path);
    size_t pos = 0;
    for (std::size_t i = 0; i < N; ++i) {
      size_t len = 1024 + (i % 5) * 100;
      ASSERT_LE(pos + len, raw.size());
      for (size_t k = 0; k < len; ++k) {
        ASSERT_EQ(raw[pos + k], (uint8_t)(i & 0xFF));
      }
      pos += len;
    }
    ASSERT_EQ(pos, raw.size());
  }
}
