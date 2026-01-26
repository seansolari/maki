
#include "maki/build/graph/archive/archive_writer.hpp"
#include "maki/core/graph/archive/archive_reader.hpp"
#include "test_common.hpp"
#include <gtest/gtest.h>
#include <vector>

class ArchiveValidation : public testing::Test {
protected:
  ArchiveValidation() : path(MakeTempPath(".char")) {}
  ~ArchiveValidation() { std::filesystem::remove(path); }
  std::string path;
};

TEST_F(ArchiveValidation, ParallelValidationAndCorruption) {
  using Writer =
      ArchiveWriter<AlignTo<4096>, NoPreallocate, FadviseSequential, CRC32>;

  Writer w(path);
  const int N = 8;
  w.reserve_chunks(N);

  std::vector<PackedResult> chunks(N);

  for (int i = 0; i < N; ++i) {
    uint8_t width = 8 + (i % 8); // vary widths
    auto vals = GenValues(5000 + 100 * i, (1ull << width) - 1);
    chunks[i] = PackRaw(vals);
    w.append_raw_packed(chunks[i].data.data(), chunks[i].data.size(),
                        chunks[i].elem_count, chunks[i].bit_width);
  }
  w.finalize();

  ArchiveReader r(path);

  auto all = r.validate_all_parallel();
  EXPECT_EQ(all.validated, (size_t)N);
  EXPECT_TRUE(all.mismatches.empty());

  // Corrupt one byte in middle of chunk 3
  {
    auto raw = ReadFile(path);
    const auto &m3 = r.meta(3);
    size_t flip_pos = (size_t)m3.file_offset + (size_t)m3.byte_len / 2;
    raw[flip_pos] ^= 0xFF;
    WriteFile(path, raw);
  }

  ArchiveReader r2(path);
  auto res = r2.validate_all_parallel();
  ASSERT_EQ(res.validated, (size_t)N);
  ASSERT_FALSE(res.mismatches.empty());
  // Expect chunk 3 to be in mismatches (could be more depending on CRC
  // collisions—but practically no).
  EXPECT_NE(std::find(res.mismatches.begin(), res.mismatches.end(), 3u),
            res.mismatches.end());
}
