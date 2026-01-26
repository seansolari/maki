
#include "maki/build/graph/archive/array_builder.hpp"
#include "test_common.hpp"
#include <filesystem>
#include <gtest/gtest.h>

// Include your policy-based writer/reader
#include "maki/build/graph/archive/archive_writer.hpp"
#include "maki/core/graph/archive/archive_reader.hpp"

class ArchiveBasic : public testing::Test {
protected:
  ArchiveBasic() : path(MakeTempPath(".char")) {}
  ~ArchiveBasic() { std::filesystem::remove(path); }
  std::string path;
};

TEST_F(ArchiveBasic, RoundTripTwoChunksAligned) {
  using Writer =
      ArchiveWriter<AlignTo<4096>, NoPreallocate, FadviseSequential, CRC32>;

  Writer w(path);
  w.reserve_chunks(2);

  // Chunk A: 16-bit values
  auto valsA = GenValues(1'000'0, (1u << 16) - 1);
  auto packedA = PackRaw(valsA);
  auto idA = w.append_raw_packed(packedA.data.data(), packedA.data.size(),
                                 packedA.elem_count, packedA.bit_width);

  // Chunk B: 20-bit values
  auto valsB = GenValues(12'345, (1u << 20) - 1);
  auto packedB = PackRaw(valsB);
  auto idB = w.append_raw_packed(packedB.data.data(), packedB.data.size(),
                                 packedB.elem_count, packedB.bit_width);

  EXPECT_EQ(idA, 0u);
  EXPECT_EQ(idB, 1u);
  w.finalize();

  ArchiveReader r(path);
  ASSERT_EQ(r.chunk_count(), 2u);

  // Check offsets & alignment
  const auto &m0 = r.meta(0);
  const auto &m1 = r.meta(1);

  EXPECT_EQ(m0.file_offset % 4096, 0u);
  EXPECT_EQ(m1.file_offset % 4096, 0u);

  EXPECT_EQ(m0.byte_len, packedA.data.size());
  EXPECT_EQ(m1.byte_len, packedB.data.size());

  // CRCs match
  auto [p0, l0] = r.chunk_region(0);
  auto [p1, l1] = r.chunk_region(1);
  EXPECT_EQ(crc32_compute(p0, l0), m0.crc32);
  EXPECT_EQ(crc32_compute(p1, l1), m1.crc32);

  // total_elems and index mapping
  EXPECT_EQ(r.total_elems(), m0.start_index + m0.elem_count + m1.elem_count);

  // Any global index in chunk 1 maps to 1
  uint64_t probe = m1.start_index + (m1.elem_count ? (m1.elem_count - 1) : 0);
  long k = r.find_chunk_by_global_index(probe);
  EXPECT_EQ(k, 1);
}

TEST_F(ArchiveBasic, NoAlignmentContiguousOffsets) {
  using Writer =
      ArchiveWriter<NoAlignment, NoPreallocate, FadviseSequential, CRC32>;

  Writer w(path);
  w.reserve_chunks(3);

  auto A = PackRaw(GenValues(1000, (1u << 8) - 1));
  auto B = PackRaw(GenValues(2000, (1u << 10) - 1));
  auto C = PackRaw(GenValues(3000, (1u << 7) - 1));

  w.append_raw_packed(A.data.data(), A.data.size(), A.elem_count, A.bit_width);
  w.append_raw_packed(B.data.data(), B.data.size(), B.elem_count, B.bit_width);
  w.append_raw_packed(C.data.data(), C.data.size(), C.elem_count, C.bit_width);
  w.finalize();

  ArchiveReader r(path);
  ASSERT_EQ(r.chunk_count(), 3u);

  const auto &m0 = r.meta(0);
  const auto &m1 = r.meta(1);
  const auto &m2 = r.meta(2);

  EXPECT_EQ(m1.file_offset, m0.file_offset + m0.byte_len);
  EXPECT_EQ(m2.file_offset, m1.file_offset + m1.byte_len);
}

TEST_F(ArchiveBasic, ReadChunk) {
  using Writer =
      ArchiveWriter<NoAlignment, NoPreallocate, FadviseSequential, CRC32>;

  Writer w(path);

  auto A = PackRaw(GenValues(1000, (1u << 8) - 1));
  auto B = PackRaw(GenValues(2000, (1u << 10) - 1));
  auto C = PackRaw(GenValues(3000, (1u << 7) - 1));

  w.append_raw_packed(A.data.data(), A.data.size(), A.elem_count, A.bit_width);
  w.append_raw_packed(B.data.data(), B.data.size(), B.elem_count, B.bit_width);
  w.append_raw_packed(C.data.data(), C.data.size(), C.elem_count, C.bit_width);
  w.finalize();

  ArchiveReader r(path);
  auto view = r.view(2);
  auto Cview = C.view();
  ASSERT_EQ(view.size(), Cview.size());
  for (std::size_t i = 0; i < view.size(); ++i) {
    ASSERT_EQ(view.get(i), Cview.get(i));
  }  
}
