
#include "maki/build/graph/archive/byte_writer.hpp"
#include "maki/build/graph/archive/vector_writer.hpp"
#include "test_common.hpp"
#include <gtest/gtest.h>


class ByteSink : public testing::Test {
protected:
  ByteSink() : path(MakeTempPath(".char")) {}
  ~ByteSink() { std::filesystem::remove(path); }
  std::string path;
};

TEST_F(ByteSink, DiskIO) {
  ByteArraySink sink(path);

  auto vA = GenVector<uint8_t>(1000, 0xFF);
  auto vB = GenVector<uint8_t>(1200, 0xFF);

  sink.write(vA);
  sink.write(vB);
  sink.finalize();

  auto raw = ReadFile(path);

  // check values
  for (size_t i = 0; i < vA.size(); ++i) {
    ASSERT_EQ(raw[i], vA[i]);
  }
  for (size_t i = 0; i < vB.size(); ++i) {
    ASSERT_EQ(raw[i + vA.size()], vB[i]);
  }
}

TEST(VectorSink, RandomVector) {
  VectorInMemorySink<uint64_t> sink;
  
  auto vA = GenVector(1000, (std::size_t)(1ull << 12) - 1);
  auto vB = GenVector(1200, (std::size_t)(1ull << 32) - 1);
  auto vC = GenVector(1100, (std::size_t)-1);

  sink.write(vA);
  sink.write(vB);
  sink.write(vC);
  sink.finalize();

  auto &res = sink.data();
  for (size_t i = 0; i < vA.size(); ++i) {
    ASSERT_EQ(res[i], vA[i]);
  }
  std::size_t rtot = vA.size();
  for (size_t i = 0; i < vB.size(); ++i) {
    ASSERT_EQ(res[i + rtot], vB[i]);
  }
  rtot += vB.size();
  for (size_t i = 0; i < vC.size(); ++i) {
    ASSERT_EQ(res[i + rtot], vC[i]);
  }
}
