
#include "maki/build/graph/archive/sdsl_writer.hpp"
#include "test_common.hpp"
#include <gtest/gtest.h>
#include <sdsl/int_vector_mapper.hpp>

template <uint8_t w> using OnDisk = SdslIntVectorOnDiskSink<w>;
template <uint8_t w> using InMem = SdslIntVectorInMemorySink<w>;

template <uint8_t w>
static sdsl::int_vector<w> MakeIV(std::vector<uint64_t> values) {
  sdsl::int_vector<w> iv(values.size(), 0);
  for (size_t i = 0; i < values.size(); ++i)
    iv[i] = values[i];
  return iv;
}

class SdslSinks : public testing::Test {
protected:
  SdslSinks() : path(MakeTempPath(".sdsl")) {}
  ~SdslSinks() { std::filesystem::remove(path); }
  std::string path;
};

TEST_F(SdslSinks, DiskIO) {
  OnDisk<12> sink(path, 0644, 2 << 20);
  auto ivA = MakeIV<12>(GenVector<uint64_t>(1000, (1u << 12) - 1)); // width 12
  auto ivB = MakeIV<12>(GenVector<uint64_t>(1200, (1u << 10) - 1)); // width 10 fits <= 12

  sink.write(ivA);
  sink.write(ivB);
  sink.finalize();

  sdsl::int_vector_mapper<12> map(path);
  ASSERT_EQ(map.size(), ivA.size() + ivB.size());
  EXPECT_EQ(map.width(), 12u);

  // spot-check a few
  EXPECT_EQ(map[0], ivA[0]);
  EXPECT_EQ(map[ivA.size() - 1], ivA[ivA.size() - 1]);
  EXPECT_EQ(map[ivA.size()], ivB[0]);
}

TEST_F(SdslSinks, RAMIO) {
  InMem<1> sink;
  auto ivA = MakeIV<1>(GenVector<uint64_t>(1000, 0b1));
  auto ivB = MakeIV<1>(GenVector<uint64_t>(1200, 0b1));

  sink.write(ivA);
  sink.write(ivB);
  sink.finalize();
  auto &res = sink.data();

  ASSERT_EQ(res.size(), ivA.size() + ivB.size());
  EXPECT_EQ(res.width(), 1);

  // check values
  for (size_t i = 0; i < ivA.size(); ++i) {
    ASSERT_EQ(res[i], ivA[i]);
  }
  for (size_t i = 0; i < ivB.size(); ++i) {
    ASSERT_EQ(res[i + ivA.size()], ivB[i]);
  }
}
