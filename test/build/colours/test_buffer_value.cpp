
#include "maki/build/graph/build_colours.hpp"
#include <cstdint>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using ::testing::ElementsAreArray;

class BufferValuePairTests : public testing::Test {
protected:
  BufferValuePairTests() : _genome(0b11111111), _edge(0b001) {}

  uint64_t _genome, _edge;
};

TEST_F(BufferValuePairTests, FromBytes) {
  uint8_t src_[2] = {0b11111001, 0b00000111};
  BufferValue kv(src_, 2);
  ASSERT_EQ(kv.colour(), _genome);
  ASSERT_EQ(kv.edge(), _edge);
}

TEST_F(BufferValuePairTests, ToBytes) {
  uint8_t dest_[2] = {0, 0};
  BufferValue kv(_genome, _edge);
  kv.flush(dest_, 2);
  ASSERT_THAT(dest_, ElementsAreArray({0b11111001, 0b00000111}));
}
