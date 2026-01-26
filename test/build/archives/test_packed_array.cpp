
#include "maki/build/graph/archive/array_builder.hpp"
#include "test_common.hpp"

TEST(PackedArrayBasic, SimpleVector) {
  ArrayBuilder arr;
  arr.reserve(3);
  arr.push(0);
  arr.push(1);
  arr.push(2);

  auto pkd = PackRaw(arr);
  auto view = pkd.view();

  ASSERT_EQ(view.size(), arr.size());
  ASSERT_EQ(view.bit_width(), 2);
  ASSERT_EQ(view.get(0), 0);
  ASSERT_EQ(view.get(1), 1);
  ASSERT_EQ(view.get(2), 2);
}

TEST(PackedArrayBasic, RandomValues) {
  std::size_t n = 100;
  auto arr = GenValues(n, (1u << 10) - 1);
  auto packedArr = PackRaw(arr);
  auto v = packedArr.view();

  ASSERT_EQ(v.size(), arr.size());
  ASSERT_LE(v.bit_width(), 10);
  for (size_t i = 0; i < n; ++i) {
    ASSERT_EQ(v.get(i), arr.get(i));
  }
}
