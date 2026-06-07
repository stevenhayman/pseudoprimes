#include "pseudoprime/modular_integer.hpp"

#include <gtest/gtest.h>

namespace {

using pseudoprime::Integer;
using pseudoprime::ModularInteger;

TEST(ModularInteger, DefaultIsIdentity) {
  ModularInteger m(Integer(7));
  EXPECT_EQ(m.value(), Integer(1));
  EXPECT_EQ(m.modulus(), Integer(7));
}

TEST(ModularInteger, NormalizesOnConstruction) {
  EXPECT_EQ(ModularInteger(Integer(10), Integer(7)).value(), Integer(3));
  EXPECT_EQ(ModularInteger(Integer(7), Integer(7)).value(), Integer(0));
  EXPECT_EQ(ModularInteger(Integer(-1), Integer(7)).value(), Integer(6));   // negative wraps up
  EXPECT_EQ(ModularInteger(Integer(-10), Integer(7)).value(), Integer(4));
}

TEST(ModularInteger, MultiplyByIntegerWraps) {
  ModularInteger m(Integer(7));
  m *= Integer(4);
  m *= Integer(5);  // 1 * 4 * 5 = 20 == 6 (mod 7)
  EXPECT_EQ(m.value(), Integer(6));
}

TEST(ModularInteger, MultiplyByNegativeInteger) {
  ModularInteger m(Integer(2), Integer(7));
  m *= Integer(-3);  // 2 * (-3) = -6 == 1 (mod 7)
  EXPECT_EQ(m.value(), Integer(1));
}

TEST(ModularInteger, MultiplyByModularAndFreeOperator) {
  ModularInteger a(Integer(3), Integer(7));
  ModularInteger b(Integer(5), Integer(7));
  ModularInteger c = a * b;  // 15 == 1 (mod 7)
  EXPECT_EQ(c.value(), Integer(1));
  a *= b;
  EXPECT_EQ(a.value(), Integer(1));
}

TEST(ModularInteger, EqualityReducesRhs) {
  ModularInteger m(Integer(3), Integer(7));
  EXPECT_TRUE(m == Integer(3));
  EXPECT_TRUE(m == Integer(10));  // 10 == 3
  EXPECT_TRUE(m == Integer(-4));  // -4 == 3
  EXPECT_FALSE(m == Integer(2));
}

}  // namespace
