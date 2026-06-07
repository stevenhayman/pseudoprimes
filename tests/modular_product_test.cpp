#include "pseudoprime/modular_product.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace {

using pseudoprime::Integer;
using pseudoprime::ModularProduct;

TEST(ModularProduct, EmptyIsIdentity) {
  const ModularProduct mp(Integer(7));
  EXPECT_EQ(mp.residue(), Integer(1));
  EXPECT_EQ(mp.modulus(), Integer(7));
  EXPECT_EQ(mp.size(), 0u);
  EXPECT_TRUE(mp.elements().empty());
}

TEST(ModularProduct, MultiplyAccumulatesResidueAndElements) {
  ModularProduct mp(Integer(7));
  mp.multiply(Integer(4), 0);
  mp.multiply(Integer(5), 2);  // 1 * 4 * 5 = 20 == 6 (mod 7)
  EXPECT_EQ(mp.residue(), Integer(6));
  EXPECT_EQ(mp.elements(), (std::vector<std::size_t>{0, 2}));
}

TEST(ModularProduct, CombineMultipliesAndConcatenates) {
  ModularProduct a(Integer(7));
  a.multiply(Integer(3), 1);  // residue 3, elements {1}
  ModularProduct b(Integer(7));
  b.multiply(Integer(5), 4);  // residue 5, elements {4}

  const ModularProduct c = a.combine(b);  // 15 == 1 (mod 7)
  EXPECT_EQ(c.residue(), Integer(1));
  EXPECT_EQ(c.elements(), (std::vector<std::size_t>{1, 4}));
  EXPECT_EQ(a.elements(), (std::vector<std::size_t>{1}));  // operands unchanged
}

}  // namespace
