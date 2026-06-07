#include "pseudoprime/large_carmichael/weight.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace {

using pseudoprime::Integer;

// L = 1260 = 2^2 * 3^2 * 5 * 7  ->  components Q = [4, 9, 5, 7].
pseudoprime::Factorization make_1260() {
  return pseudoprime::Factorization({{Integer(2), 2}, {Integer(3), 2}, {Integer(5), 1}, {Integer(7), 1}});
}

TEST(Weight, OmegaByDefinition) {
  const auto c = make_1260().prime_powers();
  EXPECT_EQ(pseudoprime::weight(Integer(1), c), 0u);    // 1 mod every Q_i
  EXPECT_EQ(pseudoprime::weight(Integer(631), c), 1u);  // != 1 only mod Q_1 = 4
  EXPECT_EQ(pseudoprime::weight(Integer(281), c), 2u);  // top nontrivial component is Q_2 = 9
  EXPECT_EQ(pseudoprime::weight(Integer(181), c), 4u);  // != 1 mod Q_4 = 7
}

}  // namespace
