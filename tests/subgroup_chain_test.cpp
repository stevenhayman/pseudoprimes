#include "pseudoprime/large_carmichael/subgroup_chain.hpp"

#include <gtest/gtest.h>

namespace {

using pseudoprime::Integer;

TEST(PrimePowerChain, PeelsModuliInOrder) {
  // ghat = 12 = 2^2 * 3  ->  moduli M_1..M_t = 2, 4, 12.
  const pseudoprime::Factorization ghat({{Integer(2), 2}, {Integer(3), 1}});
  const auto chain = pseudoprime::PrimePowerChainStrategy{}.build(ghat);
  ASSERT_EQ(chain->levels(), 3u);

  // coordinate(a, i) = a mod M_{i+1}.
  EXPECT_EQ(chain->coordinate(Integer(7), 0), Integer(1));  // 7 mod 2
  EXPECT_EQ(chain->coordinate(Integer(7), 1), Integer(3));  // 7 mod 4
  EXPECT_EQ(chain->coordinate(Integer(7), 2), Integer(7));  // 7 mod 12
}

TEST(PrimePowerChain, InverseKeyIsModularInverse) {
  const pseudoprime::Factorization ghat({{Integer(2), 2}, {Integer(3), 1}});
  const auto chain = pseudoprime::PrimePowerChainStrategy{}.build(ghat);

  // 3 * 3 == 1 (mod 4); 7 * 7 == 1 (mod 12); 5 * 5 == 1 (mod 12).
  EXPECT_EQ(chain->inverse_key(Integer(3), 1), Integer(3));
  EXPECT_EQ(chain->inverse_key(Integer(7), 2), Integer(7));
  EXPECT_EQ(chain->inverse_key(Integer(5), 2), Integer(5));
}

TEST(BalancedCongruenceChain, CoarserThanPrimePower) {
  // ghat = 2^6 * 3^2 * 5 * 7: prime-power peels every factor (many levels);
  // the balanced chain groups them into ~2^l-sized steps (fewer levels).
  const pseudoprime::Factorization ghat(
      {{Integer(2), 6}, {Integer(3), 2}, {Integer(5), 1}, {Integer(7), 1}});
  const auto prime_power = pseudoprime::PrimePowerChainStrategy{}.build(ghat);
  const auto balanced = pseudoprime::BalancedCongruenceChainStrategy{}.build(ghat);
  EXPECT_GE(balanced->levels(), 1u);
  EXPECT_LT(balanced->levels(), prime_power->levels());
}

TEST(BalancedCongruenceChain, IdentityAndInverseKeyAreConsistent) {
  const pseudoprime::Factorization ghat(
      {{Integer(2), 6}, {Integer(3), 2}, {Integer(5), 1}, {Integer(7), 1}});
  const auto chain = pseudoprime::BalancedCongruenceChainStrategy{}.build(ghat);
  for (std::size_t i = 0; i < chain->levels(); ++i) {
    EXPECT_EQ(chain->coordinate(Integer(1), i), Integer(1));  // identity is the identity coset
    const Integer c = chain->coordinate(Integer(11), i);      // 11 is a unit mod ghat
    EXPECT_EQ(chain->inverse_key(chain->inverse_key(c, i), i), c);  // inverse is an involution
  }
}

}  // namespace
