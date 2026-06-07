#include "pseudoprime/large_carmichael/construction.hpp"

#include <gtest/gtest.h>
#include <initializer_list>
#include <memory>
#include <utility>
#include <vector>

#include "pseudoprime/carmichael.hpp"
#include "pseudoprime/large_carmichael/subgroup_chain.hpp"

namespace {

using pseudoprime::Integer;

pseudoprime::Factorization L_from(std::initializer_list<std::pair<long, unsigned long>> ps) {
  std::vector<pseudoprime::PrimePower> f;
  for (const auto& pe : ps) f.push_back({Integer(pe.first), pe.second});
  return pseudoprime::Factorization(f);
}

std::vector<pseudoprime::Factorization> candidate_Ls() {
  return {
      L_from({{2, 2}, {3, 2}, {5, 1}, {7, 1}}),          // 1260
      L_from({{2, 3}, {3, 2}, {5, 1}, {7, 1}}),          // 2520
      L_from({{2, 4}, {3, 2}, {5, 1}, {7, 1}}),          // 5040
      L_from({{2, 3}, {3, 2}, {5, 1}, {7, 1}, {11, 1}}), // 27720
      L_from({{2, 4}, {3, 2}, {5, 1}, {7, 1}, {11, 1}, {13, 1}}),
  };
}

// Every returned result must be a genuine Carmichael number whose product is
// congruent to 1 modulo L, and at least one candidate L must succeed.
TEST(LargeCarmichael, ConstructsValidCarmichael) {
  bool any = false;
  for (const auto& L : candidate_Ls()) {
    const auto result = pseudoprime::construct_large_carmichael(L);
    if (!result) continue;
    any = true;
    EXPECT_TRUE(pseudoprime::is_carmichael(*result)) << L.value();
    EXPECT_EQ(result->value() % L.value(), Integer(1)) << L.value();
    EXPECT_GE(result->num_distinct_primes(), 3u) << L.value();
  }
  EXPECT_TRUE(any) << "no candidate L produced a Carmichael number";
}

// The balanced chain is injectable through the engine; any result it does
// produce must still be a valid Carmichael number (it may produce none at these
// small sizes, which is acceptable).
TEST(LargeCarmichael, BalancedChainProducesValidOrNone) {
  pseudoprime::BalancedCongruenceChainStrategy balanced;
  for (const auto& L : candidate_Ls()) {
    const auto result = pseudoprime::construct_large_carmichael(L, balanced);
    if (!result) continue;
    EXPECT_TRUE(pseudoprime::is_carmichael(*result)) << L.value();
    EXPECT_EQ(result->value() % L.value(), Integer(1)) << L.value();
  }
}

// Triple (k=3) pairing must preserve correctness: any number it constructs is
// still a genuine Carmichael number with product 1 mod L.
TEST(LargeCarmichael, TriplePairingProducesValidOrNone) {
  for (const auto& L : candidate_Ls()) {
    const auto result = pseudoprime::construct_large_carmichael(
        L, pseudoprime::PrimePowerChainStrategy{}, pseudoprime::LayeredPrimalityTest{},
        /*max_pairing=*/3);
    if (!result) continue;
    EXPECT_TRUE(pseudoprime::is_carmichael(*result)) << L.value();
    EXPECT_EQ(result->value() % L.value(), Integer(1)) << L.value();
  }
}

// Quadruple (k=4) pairing exercises the generalized advance() (groups > 3 for the
// distinguished element). Any number it constructs must still be a genuine
// Carmichael number with product 1 mod L.
TEST(LargeCarmichael, QuadruplePairingProducesValidOrNone) {
  for (const auto& L : candidate_Ls()) {
    const auto result = pseudoprime::construct_large_carmichael(
        L, pseudoprime::PrimePowerChainStrategy{}, pseudoprime::LayeredPrimalityTest{},
        /*max_pairing=*/4);
    if (!result) continue;
    EXPECT_TRUE(pseudoprime::is_carmichael(*result)) << L.value();
    EXPECT_EQ(result->value() % L.value(), Integer(1)) << L.value();
  }
}

// The subgroup-chain strategy is injectable: a user-defined strategy that
// forwards to the prime-power chain yields the same outcome as the default.
TEST(LargeCarmichael, StrategyIsInjectable) {
  struct ForwardingStrategy : pseudoprime::SubgroupChainStrategy {
    std::unique_ptr<pseudoprime::SubgroupChain> build(
        const pseudoprime::Factorization& ghat) const override {
      return pseudoprime::PrimePowerChainStrategy{}.build(ghat);
    }
  };

  for (const auto& L : candidate_Ls()) {
    const auto a = pseudoprime::construct_large_carmichael(L);
    const auto b = pseudoprime::construct_large_carmichael(L, ForwardingStrategy{});
    EXPECT_EQ(a.has_value(), b.has_value()) << L.value();
    if (a && b) EXPECT_EQ(a->value(), b->value()) << L.value();
  }
}

}  // namespace
