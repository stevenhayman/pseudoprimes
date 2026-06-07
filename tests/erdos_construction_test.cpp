#include "pseudoprime/carmichael.hpp"
#include "pseudoprime/erdos_construction.hpp"
#include "pseudoprime/subset_product.hpp"

#include <algorithm>
#include <gtest/gtest.h>
#include <vector>

namespace {

std::vector<pseudoprime::Integer> prime_set_values(const pseudoprime::Factorization& L) {
  std::vector<pseudoprime::Integer> out;
  for (const auto& p : pseudoprime::generate_prime_set(L)) out.push_back(p);
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<pseudoprime::Integer> construct_all(const pseudoprime::Factorization& L) {
  pseudoprime::BruteForceSolver solver;
  std::vector<pseudoprime::Integer> out;
  for (const pseudoprime::ConstructedCarmichael& cc : pseudoprime::construct_carmichael(L, solver)) {
    EXPECT_TRUE(pseudoprime::is_carmichael(cc.factorization)) << cc.value;
    EXPECT_EQ(cc.factorization.value(), cc.value);
    out.push_back(cc.value);
  }
  std::sort(out.begin(), out.end());
  return out;
}

TEST(Erdos, PrimeSetForL80) {
  const pseudoprime::Factorization L({{pseudoprime::Integer(2), 4}, {pseudoprime::Integer(5), 1}});  // 80
  EXPECT_EQ(prime_set_values(L),
            (std::vector<pseudoprime::Integer>{3, 11, 17, 41}));
}

TEST(Erdos, ConstructsSmallestCarmichael) {
  const pseudoprime::Factorization L({{pseudoprime::Integer(2), 4}, {pseudoprime::Integer(5), 1}});  // 80
  EXPECT_EQ(construct_all(L), (std::vector<pseudoprime::Integer>{561}));
}

TEST(Erdos, PrimeSetForL48) {
  const pseudoprime::Factorization L({{pseudoprime::Integer(2), 4}, {pseudoprime::Integer(3), 1}});  // 48
  EXPECT_EQ(prime_set_values(L),
            (std::vector<pseudoprime::Integer>{5, 7, 13, 17}));
}

TEST(Erdos, ConstructsCarmichael1105) {
  const pseudoprime::Factorization L({{pseudoprime::Integer(2), 4}, {pseudoprime::Integer(3), 1}});  // 48
  EXPECT_EQ(construct_all(L), (std::vector<pseudoprime::Integer>{1105}));
}

TEST(Erdos, FindFirstStopsEarly) {
  const pseudoprime::Factorization L({{pseudoprime::Integer(2), 4}, {pseudoprime::Integer(5), 1}});
  pseudoprime::BruteForceSolver solver;
  int count = 0;
  for ([[maybe_unused]] const auto& cc : pseudoprime::construct_carmichael(L, solver)) {
    ++count;
    break;  // stop after first
  }
  EXPECT_EQ(count, 1);
}

// The order-2 prime set P(2,L) keeps only primes p with (p^2 - 1) | L, so it is a
// subset of the order-1 set P(L). For L = 48: P(1) = {5,7,13,17}; of these only
// 5 (24|48) and 7 (48|48) survive (13^2-1 = 168, 17^2-1 = 288 do not divide 48).
TEST(Erdos, OrderTwoPrimeSetFiltersByPSquaredMinusOne) {
  const pseudoprime::Factorization L({{pseudoprime::Integer(2), 4}, {pseudoprime::Integer(3), 1}});
  const pseudoprime::Integer Lval = L.value();  // 48

  std::vector<pseudoprime::Integer> order2;
  pseudoprime::LayeredPrimalityTest primality;
  for (const auto& p : pseudoprime::generate_prime_set(L, primality, 2)) order2.push_back(p);
  std::sort(order2.begin(), order2.end());
  EXPECT_EQ(order2, (std::vector<pseudoprime::Integer>{5, 7}));

  // Subset of P(1,L), and each survivor genuinely satisfies (p^2 - 1) | L.
  const std::vector<pseudoprime::Integer> order1 = prime_set_values(L);
  for (const auto& p : order2) {
    EXPECT_NE(std::find(order1.begin(), order1.end(), p), order1.end());
    EXPECT_EQ(Lval % (p * p - 1), 0);
  }
}

// Anything the order-2 construction emits must be a rigid Carmichael number of
// order 2. (Brute-force-feasible small L do not yield an order-2 number -- the
// paper needs richly composite L with |P(2,L)| ~ 45+ -- so this checks the
// invariant on the emitted stream, which is empty here, without crashing.)
TEST(Erdos, OrderTwoConstructionYieldsRigidOrderTwo) {
  const pseudoprime::Factorization L({{pseudoprime::Integer(2), 4}, {pseudoprime::Integer(3), 1}});
  pseudoprime::BruteForceSolver solver;
  std::size_t count = 0;
  for (const pseudoprime::ConstructedCarmichael& cc :
       pseudoprime::construct_carmichael(L, solver, /*min_factors=*/3, /*order=*/2)) {
    EXPECT_TRUE(pseudoprime::is_rigid_carmichael(cc.factorization, 2)) << cc.value;
    ++count;
  }
  EXPECT_EQ(count, 0u);  // P(2,48) = {5,7} cannot form a >=3-prime subset product
}

}  // namespace
