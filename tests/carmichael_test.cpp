#include "pseudoprime/carmichael.hpp"

#include <gtest/gtest.h>
#include <initializer_list>
#include <vector>

namespace {

using pseudoprime::Factorization;
using pseudoprime::Integer;
using pseudoprime::PrimePower;

Factorization squarefree(std::initializer_list<long> primes) {
  std::vector<Integer> v;
  for (long p : primes) v.push_back(Integer(p));
  return Factorization::from_primes(v);
}

TEST(IsCarmichael, AcceptsKnownCarmichaelNumbers) {
  EXPECT_TRUE(pseudoprime::is_carmichael(squarefree({3, 11, 17})));      // 561
  EXPECT_TRUE(pseudoprime::is_carmichael(squarefree({5, 13, 17})));      // 1105
  EXPECT_TRUE(pseudoprime::is_carmichael(squarefree({7, 11, 13, 41})));  // 41041
}

TEST(IsCarmichael, RejectsPrimeAndTwoPrimeProducts) {
  EXPECT_FALSE(pseudoprime::is_carmichael(squarefree({7})));     // prime: < 2 factors
  EXPECT_FALSE(pseudoprime::is_carmichael(squarefree({3, 5})));  // 15: (5-1) does not divide 14
}

TEST(IsCarmichael, RejectsEvenComposite) {
  EXPECT_FALSE(pseudoprime::is_carmichael(squarefree({2, 3})));  // 6: (3-1) does not divide 5
}

TEST(IsCarmichael, RejectsNonSquarefree) {
  // 9 = 3^2: prime power, not squarefree.
  EXPECT_FALSE(pseudoprime::is_carmichael(Factorization::from_prime_power(Integer(3), 2)));
  // 3^2 * 11 * 17: has the right primes but is not squarefree.
  const Factorization f({{Integer(3), 2}, {Integer(11), 1}, {Integer(17), 1}});
  EXPECT_FALSE(pseudoprime::is_carmichael(f));
}

// At order 1 the general and rigid predicates coincide with ordinary Korselt.
TEST(IsCarmichael, RigidEqualsGeneralAtOrderOne) {
  for (const Factorization& n :
       {squarefree({3, 11, 17}), squarefree({5, 13, 17}), squarefree({7, 11, 13, 41})}) {
    EXPECT_TRUE(pseudoprime::is_carmichael(n, 1));
    EXPECT_TRUE(pseudoprime::is_rigid_carmichael(n, 1));
    EXPECT_TRUE(pseudoprime::is_rigid_carmichael(n));  // default order == 1
  }
}

// Howe's example (Sec. 1) of a rigid Carmichael number of order 2.
TEST(IsCarmichael, AcceptsKnownRigidOrderTwo) {
  const Factorization n = squarefree({17, 31, 41, 43, 89, 97, 167, 331});
  EXPECT_TRUE(pseudoprime::is_carmichael(n, 1));        // also ordinary Carmichael
  EXPECT_TRUE(pseudoprime::is_carmichael(n, 2));        // order 2, general
  EXPECT_TRUE(pseudoprime::is_rigid_carmichael(n, 2));  // and rigid
}

// An ordinary Carmichael number that is not a Carmichael number of order 2.
TEST(IsCarmichael, RejectsOrdinaryCarmichaelAtOrderTwo) {
  const Factorization n = squarefree({3, 11, 17});  // 561
  EXPECT_TRUE(pseudoprime::is_carmichael(n, 1));
  EXPECT_FALSE(pseudoprime::is_carmichael(n, 2));
  EXPECT_FALSE(pseudoprime::is_rigid_carmichael(n, 2));
}

// Howe's smallest non-rigid Carmichael number of order 2 (Sec. 5, built from the
// distinguished prime p_0 = 1153): it satisfies the general order-2 criterion but
// is NOT rigid -- exactly the general-vs-rigid distinction.
TEST(IsCarmichael, DistinguishesNonRigidOrderTwo) {
  const Factorization n = squarefree({23, 67, 71, 89, 109, 113, 191, 199, 233, 239, 271, 307, 373,
                                      419, 521, 911, 929, 1153, 1217, 1429, 2089, 2729, 23561});
  EXPECT_TRUE(pseudoprime::is_carmichael(n, 2));         // general: passes
  EXPECT_FALSE(pseudoprime::is_rigid_carmichael(n, 2));  // rigid: fails (p_0 = 1153)
}

}  // namespace
