#include "pseudoprime/factorization.hpp"

#include <gtest/gtest.h>

namespace {

// L = 80 = 2^4 * 5
pseudoprime::Factorization make_80() {
  return pseudoprime::Factorization({{pseudoprime::Integer(2), 4}, {pseudoprime::Integer(5), 1}});
}

TEST(Factorization, Value) {
  EXPECT_EQ(make_80().value(), pseudoprime::Integer(80));
  EXPECT_EQ(pseudoprime::Factorization{}.value(), pseudoprime::Integer(1));
}

TEST(Factorization, Squarefree) {
  EXPECT_FALSE(make_80().is_squarefree());
  EXPECT_TRUE(pseudoprime::Factorization::from_primes({pseudoprime::Integer(3), pseudoprime::Integer(11),
                                              pseudoprime::Integer(17)})
                  .is_squarefree());
}

TEST(Factorization, DivisorsOf80) {
  const auto divs = make_80().divisors();
  const std::vector<pseudoprime::Integer> expected = {1, 2, 4, 5, 8, 10, 16, 20, 40, 80};
  ASSERT_EQ(divs.size(), expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) EXPECT_EQ(divs[i], expected[i]);
}

TEST(Factorization, DivisorsFactoredReconstruct) {
  const auto facs = make_80().divisors_factored();
  EXPECT_EQ(facs.size(), 10u);  // (4+1)*(1+1)
  for (const auto& f : facs) {
    // Each factored divisor reconstructs to a genuine divisor of 80.
    EXPECT_EQ(pseudoprime::Integer(80) % f.value(), pseudoprime::Integer(0));
  }
}

TEST(Factorization, EulerPhiAndCarmichaelLambda) {
  EXPECT_EQ(make_80().euler_phi(), pseudoprime::Integer(32));        // phi(80)
  EXPECT_EQ(make_80().carmichael_lambda(), pseudoprime::Integer(4));  // lambda(80)
}

TEST(Factorization, FactoredGcd) {
  // gcd(2^4 * 5, 2^2 * 3) = 2^2 = 4, by min-exponent intersection.
  const pseudoprime::Factorization a({{pseudoprime::Integer(2), 4}, {pseudoprime::Integer(5), 1}});
  const pseudoprime::Factorization b({{pseudoprime::Integer(2), 2}, {pseudoprime::Integer(3), 1}});
  const pseudoprime::Factorization g = pseudoprime::gcd(a, b);
  EXPECT_EQ(g.value(), pseudoprime::Integer(4));

  // Cross-check against the ordinary integer gcd.
  pseudoprime::Integer expected;
  mpz_gcd(expected.get_mpz_t(), a.value().get_mpz_t(), b.value().get_mpz_t());
  EXPECT_EQ(g.value(), expected);
}

TEST(Factorization, FromPrimePower) {
  EXPECT_EQ(pseudoprime::Factorization::from_prime_power(pseudoprime::Integer(7), 3).value(),
            pseudoprime::Integer(343));
  // Exponent 0 collapses to the empty product (== 1).
  EXPECT_TRUE(pseudoprime::Factorization::from_prime_power(pseudoprime::Integer(7), 0).is_one());
}

TEST(Factorization, FromPrimesMergesDuplicates) {
  const auto f = pseudoprime::Factorization::from_primes(
      {pseudoprime::Integer(2), pseudoprime::Integer(3), pseudoprime::Integer(2)});
  EXPECT_EQ(f.value(), pseudoprime::Integer(12));  // 2^2 * 3
  EXPECT_EQ(f.multiplicity(pseudoprime::Integer(2)), 2u);
  EXPECT_FALSE(f.is_squarefree());
}

TEST(Factorization, PrimesAndMultiplicity) {
  const pseudoprime::Factorization f({{pseudoprime::Integer(2), 4}, {pseudoprime::Integer(5), 1}});
  EXPECT_EQ(f.primes(), (std::vector<pseudoprime::Integer>{2, 5}));
  EXPECT_EQ(f.multiplicity(pseudoprime::Integer(2)), 4u);
  EXPECT_EQ(f.multiplicity(pseudoprime::Integer(5)), 1u);
  EXPECT_EQ(f.multiplicity(pseudoprime::Integer(3)), 0u);  // absent prime
  EXPECT_EQ(f.num_distinct_primes(), 2u);
  EXPECT_FALSE(f.is_one());
}

TEST(Factorization, CarmichaelLambdaPowersOfTwo) {
  // lambda(2)=1, lambda(4)=2, lambda(8)=2, lambda(16)=4, lambda(32)=8.
  auto lambda2 = [](unsigned long e) {
    return pseudoprime::Factorization::from_prime_power(pseudoprime::Integer(2), e)
        .carmichael_lambda();
  };
  EXPECT_EQ(lambda2(1), pseudoprime::Integer(1));
  EXPECT_EQ(lambda2(2), pseudoprime::Integer(2));
  EXPECT_EQ(lambda2(3), pseudoprime::Integer(2));  // the 2^e (e>=3) branch
  EXPECT_EQ(lambda2(4), pseudoprime::Integer(4));
  EXPECT_EQ(lambda2(5), pseudoprime::Integer(8));
}

TEST(Factorization, EulerPhiOddPrimePower) {
  // phi(27) = 27 - 9 = 18.
  EXPECT_EQ(pseudoprime::Factorization::from_prime_power(pseudoprime::Integer(3), 3).euler_phi(),
            pseudoprime::Integer(18));
}

TEST(Factorization, PrimePowers) {
  // 80 = 2^4 * 5 -> components [16, 5], ascending by prime.
  EXPECT_EQ(make_80().prime_powers(), (std::vector<pseudoprime::Integer>{16, 5}));
}

TEST(Factorization, GcdOfCoprimeIsOne) {
  const pseudoprime::Factorization a({{pseudoprime::Integer(3), 1}});
  const pseudoprime::Factorization b({{pseudoprime::Integer(5), 1}});
  EXPECT_TRUE(pseudoprime::gcd(a, b).is_one());
}

}  // namespace
