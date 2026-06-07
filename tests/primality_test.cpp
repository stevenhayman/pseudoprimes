#include "pseudoprime/order.hpp"
#include "pseudoprime/primality.hpp"

#include <gtest/gtest.h>
#include <stdexcept>

namespace {

// The primality tests require the factorization of p - 1.
pseudoprime::Factorization fac_of(unsigned long n) {
  // Small helper: trial-factor n (only used to build test inputs).
  std::vector<pseudoprime::Integer> primes;
  for (unsigned long d = 2; d * d <= n; ++d) {
    while (n % d == 0) {
      primes.push_back(pseudoprime::Integer(d));
      n /= d;
    }
  }
  if (n > 1) primes.push_back(pseudoprime::Integer(n));
  return pseudoprime::Factorization::from_primes(primes);
}

TEST(FellowsKoblitz, CertifiesPrimes) {
  EXPECT_TRUE(pseudoprime::is_prime_fk(pseudoprime::Integer(17), fac_of(16)));
  EXPECT_TRUE(pseudoprime::is_prime_fk(pseudoprime::Integer(41), fac_of(40)));
  EXPECT_TRUE(pseudoprime::is_prime_fk(pseudoprime::Integer(101), fac_of(100)));
  EXPECT_TRUE(pseudoprime::is_prime_fk(pseudoprime::Integer(7919), fac_of(7918)));
  EXPECT_TRUE(pseudoprime::is_prime_fk(pseudoprime::Integer(2147483647),
                                       fac_of(2147483646)));  // 2^31 - 1
}

TEST(FellowsKoblitz, RejectsComposites) {
  EXPECT_FALSE(pseudoprime::is_prime_fk(pseudoprime::Integer(21), fac_of(20)));
  EXPECT_FALSE(pseudoprime::is_prime_fk(pseudoprime::Integer(25), fac_of(24)));
  // Carmichael numbers pass every Fermat test for coprime bases; the test must
  // still reject them (via the gcd witness or the lcm threshold).
  EXPECT_FALSE(pseudoprime::is_prime_fk(pseudoprime::Integer(561), fac_of(560)));
  EXPECT_FALSE(pseudoprime::is_prime_fk(pseudoprime::Integer(1105), fac_of(1104)));
  EXPECT_FALSE(pseudoprime::is_prime_fk(pseudoprime::Integer(41041), fac_of(41040)));
}

TEST(FellowsKoblitz, SelfWitnessedFactorDividesN) {
  // When a factor is returned on a composite verdict, it must be nontrivial.
  for (unsigned long n : {15u, 21u, 91u, 561u, 1105u}) {
    const auto cert = pseudoprime::is_prime_fellows_koblitz(pseudoprime::Integer(n),
                                                            fac_of(n - 1));
    EXPECT_EQ(cert.verdict, pseudoprime::Primality::Composite) << n;
    if (cert.factor) {
      EXPECT_GT(*cert.factor, pseudoprime::Integer(1)) << n;
      EXPECT_LT(*cert.factor, pseudoprime::Integer(n)) << n;
      EXPECT_EQ(pseudoprime::Integer(n) % *cert.factor, pseudoprime::Integer(0)) << n;
    }
  }
}

TEST(FellowsKoblitz, HandlesSmallAndEdgeN) {
  EXPECT_TRUE(pseudoprime::is_prime_fk(pseudoprime::Integer(2), fac_of(1)));
  EXPECT_TRUE(pseudoprime::is_prime_fk(pseudoprime::Integer(3), fac_of(2)));
  EXPECT_TRUE(pseudoprime::is_prime_fk(pseudoprime::Integer(5), fac_of(4)));
  EXPECT_TRUE(pseudoprime::is_prime_fk(pseudoprime::Integer(7), fac_of(6)));
  EXPECT_FALSE(pseudoprime::is_prime_fk(pseudoprime::Integer(9), fac_of(8)));   // 3^2
  EXPECT_FALSE(pseudoprime::is_prime_fk(pseudoprime::Integer(15), fac_of(14)));
}

TEST(FellowsKoblitz, ThrowsOnWrongFactorization) {
  // 10 is not the factorization of 17 - 1 = 16 -> misuse must throw.
  EXPECT_THROW(pseudoprime::is_prime_fellows_koblitz(pseudoprime::Integer(17), fac_of(10)),
               std::invalid_argument);
}

TEST(Lucas, CertifiesPrimes) {
  EXPECT_EQ(pseudoprime::lucas_certifies_prime(pseudoprime::Integer(17), fac_of(16)),
            pseudoprime::Primality::Prime);
  EXPECT_EQ(pseudoprime::lucas_certifies_prime(pseudoprime::Integer(41), fac_of(40)),
            pseudoprime::Primality::Prime);
  EXPECT_EQ(
      pseudoprime::lucas_certifies_prime(pseudoprime::Integer(2147483647), fac_of(2147483646)),
      pseudoprime::Primality::Prime);
}

TEST(Lucas, InconclusiveForComposite) {
  // No witness exists for a composite, so the result is Inconclusive -- NOT a
  // proof of compositeness; the layered test resolves these via Fermat/FK.
  EXPECT_EQ(pseudoprime::lucas_certifies_prime(pseudoprime::Integer(21), fac_of(20)),
            pseudoprime::Primality::Inconclusive);
  EXPECT_EQ(pseudoprime::lucas_certifies_prime(pseudoprime::Integer(561), fac_of(560)),
            pseudoprime::Primality::Inconclusive);
}

TEST(LayeredPrimalityTest, MatchesFellowsKoblitz) {
  const pseudoprime::LayeredPrimalityTest layered;
  const pseudoprime::FellowsKoblitzPrimalityTest fk;
  for (unsigned long n : {2u, 3u, 5u, 7u, 9u, 15u, 17u, 21u, 25u, 41u, 101u, 561u, 1105u, 7919u,
                          41041u}) {
    const pseudoprime::Factorization nm1 = fac_of(n - 1);
    EXPECT_EQ(layered.is_prime(pseudoprime::Integer(n), nm1),
              fk.is_prime(pseudoprime::Integer(n), nm1))
        << n;
  }
  EXPECT_TRUE(layered.is_prime(pseudoprime::Integer(2147483647), fac_of(2147483646)));
}

// Build a fully-factored partial part F from explicit prime factors (with
// multiplicity). Used to feed BLS only PART of the factorization of n - 1.
pseudoprime::Factorization partial(std::initializer_list<unsigned long> primes) {
  std::vector<pseudoprime::Integer> v;
  for (unsigned long p : primes) v.push_back(pseudoprime::Integer(p));
  return pseudoprime::Factorization::from_primes(v);
}

TEST(BLS, CertifiesPrimesFromPartialFactorization) {
  // 1009 prime; F = 2^4 = 16 is only PART of 1008 = 2^4*3^2*7. R = 63, gcd = 1,
  // and 16 < sqrt(1009) so this exercises the cube-root + discriminant branch.
  EXPECT_EQ(pseudoprime::bls_certifies_prime(pseudoprime::Integer(1009), partial({2, 2, 2, 2}))
                .verdict,
            pseudoprime::Primality::Prime);

  // 2^31 - 1 prime; 2^31 - 2 = 2*3^2*7*11*31*151*331. Take F = 11*331 = 3641,
  // coprime to the cofactor, with N^(1/3) (~1290) < F < sqrt(N) (~46341):
  // cube-root branch on a large prime.
  EXPECT_EQ(
      pseudoprime::bls_certifies_prime(pseudoprime::Integer(2147483647), partial({11, 331}))
          .verdict,
      pseudoprime::Primality::Prime);
}

TEST(BLS, CertifiesPrimesPocklingtonBranch) {
  // F > sqrt(N): settled without the discriminant.
  EXPECT_EQ(pseudoprime::bls_certifies_prime(pseudoprime::Integer(17), partial({2, 2, 2, 2}))
                .verdict,
            pseudoprime::Primality::Prime);  // F = 16 = 17 - 1
  EXPECT_EQ(
      pseudoprime::bls_certifies_prime(pseudoprime::Integer(41), partial({2, 2, 2})).verdict,
      pseudoprime::Primality::Prime);  // F = 8 > sqrt(41), R = 5 coprime
}

TEST(BLS, NeverCertifiesComposites) {
  // Composites must never be certified prime (Inconclusive or Composite is fine).
  EXPECT_NE(
      pseudoprime::bls_certifies_prime(pseudoprime::Integer(561), partial({2, 2, 2, 2})).verdict,
      pseudoprime::Primality::Prime);  // Carmichael, F = 16 | 560
  EXPECT_NE(pseudoprime::bls_certifies_prime(pseudoprime::Integer(1105), partial({2, 2, 2, 2}))
                .verdict,
            pseudoprime::Primality::Prime);  // Carmichael, F = 16 | 1104
}

TEST(BLS, SelfWitnessedFactorDividesN) {
  // 561 = 3*11*17 (Carmichael), F = 16 | 560. Base 2 passes Fermat but yields
  // gcd == n; base 3 shares a factor with n, so the search surfaces it.
  const auto cert =
      pseudoprime::bls_certifies_prime(pseudoprime::Integer(561), partial({2, 2, 2, 2}));
  EXPECT_EQ(cert.verdict, pseudoprime::Primality::Composite);
  ASSERT_TRUE(cert.factor.has_value());
  EXPECT_GT(*cert.factor, pseudoprime::Integer(1));
  EXPECT_LT(*cert.factor, pseudoprime::Integer(561));
  EXPECT_EQ(pseudoprime::Integer(561) % *cert.factor, pseudoprime::Integer(0));
}

TEST(BLS, InconclusiveIsExplicit) {
  // True prime, but the witness search bound is too small to find one.
  EXPECT_EQ(
      pseudoprime::bls_certifies_prime(pseudoprime::Integer(31), partial({2, 3}), 1).verdict,
      pseudoprime::Primality::Inconclusive);

  // True prime, but F = 7 is too small (F^3 = 343 < 1009) to certify anything.
  EXPECT_EQ(pseudoprime::bls_certifies_prime(pseudoprime::Integer(1009), partial({7})).verdict,
            pseudoprime::Primality::Inconclusive);
}

TEST(BLS, ThrowsOnInvalidF) {
  // F must divide n - 1: 4 does not divide 30.
  EXPECT_THROW(pseudoprime::bls_certifies_prime(pseudoprime::Integer(31), partial({2, 2})),
               std::invalid_argument);
  // F and (n-1)/F must be coprime: F = 2 shares a factor with 1008/2 = 504.
  EXPECT_THROW(pseudoprime::bls_certifies_prime(pseudoprime::Integer(1009), partial({2})),
               std::invalid_argument);
}

TEST(Order, MultiplicativeOrderModPrime) {
  // (Z/7)* has order 6 = 2 * 3.
  const pseudoprime::Factorization six = fac_of(6);
  EXPECT_EQ(pseudoprime::multiplicative_order(pseudoprime::Integer(2), pseudoprime::Integer(7), six),
            pseudoprime::Integer(3));
  EXPECT_EQ(pseudoprime::multiplicative_order(pseudoprime::Integer(3), pseudoprime::Integer(7), six),
            pseudoprime::Integer(6));
  EXPECT_EQ(pseudoprime::multiplicative_order(pseudoprime::Integer(1), pseudoprime::Integer(7), six),
            pseudoprime::Integer(1));
}

TEST(Order, MultiplicativeOrderCompositeModulus) {
  // ord_15(2): 2, 4, 8, 16 == 1 -> order 4. Group order phi(15) = 8 = 2^3.
  EXPECT_EQ(
      pseudoprime::multiplicative_order(pseudoprime::Integer(2), pseudoprime::Integer(15), fac_of(8)),
      pseudoprime::Integer(4));
}

}  // namespace
