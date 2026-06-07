#include "pseudoprime/prime_set_strategy.hpp"

#include <algorithm>
#include <cstdint>
#include <gtest/gtest.h>
#include <vector>

#include "pseudoprime/factorization.hpp"
#include "pseudoprime/integer.hpp"
#include "pseudoprime/primality.hpp"

namespace {

using pseudoprime::DivisorPrimeSetStrategy;
using pseudoprime::Factorization;
using pseudoprime::Integer;
using pseudoprime::LayeredPrimalityTest;
using pseudoprime::PrimeSetStrategy;
using pseudoprime::SievePrimeSetStrategy;

// Collect a strategy's prime set, canonicalized by ascending prime.
std::vector<Integer> collect(const PrimeSetStrategy& strategy, const Factorization& L,
                             unsigned long order) {
  const LayeredPrimalityTest primality;
  std::vector<Integer> out;
  for (const Integer& p : strategy.generate(L, primality, order)) out.push_back(p);
  std::sort(out.begin(), out.end());
  return out;
}

// Representative moduli: small hand-checked, plus a few highly composite L
// (still 64-bit so the cap == 0 "complete" sieve covers every divisor).
std::vector<Factorization> test_moduli() {
  return {
      Factorization({{Integer(2), 4}, {Integer(3), 1}}),                  // 48
      Factorization({{Integer(2), 4}, {Integer(5), 1}}),                  // 80
      Factorization({{Integer(2), 3}, {Integer(3), 2}, {Integer(5), 1}, {Integer(7), 1}}),  // 2520
      Factorization({{Integer(2), 4}, {Integer(3), 2}, {Integer(5), 1}, {Integer(7), 1}}),  // 5040
      Factorization({{Integer(2), 3}, {Integer(3), 2}, {Integer(5), 1}, {Integer(7), 1},
                     {Integer(11), 1}}),                                  // 27720
      Factorization({{Integer(2), 4}, {Integer(3), 2}, {Integer(5), 1}, {Integer(7), 1},
                     {Integer(11), 1}}),                                  // 55440
  };
}

// With cap == 0 (complete: sweep every divisor d <= L) the sieve strategy must
// produce a result byte-identical to the divisor strategy, for every tested L
// and order 1..3. This is the 100%-certainty guarantee, and also confirms the
// generate_prime_set refactor onto DivisorPrimeSetStrategy is behaviour-preserving.
TEST(SievePrimeSet, MatchesDivisorStrategyComplete) {
  const DivisorPrimeSetStrategy divisor;
  const SievePrimeSetStrategy sieve(/*cap=*/0, /*window=*/256);  // tiny window to exercise segmenting

  for (const Factorization& L : test_moduli()) {
    for (unsigned long order = 1; order <= 3; ++order) {
      const auto expected = collect(divisor, L, order);
      const auto actual = collect(sieve, L, order);
      EXPECT_EQ(expected, actual) << "L = " << L.value() << ", order = " << order;
    }
  }
}

// Known small prime sets, computed via the sieve directly.
TEST(SievePrimeSet, KnownSmallPrimeSets) {
  const SievePrimeSetStrategy sieve;  // default cap (complete) and window

  const Factorization L80({{Integer(2), 4}, {Integer(5), 1}});  // 80
  EXPECT_EQ(collect(sieve, L80, 1), (std::vector<Integer>{3, 11, 17, 41}));

  const Factorization L48({{Integer(2), 4}, {Integer(3), 1}});  // 48
  EXPECT_EQ(collect(sieve, L48, 1), (std::vector<Integer>{5, 7, 13, 17}));
  EXPECT_EQ(collect(sieve, L48, 2), (std::vector<Integer>{5, 7}));
}

// A bounded cap returns exactly the members of P(m,L) with p <= cap + 1: the
// divisor-strategy result filtered to that bound, nothing spurious or missing.
TEST(SievePrimeSet, BoundedCapIsExactSubset) {
  const Factorization L({{Integer(2), 3}, {Integer(3), 2}, {Integer(5), 1}, {Integer(7), 1},
                         {Integer(11), 1}});  // 27720
  const DivisorPrimeSetStrategy divisor;

  for (std::uint64_t cap : {std::uint64_t(50), std::uint64_t(200), std::uint64_t(1000)}) {
    const SievePrimeSetStrategy sieve(cap, /*window=*/64);
    const auto full = collect(divisor, L, 1);

    std::vector<Integer> expected;
    for (const Integer& p : full)
      if (p <= Integer(static_cast<unsigned long>(cap)) + 1) expected.push_back(p);

    EXPECT_EQ(collect(sieve, L, 1), expected) << "cap = " << cap;
  }
}

// provable_cap(L, m) = ceil((L+1)^(1/m)): the smallest cap that misses nothing,
// derived purely from L and m. Check the defining inequality r^m >= L+1 with
// (r-1)^m < L+1 (r is the exact ceiling), for several L and orders.
TEST(SievePrimeSet, ProvableCapIsExactCeiling) {
  using pseudoprime::provable_cap;
  const std::vector<Factorization> Ls = test_moduli();
  for (const Factorization& L : Ls) {
    const Integer base = L.value() + 1;
    for (unsigned long m = 1; m <= 6; ++m) {
      const Integer r = provable_cap(L, m);
      Integer rm;
      mpz_pow_ui(rm.get_mpz_t(), r.get_mpz_t(), m);
      EXPECT_GE(rm, base) << "L=" << L.value() << " m=" << m;  // r^m >= L+1
      if (r > 1) {
        Integer rm1;
        const Integer rprev = r - 1;
        mpz_pow_ui(rm1.get_mpz_t(), rprev.get_mpz_t(), m);
        EXPECT_LT(rm1, base) << "L=" << L.value() << " m=" << m;  // (r-1)^m < L+1
      }
    }
    EXPECT_EQ(provable_cap(L, 1), L.value() + 1);  // order 1: exactly L+1
  }
}

// With the default cap == 0 (now the per-order ceiling), the sieve is still
// byte-identical to the full divisor strategy for every order: no order-m member
// can have d beyond ceil((L+1)^(1/m)).
TEST(SievePrimeSet, DefaultCapMatchesDivisorPerOrder) {
  const DivisorPrimeSetStrategy divisor;
  const SievePrimeSetStrategy sieve;  // cap == 0 -> provable_cap per order
  for (const Factorization& L : test_moduli()) {
    for (unsigned long order = 1; order <= 3; ++order) {
      EXPECT_EQ(collect(sieve, L, order), collect(divisor, L, order))
          << "L = " << L.value() << ", order = " << order;
    }
  }
}

// The sieve yields in ascending prime order without an external sort.
TEST(SievePrimeSet, EmitsAscending) {
  const Factorization L({{Integer(2), 4}, {Integer(3), 2}, {Integer(5), 1}, {Integer(7), 1}});  // 5040
  const SievePrimeSetStrategy sieve(/*cap=*/0, /*window=*/128);
  const LayeredPrimalityTest primality;

  Integer prev = 0;
  for (const Integer& p : sieve.generate(L, primality, 1)) {
    EXPECT_LT(prev, p);
    prev = p;
  }
}

}  // namespace
