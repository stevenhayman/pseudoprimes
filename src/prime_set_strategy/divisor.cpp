#include "pseudoprime/prime_set_strategy/divisor.hpp"

#include "detail.hpp"

namespace pseudoprime {

namespace detail {

bool satisfies_higher_order(const Integer& p, const Integer& L_value, unsigned long order) {
  Integer power = p;  // p^1
  for (unsigned long r = 2; r <= order; ++r) {
    power *= p;  // p^r
    // p^r - 1 can exceed L; then it cannot divide L and the modulo rejects it.
    if (L_value % (power - 1) != 0) return false;
  }
  return true;
}

}  // namespace detail

std::generator<Integer> DivisorPrimeSetStrategy::generate(
    const Factorization& L, const PrimalityTest& primality, unsigned long order) const {
  const Integer L_value = L.value();

  // Stream the divisors of L (each carrying its own factorization) and test
  // p = divisor + 1, so we never hold all divisors at once -- only the kept primes.
  for (const Factorization& divisor_factorization : L.divisors_lazy()) {
    const Integer divisor = divisor_factorization.value();
    const Integer p = divisor + 1;

    // Carmichael factors are odd primes, so p must be odd (i.e. divisor even).
    if (mpz_even_p(p.get_mpz_t())) continue;
    // Require p does not divide L (else the product can never be 1 mod L).
    if (L_value % p == 0) continue;
    // Certify primality deterministically using the free factorization of
    // p - 1 = divisor, via the configured (deterministic) primality test.
    if (!primality.is_prime(p, divisor_factorization)) continue;
    // Order-m membership: also require (p^r - 1) | L for r = 2..order.
    if (order > 1 && !detail::satisfies_higher_order(p, L_value, order)) continue;

    co_yield p;
  }
}

}  // namespace pseudoprime
