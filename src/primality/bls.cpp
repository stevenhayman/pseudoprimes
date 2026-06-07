#include "pseudoprime/primality/bls.hpp"

#include <stdexcept>
#include <vector>

namespace pseudoprime {

namespace {

PrimalityCertificate prime() { return {Primality::Prime, std::nullopt}; }
PrimalityCertificate composite() { return {Primality::Composite, std::nullopt}; }
PrimalityCertificate composite_with(const Integer& factor) {
  return {Primality::Composite, factor};
}
PrimalityCertificate inconclusive() { return {Primality::Inconclusive, std::nullopt}; }

}  // namespace

PrimalityCertificate bls_certifies_prime(const Integer& n, const Factorization& f,
                                         unsigned long max_base) {
  if (n < 2) return composite();
  if (n == 2 || n == 3) return prime();
  if (mpz_even_p(n.get_mpz_t())) return composite();

  const Integer nm1 = n - 1;
  const Integer F = f.value();
  if (F <= 1) {
    throw std::invalid_argument("bls_certifies_prime: F must be > 1");
  }
  if (nm1 % F != 0) {
    throw std::invalid_argument("bls_certifies_prime: F must divide n - 1");
  }
  const Integer R = nm1 / F;
  if (gcd(F, R) != 1) {
    throw std::invalid_argument("bls_certifies_prime: F and (n-1)/F must be coprime");
  }

  // The factored part must exceed n^(1/3) for the test to apply.
  if (F * F * F <= n) return inconclusive();

  // For each prime q | F, search a base a in [2, max_base] that witnesses q.
  const std::vector<Integer> qs = f.primes();
  for (const Integer& q : qs) {
    bool witnessed = false;
    for (unsigned long ai = 2; ai <= max_base; ++ai) {
      const Integer a(ai);

      // Fermat condition. A failure proves compositeness; surface a shared
      // factor when the base exposes one (self-witnessing).
      if (mod_pow(a, nm1, n) != 1) {
        const Integer g = gcd(a, n);
        if (g > 1 && g < n) return composite_with(g);
        return composite();
      }

      const Integer g = gcd(mod_pow(a, nm1 / q, n) - 1, n);
      if (g == 1) {  // a is a valid witness for q
        witnessed = true;
        break;
      }
      if (g > 1 && g < n) return composite_with(g);  // nontrivial gcd reveals a factor
      // g == n: this base says nothing about q; try the next base.
    }
    if (!witnessed) return inconclusive();  // bound exhausted for this prime
  }

  // Pocklington branch: F > sqrt(N) settles primality outright.
  if (F * F > n) return prime();

  // Cube-root branch: N^(1/3) < F <= sqrt(N). Write R = 2*m*F + s and test the
  // discriminant of the only possible nontrivial factorization of N.
  const Integer two_F = 2 * F;
  const Integer m = R / two_F;
  const Integer s = R % two_F;
  const Integer disc = s * s - 8 * m;
  if (m == 0 || disc < 0 || !is_perfect_square(disc)) return prime();
  return composite();
}

}  // namespace pseudoprime
