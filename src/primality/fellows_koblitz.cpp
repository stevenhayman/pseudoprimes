#include "pseudoprime/primality/fellows_koblitz.hpp"

#include <stdexcept>

#include "pseudoprime/order.hpp"

namespace pseudoprime {

namespace {

PrimalityCertificate prime() { return {Primality::Prime, std::nullopt}; }
PrimalityCertificate composite() { return {Primality::Composite, std::nullopt}; }
PrimalityCertificate composite_with(const Integer& factor) {
  return {Primality::Composite, factor};
}

}  // namespace

PrimalityCertificate is_prime_fellows_koblitz(const Integer& n,
                                              const Factorization& n_minus_one) {
  if (n < 2) return composite();
  if (n == 2 || n == 3) return prime();
  if (mpz_even_p(n.get_mpz_t())) return composite();

  const Integer nm1 = n - 1;
  if (n_minus_one.value() != nm1) {
    throw std::invalid_argument(
        "is_prime_fellows_koblitz: n_minus_one must be the factorization of n-1");
  }

  // Base bound B = (log2 n)^2, capped at n-1 so bases stay in [2, n-1] for tiny n.
  unsigned long bits = mpz_sizeinbase(n.get_mpz_t(), 2);
  unsigned long B = bits * bits;
  if (nm1 < Integer(B)) B = mpz_get_ui(nm1.get_mpz_t());

  const std::vector<Integer> qs = n_minus_one.primes();
  Integer h = 1;  // lcm of the element orders seen so far

  for (unsigned long ai = 2; ai <= B; ++ai) {
    const Integer a(ai);

    // (2) Fermat test. A failure proves compositeness; if a shares a factor
    //     with n, surface it (self-witnessing).
    if (mod_pow(a, nm1, n) != 1) {
      const Integer g = gcd(a, n);
      if (g > 1 && g < n) return composite_with(g);
      return composite();
    }

    // (3) exact order of a modulo n, from the known factorization of n-1.
    const Integer ord = multiplicative_order(a, n, n_minus_one);

    // (4)-(5) for each prime q | ord, a nontrivial gcd reveals a factor of n.
    for (const Integer& q : qs) {
      if (ord % q != 0) continue;
      const Integer g = gcd(mod_pow(a, ord / q, n) - 1, n);
      if (g > 1 && g < n) return composite_with(g);
    }

    // (6) accumulate h = lcm(h, ord).
    h = lcm(h, ord);
  }

  // (7) h >= sqrt(n)  <=>  h*h >= n. If so every prime factor exceeds sqrt(n),
  //     so n is prime; otherwise the smooth-root argument makes n composite.
  if (h * h < n) return composite();
  return prime();
}

}  // namespace pseudoprime
