#include "pseudoprime/primality/lucas.hpp"

#include <vector>

namespace pseudoprime {

Primality lucas_certifies_prime(const Integer& p, const Factorization& p_minus_one,
                                unsigned long max_base) {
  if (p < 2) return Primality::Inconclusive;
  const Integer pm1 = p - 1;
  const std::vector<Integer> qs = p_minus_one.primes();

  for (unsigned long ai = 2; ai <= max_base; ++ai) {
    const Integer a(ai);
    if (mod_pow(a, pm1, p) != 1) continue;  // not even a Fermat base; cannot be a primitive root

    bool primitive = true;
    for (const Integer& q : qs) {
      if (gcd(mod_pow(a, pm1 / q, p) - 1, p) != 1) {  // a has order < p-1 modulo some prime factor
        primitive = false;
        break;
      }
    }
    if (primitive) return Primality::Prime;  // a is a primitive root: proof that p is prime
  }
  return Primality::Inconclusive;  // no witness within the bound (NOT a proof of compositeness)
}

}  // namespace pseudoprime
