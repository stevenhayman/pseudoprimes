#include "pseudoprime/carmichael.hpp"

namespace pseudoprime {

namespace {

// Is n congruent to some power p^i (i >= 0) modulo p^r - 1?
//
// Since gcd(p, p^r - 1) = 1 and p has multiplicative order exactly r modulo
// p^r - 1 (p^r == 1 there, and p^j - 1 < p^r - 1 for j < r), the reachable
// residues { p^i mod (p^r - 1) } are exactly the r elements
// { p^0, p^1, ..., p^{r-1} } mod (p^r - 1). So this is a finite membership test
// -- no discrete logarithm needed. (For r = 1 the modulus is p - 1 and the only
// residue is 1, recovering the ordinary Korselt condition n == 1 (mod p - 1).)
bool is_power_of_p_mod(const Integer& n, const Integer& p, unsigned long r) {
  Integer modulus = 1;
  for (unsigned long k = 0; k < r; ++k) modulus *= p;  // modulus = p^r
  modulus -= 1;                                        // modulus = p^r - 1

  const Integer target = n % modulus;
  Integer power = 1 % modulus;  // p^0, reduced (handles p^1 - 1 == 1, i.e. p = 2)
  for (unsigned long i = 0; i < r; ++i) {
    if (power == target) return true;
    power = (power * p) % modulus;
  }
  return false;
}

}  // namespace

bool is_carmichael(const Factorization& n, unsigned long order) {
  // Must be composite (>= 2 distinct primes) and squarefree.
  if (n.num_distinct_primes() < 2) return false;
  if (!n.is_squarefree()) return false;

  const Integer value = n.value();
  for (const PrimePower& pe : n.factors()) {
    for (unsigned long r = 1; r <= order; ++r) {
      if (!is_power_of_p_mod(value, pe.prime, r)) return false;
    }
  }
  return true;
}

bool is_rigid_carmichael(const Factorization& n, unsigned long order) {
  // Must be composite (>= 2 distinct primes) and squarefree.
  if (n.num_distinct_primes() < 2) return false;
  if (!n.is_squarefree()) return false;

  const Integer value = n.value();
  const Integer n_minus_one = value - 1;
  for (const PrimePower& pe : n.factors()) {
    const Integer& p = pe.prime;
    Integer modulus = 1;
    for (unsigned long r = 1; r <= order; ++r) {
      modulus *= p;  // modulus = p^r
      // Rigid condition: n == 1 (mod p^r - 1), i.e. (p^r - 1) | (n - 1).
      if (n_minus_one % (modulus - 1) != 0) return false;
    }
  }
  return true;
}

}  // namespace pseudoprime
