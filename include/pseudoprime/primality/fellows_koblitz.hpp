#pragma once

#include "pseudoprime/factorization.hpp"
#include "pseudoprime/primality/result.hpp"

namespace pseudoprime {

// Deterministic, unconditional primality test for `n`, given the complete
// factorization of `n - 1` (Bach & Shallit, Algorithm FELLOWS-KOBLITZ,
// Theorem 9.1.3; the deterministic n-1 test attributed to Pomerance).
//
// Sweeps bases a in [2, (log2 n)^2], computing each ord_n(a) from the known
// factorization of n-1; if h = lcm of those orders satisfies h >= sqrt(n) then
// every prime factor of n exceeds sqrt(n), forcing n prime; if h < sqrt(n) a
// smooth-number counting argument forces n composite. No randomness and no GRH.
// Cost: O((log n)^6 / log log n) bit operations. The first base (a = 2) is just
// a base-2 Fermat test, so most composites are rejected immediately.
//
// `n_minus_one` must factor `n - 1`; otherwise std::invalid_argument is thrown.
PrimalityCertificate is_prime_fellows_koblitz(const Integer& n,
                                              const Factorization& n_minus_one);

// Boolean convenience wrapper.
inline bool is_prime_fk(const Integer& n, const Factorization& n_minus_one) {
  return is_prime_fellows_koblitz(n, n_minus_one).verdict == Primality::Prime;
}

}  // namespace pseudoprime
