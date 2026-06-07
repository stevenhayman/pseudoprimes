#pragma once

#include "pseudoprime/factorization.hpp"
#include "pseudoprime/primality/result.hpp"

namespace pseudoprime {

// Brillhart-Lehmer-Selfridge (BLS) N-1 primality test (BLS 1975, Theorems 5 and
// 11; "Theorem 3.1").
//
// Unlike lucas_certifies_prime and is_prime_fellows_koblitz -- which need the
// COMPLETE factorization of n-1 -- BLS proves primality from a PARTIAL factored
// part. `f` is a fully-factored divisor F of n-1 with gcd(F, (n-1)/F) == 1; the
// cofactor R = (n-1)/F need not be factored. This is BLS's defining advantage:
// a primality proof is possible as soon as the easy-to-factor part of n-1 is
// large enough, even when the rest is intractable.
//
// Let N = n, N - 1 = F * R. For each prime q | F we seek a witness a_q with
// a_q^(N-1) == 1 (mod N) and gcd(a_q^((N-1)/q) - 1, N) == 1. Given such
// witnesses:
//   * if F > sqrt(N)            -> N is prime (Pocklington branch); else
//   * if F > N^(1/3): write R = 2*m*F + s (0 <= s < 2F) and set D = s^2 - 8*m;
//     if m == 0 or D is not a perfect square (in particular if D < 0) then N is
//     prime, otherwise N is composite.
//
// Returns:
//   * Primality::Prime       -- a proof that n is prime.
//   * Primality::Composite   -- a proof that n is composite (a Fermat or gcd
//                               witness failed, or the discriminant is square).
//                               May carry a nontrivial factor (self-witnessing).
//   * Primality::Inconclusive-- not proven: F is too small (F^3 <= n) to certify,
//                               or no witness was found within [2, max_base].
//                               NOT evidence of compositeness.
//
// Throws std::invalid_argument if F does not divide n-1, or gcd(F, (n-1)/F) != 1,
// or F <= 1.
PrimalityCertificate bls_certifies_prime(const Integer& n, const Factorization& f,
                                         unsigned long max_base = 256);

}  // namespace pseudoprime
