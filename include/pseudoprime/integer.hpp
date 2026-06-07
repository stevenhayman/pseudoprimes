#pragma once

#include <gmpxx.h>

#include <optional>

namespace pseudoprime {

// Arbitrary-precision integer type used throughout the library.
using Integer = mpz_class;

// Thin, inline wrappers over GMP's C API so call sites read as mathematics
// rather than out-parameter plumbing. Each compiles to the same code as the
// raw mpz_* call, so there is no performance cost.

inline Integer pow_ui(const Integer& base, unsigned long exponent) {
  Integer r;
  mpz_pow_ui(r.get_mpz_t(), base.get_mpz_t(), exponent);
  return r;
}

inline Integer mod_pow(const Integer& base, const Integer& exp, const Integer& mod) {
  Integer r;
  mpz_powm(r.get_mpz_t(), base.get_mpz_t(), exp.get_mpz_t(), mod.get_mpz_t());
  return r;
}

inline Integer gcd(const Integer& a, const Integer& b) {
  Integer r;
  mpz_gcd(r.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  return r;
}

inline Integer lcm(const Integer& a, const Integer& b) {
  Integer r;
  mpz_lcm(r.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  return r;
}

// True iff n is a perfect square (n = k*k for some integer k >= 0). Negative
// inputs are never perfect squares, so this returns false for them.
inline bool is_perfect_square(const Integer& n) {
  return mpz_perfect_square_p(n.get_mpz_t()) != 0;
}

// The inverse of a modulo m, or nullopt if a is not invertible modulo m.
inline std::optional<Integer> mod_inverse(const Integer& a, const Integer& m) {
  Integer r;
  if (!mpz_invert(r.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t())) return std::nullopt;
  return r;
}

}  // namespace pseudoprime
