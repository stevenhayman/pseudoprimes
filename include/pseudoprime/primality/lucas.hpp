#pragma once

#include "pseudoprime/factorization.hpp"
#include "pseudoprime/primality/result.hpp"

namespace pseudoprime {

// Search bases a in [2, max_base] for a primitive root modulo p, using the
// factorization of p - 1: a witness satisfies a^(p-1) == 1 (mod p) and
// gcd(a^((p-1)/q) - 1, p) == 1 for every prime q | (p-1).
//
// Returns Primality::Prime iff a witness is found -- a complete, deterministic
// PROOF that p is prime (Lehmer's converse of Fermat). Otherwise returns
// Primality::Inconclusive: no witness was found within the bound (p might be a
// prime whose least primitive root exceeds max_base), which is NOT a proof of
// compositeness. This test never returns Primality::Composite; a caller must
// resolve an Inconclusive result by other means (e.g. is_prime_fellows_koblitz).
//
// `p_minus_one` must factor p - 1.
Primality lucas_certifies_prime(const Integer& p, const Factorization& p_minus_one,
                                unsigned long max_base = 256);

}  // namespace pseudoprime
