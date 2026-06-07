#pragma once

#include <cstddef>
#include <generator>
#include <vector>

#include "pseudoprime/factorization.hpp"

namespace pseudoprime {

// Generators of richly composite moduli L, produced directly as factorizations
// (born factored -- no factoring needed). These feed straight into the Erdos
// path: the prime set P(m,L) is built from the divisors of L, so |P(m,L)| is
// bounded by d(L). The numbers below are exactly the L that maximize the divisor
// count for their size, which makes them the natural moduli for studying how the
// prime set grows.

// Build 2^e[0] * 3^e[1] * 5^e[2] * ... from an exponent vector over the first
// e.size() primes. Trailing/interior zero exponents are dropped. This is the
// low-level constructor; the exponents need not be non-increasing.
Factorization factorization_from_exponents(const std::vector<unsigned long>& exponents);

// Superior highly composite numbers (Ramanujan; OEIS A002201). For a parameter
// eps > 0 the exponent of prime p is floor(1/(p^eps - 1)); as eps decreases the
// sequence 2, 6, 12, 60, 120, 360, 2520, 5040, ... is swept out, each term
// differing from the previous by a single prime factor. They are streamed in
// increasing value order, born factored, by popping the largest "critical eps"
// event from a per-prime queue. count == 0 streams unboundedly (the caller
// breaks); otherwise at most `count` terms are produced. Every SHCN is a highly
// composite number, so this is a (sparse) subsequence of highly_composite below.
std::generator<Factorization> superior_highly_composite(std::size_t count = 0);

// Candidate highly composite numbers: integers 2^a1 * 3^a2 * ... over an initial
// run of primes with NON-INCREASING exponents a1 >= a2 >= ... >= 1 (Ramanujan's
// necessary form for a highly composite number). Streamed in increasing value
// order. The condition is necessary but not sufficient, so this is a superset of
// the true HCN sequence; every highly composite number appears here.
std::generator<Factorization> candidate_highly_composite();

// True highly composite numbers (OEIS A002182): 1, 2, 4, 6, 12, 24, 36, 48, ...
// -- the integers n with d(n) greater than d(m) for every m < n. Obtained by
// filtering candidate_highly_composite (which is in value order) to the running
// records of the divisor count d(n); since every record-setter is of candidate
// form, this filter reproduces the HCN sequence exactly, with no separate search.
std::generator<Factorization> highly_composite();

}  // namespace pseudoprime
