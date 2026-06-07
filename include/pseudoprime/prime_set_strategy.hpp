#pragma once

#include <cstddef>
#include <cstdint>
#include <generator>

#include "pseudoprime/factorization.hpp"
#include "pseudoprime/integer.hpp"
#include "pseudoprime/primality.hpp"

namespace pseudoprime {

// Strategy for producing the order-m Erdos prime set
//   P(m,L) = { odd primes p : p does not divide L, and (p^r - 1) | L for all
//              1 <= r <= order }.
//
// Each yielded value is a member prime p: an odd prime with (p-1) | L and p not
// dividing L (plus the higher-order conditions for order > 1). Pluggable so
// callers can swap the generation algorithm, mirroring the project's other
// strategy interfaces (PrimalityTest, SubsetProductSolver, SubgroupChainStrategy).
// Every implementation is EXACT: each yielded p is proven to satisfy all order-m
// conditions (deterministic primality, no probable-prime answers). The interface
// does NOT promise any particular emission order; callers that need a canonical
// order sort the collected result.
class PrimeSetStrategy {
 public:
  virtual ~PrimeSetStrategy() = default;
  virtual std::generator<Integer> generate(const Factorization& L,
                                           const PrimalityTest& primality,
                                           unsigned long order) const = 0;
};

// Default strategy: stream every divisor d of L (born factored, via
// Factorization::divisors_lazy) and test p = d + 1. Exact and complete -- this
// is the reference behaviour that generate_prime_set has always had, and it is
// the oracle the sieve strategy is checked against.
class DivisorPrimeSetStrategy : public PrimeSetStrategy {
 public:
  std::generator<Integer> generate(const Factorization& L, const PrimalityTest& primality,
                                   unsigned long order) const override;
};

// Sieve strategy: enumerate the divisor-candidates p = d + 1 by *value window*
// using a segmented factorizing sieve over supp(L) = { primes dividing L },
// rather than by divisor-tree order. A segment of consecutive d-values is sieved
// by the small primes of L (with their exponent caps) to identify, exactly, the
// d in that window that divide L; each such divisor is then confirmed with the
// same deterministic primality and (p^r - 1) | L checks as the divisor strategy.
//
// The sieve is bounded by a value cap B on the candidate prime p (equivalently
// on d = p - 1): its cost scales with the swept range and the divisors found in
// it, INDEPENDENT of the total divisor count d(L), so it reaches much larger /
// richer L than full divisor enumeration. It yields in ascending p order.
//
// With cap == 0 (or any cap at/above the per-order ceiling) it is byte-identical
// to DivisorPrimeSetStrategy: no order-m member can have p^m - 1 > L, so every
// member has p <= (L+1)^(1/m) and the bounded sweep still finds them all. With a
// smaller cap it returns exactly the members with p <= cap + 1 (exact for what it
// returns; "complete up to the cap").
//
// Construction parameters:
//   cap    -- upper bound B on candidate values d = p - 1. 0 means "complete":
//             cap at the exact per-order ceiling provable_cap(L, order) =
//             ceil((L+1)^(1/order)), beyond which no member can exist. The sieve
//             runs in 64-bit words, so this requires that ceiling to fit in 64
//             bits (for order 1 it is L + 1); otherwise pass an explicit cap.
//   window -- segment width: number of candidate values sieved per pass. Tunes
//             the memory/cache footprint (one machine word per slot).
class SievePrimeSetStrategy : public PrimeSetStrategy {
 public:
  explicit SievePrimeSetStrategy(std::uint64_t cap = 0, std::size_t window = std::size_t{1} << 18)
      : cap_(cap), window_(window == 0 ? 1 : window) {}

  std::generator<Integer> generate(const Factorization& L, const PrimalityTest& primality,
                                   unsigned long order) const override;

 private:
  std::uint64_t cap_;
  std::size_t window_;
};

// The exact per-order ceiling on candidate values for the sieve strategy:
//   provable_cap(L, m) = ceil((L+1)^(1/m)).
// Since an order-m member p satisfies p^m - 1 <= L, every member has
// p <= (L+1)^(1/m), so sweeping d = p - 1 up to this value is COMPLETE -- it is
// the smallest cap derived purely from L and m that misses nothing (no tuning
// constant). It can be large (for order 1 it is L + 1); members typically
// concentrate far below it, so callers trading completeness for speed may pass a
// smaller explicit cap to SievePrimeSetStrategy.
Integer provable_cap(const Factorization& L, unsigned long order);

}  // namespace pseudoprime
