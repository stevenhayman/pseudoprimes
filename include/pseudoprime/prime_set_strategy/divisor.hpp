#pragma once

#include <generator>

#include "pseudoprime/prime_set_strategy/prime_set_strategy.hpp"

namespace pseudoprime {

// Default strategy: stream every divisor d of L (born factored, via
// Factorization::divisors_lazy) and test p = d + 1. Exact and complete -- this
// is the reference behaviour that generate_prime_set has always had, and it is
// the oracle the sieve strategy is checked against.
class DivisorPrimeSetStrategy : public PrimeSetStrategy {
 public:
  std::generator<Integer> generate(const Factorization& L, const PrimalityTest& primality,
                                   unsigned long order) const override;
};

}  // namespace pseudoprime
