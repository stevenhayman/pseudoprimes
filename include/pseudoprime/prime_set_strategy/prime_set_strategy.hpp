#pragma once

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

}  // namespace pseudoprime
