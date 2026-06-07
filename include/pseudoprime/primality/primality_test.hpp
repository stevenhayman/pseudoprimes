#pragma once

#include "pseudoprime/factorization.hpp"

namespace pseudoprime {

// Strategy for deciding the primality of n given the factorization of n - 1.
// Every implementation is fully deterministic (no probable-prime tests), so a
// `true`/`false` answer is a proof. Pluggable so the construction can trade
// speed for rigor, and so future tests (e.g. one consulting a precomputed prime
// table) can drop in behind the same interface.
class PrimalityTest {
 public:
  virtual ~PrimalityTest() = default;
  virtual bool is_prime(const Integer& n, const Factorization& n_minus_one) const = 0;
};

// The full Fellows-Koblitz test (unconditional, but sweeps ~(log n)^2 bases).
class FellowsKoblitzPrimalityTest : public PrimalityTest {
 public:
  bool is_prime(const Integer& n, const Factorization& n_minus_one) const override;
};

// Layered, deterministic, and ~1000x faster than full Fellows-Koblitz on typical
// inputs: a base-2 Fermat test rejects almost all composites in one step; a
// single-witness Lucas search proves almost all primes; the full Fellows-Koblitz
// test is invoked only for the rare residue (base-2 Fermat pseudoprimes, or
// primes whose least primitive root exceeds the Lucas bound).
class LayeredPrimalityTest : public PrimalityTest {
 public:
  explicit LayeredPrimalityTest(unsigned long lucas_max_base = 256)
      : lucas_max_base_(lucas_max_base) {}
  bool is_prime(const Integer& n, const Factorization& n_minus_one) const override;

 private:
  unsigned long lucas_max_base_;
};

}  // namespace pseudoprime
