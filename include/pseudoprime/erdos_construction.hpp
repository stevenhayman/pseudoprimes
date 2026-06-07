#pragma once

#include <cstddef>
#include <generator>
#include <vector>

#include "pseudoprime/factorization.hpp"
#include "pseudoprime/primality.hpp"
#include "pseudoprime/prime_set_strategy.hpp"  // PrimeSetStrategy
#include "pseudoprime/subset_product.hpp"

namespace pseudoprime {

// Build the order-m Erdos prime set
//   P(m,L) = { odd primes p : p does not divide L, and (p^r - 1) | L for all
//              1 <= r <= order }
// from the factored modulus L (Howe, "Higher-order Carmichael numbers", Sec. 3).
// Each divisor d of L gives a candidate p = d + 1 whose primality is certified
// using the factorization of p - 1 = d (available for free as a sub-factorization
// of L); this already secures (p - 1) | L, so for order > 1 only the extra
// conditions (p^r - 1) | L for r = 2..order are checked. order == 1 is the
// ordinary Erdos prime set P(L). The members are exactly those primes from which
// a subset-product == 1 (mod L) yields a rigid Carmichael number of order m.
std::vector<Integer> generate_prime_set(
    const Factorization& L, const PrimalityTest& primality = LayeredPrimalityTest{},
    unsigned long order = 1);

// A constructed Carmichael number, born factored: its prime factors are exactly
// the subset chosen by the solver.
struct ConstructedCarmichael {
  Factorization factorization;
  Integer value;
};

// Run the Erdos construction over the factored modulus L: build P(m,L), then
// yield each number whose prime factors form a subset of P(m,L) with product
// congruent to 1 modulo L (size >= min_factors), as found by `solver`. With
// order == 1 these are ordinary Carmichael numbers; with order > 1 they are
// rigid Carmichael numbers of order m (Howe's C(m,L)). Lazy: the caller picks
// find-first (break) or enumerate-all by how far it iterates. `max_solutions`
// (0 == all) is also passed to the solver as an internal cap.
std::generator<ConstructedCarmichael> construct_carmichael(const Factorization& L,
                                                           const SubsetProductSolver& solver,
                                                           std::size_t min_factors = 3,
                                                           unsigned long order = 1,
                                                           std::size_t max_solutions = 0);

}  // namespace pseudoprime
