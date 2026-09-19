#pragma once

#include <generator>

#include "pseudoprime/subset_product/subset_product_solver.hpp"

namespace pseudoprime {

// Meet-in-the-middle. Splits the elements into two halves, tabulates the modular
// product of every subset of the first half, then for each subset of the second
// half looks up the complementary residue target * (prod_2)^{-1} (mod modulus).
// Time and space O(2^(n/2)) instead of O(2^n).
//
// Requires every element to be a unit modulo `modulus` (so half-products are
// invertible); this holds in the Erdos construction, where elements are primes
// not dividing L. Throws std::invalid_argument otherwise.
class MeetInMiddleSolver : public SubsetProductSolver {
 public:
  std::generator<SubsetProductSolution> solve(const SubsetProductProblem& problem) const override;
};

}  // namespace pseudoprime
