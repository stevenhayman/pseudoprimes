#pragma once

#include <generator>

#include "pseudoprime/subset_product/subset_product_solver.hpp"

namespace pseudoprime {

// Exhaustive depth-first search over include/exclude choices, maintaining a
// running product modulo `modulus`. O(2^n); the right starting point and a
// correctness baseline for the faster strategies.
class BruteForceSolver : public SubsetProductSolver {
 public:
  std::generator<SubsetProductSolution> solve(const SubsetProductProblem& problem) const override;
};

}  // namespace pseudoprime
