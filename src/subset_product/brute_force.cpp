#include "pseudoprime/subset_product/brute_force.hpp"

#include <cstddef>
#include <generator>
#include <vector>

#include "pseudoprime/modular_integer.hpp"

namespace pseudoprime {

namespace {

// DFS over include/exclude choices for elements[i..n). `acc` is the running
// product modulo the problem modulus (passed by value, so branches don't
// interfere); `chosen` is the current subset, shared by reference (push/pop
// bracket each recursive descent).
std::generator<SubsetProductSolution> search(const SubsetProductProblem& problem, std::size_t i,
                                             ModularInteger acc, const Integer& target,
                                             std::vector<std::size_t>& chosen) {
  if (i == problem.elements.size()) {
    if (chosen.size() >= problem.min_size && acc == target) {
      SubsetProductSolution sol;
      sol.indices = chosen;
      sol.product = 1;
      for (std::size_t idx : chosen) {
        sol.labels.push_back(problem.labels[idx]);
        sol.product *= problem.labels[idx];
      }
      co_yield sol;
    }
    co_return;
  }

  // Exclude elements[i].
  co_yield std::ranges::elements_of(search(problem, i + 1, acc, target, chosen));

  // Include elements[i].
  chosen.push_back(i);
  acc *= problem.elements[i];
  co_yield std::ranges::elements_of(search(problem, i + 1, acc, target, chosen));
  chosen.pop_back();
}

}  // namespace

std::generator<SubsetProductSolution> BruteForceSolver::solve(
    const SubsetProductProblem& problem) const {
  Integer target = problem.target % problem.modulus;
  if (target < 0) target += problem.modulus;
  std::vector<std::size_t> chosen;
  // Cap the number of solutions if requested (0 == enumerate all). The DFS is a
  // lazy generator, so stopping here also stops the search.
  std::size_t count = 0;
  for (auto&& sol : search(problem, 0, ModularInteger(problem.modulus), target, chosen)) {
    co_yield std::move(sol);
    if (problem.max_solutions != 0 && ++count >= problem.max_solutions) co_return;
  }
}

}  // namespace pseudoprime
