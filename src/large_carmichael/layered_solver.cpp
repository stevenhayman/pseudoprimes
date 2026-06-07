#include "pseudoprime/large_carmichael/layered_solver.hpp"

#include <algorithm>
#include <cstddef>
#include <generator>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

#include "pseudoprime/large_carmichael/pairing.hpp"
#include "pseudoprime/modular_integer.hpp"
#include "pseudoprime/modular_product.hpp"

namespace pseudoprime {

namespace {

// DFS over include/exclude choices for survivors[i..n). `acc` is the running
// product modulo the problem modulus (by value, so branches don't interfere);
// `chosen` is the current list of ORIGINAL element indices (each survivor
// contributes its whole index list), shared by reference with push/pop
// bracketing. Mirrors BruteForceSolver's search, but over paired survivors.
std::generator<SubsetProductSolution> residual_search(
    const std::vector<ModularProduct>& survivors, std::size_t i, ModularInteger acc,
    const Integer& target, std::size_t min_size, std::vector<std::size_t>& chosen,
    const std::vector<Integer>& labels) {
  if (i == survivors.size()) {
    if (chosen.size() >= min_size && acc == target) {
      SubsetProductSolution sol;
      sol.indices = chosen;
      std::sort(sol.indices.begin(), sol.indices.end());
      sol.product = 1;
      for (std::size_t idx : sol.indices) {
        sol.labels.push_back(labels[idx]);
        sol.product *= labels[idx];
      }
      co_yield sol;
    }
    co_return;
  }

  // Exclude survivor i.
  co_yield std::ranges::elements_of(
      residual_search(survivors, i + 1, acc, target, min_size, chosen, labels));

  // Include survivor i: fold in its residue and all of its original indices.
  acc *= survivors[i].residue();
  const std::size_t added = survivors[i].elements().size();
  for (std::size_t e : survivors[i].elements()) chosen.push_back(e);
  co_yield std::ranges::elements_of(
      residual_search(survivors, i + 1, acc, target, min_size, chosen, labels));
  for (std::size_t k = 0; k < added; ++k) chosen.pop_back();
}

}  // namespace

LayeredSubsetProductSolver::LayeredSubsetProductSolver(const Factorization& modulus,
                                                       const SubgroupChainStrategy& chain_strategy,
                                                       std::size_t brute_threshold,
                                                       std::size_t max_pairing)
    : modulus_(modulus.value()),
      chain_(chain_strategy.build(modulus)),
      brute_threshold_(brute_threshold),
      max_pairing_(max_pairing) {}

std::generator<SubsetProductSolution> LayeredSubsetProductSolver::solve(
    const SubsetProductProblem& problem) const {
  if (problem.modulus != modulus_) {
    throw std::invalid_argument(
        "LayeredSubsetProductSolver: problem.modulus differs from the factored modulus it was "
        "constructed with");
  }
  const Integer& L = problem.modulus;
  Integer target = problem.target % L;
  if (target < 0) target += L;

  // Initial pool: one ModularProduct per element, tagged with its index.
  std::vector<ModularProduct> pool;
  pool.reserve(problem.elements.size());
  for (std::size_t i = 0; i < problem.elements.size(); ++i) {
    pool.emplace_back(ModularInteger(problem.elements[i], L), std::vector<std::size_t>{i});
  }

  // Phase 1: pair the pool down level by level until it is small enough for an
  // exact search, or the chain is exhausted (the default "pair all levels").
  for (std::size_t level = 0; level < chain_->levels(); ++level) {
    if (pool.size() <= brute_threshold_) break;
    pool = detail::pair_pool(pool, *chain_, level, max_pairing_);
  }

  // Phase 2: enumerate survivor sub-collections whose product == target,
  // honoring max_solutions (0 == enumerate all of the reachable solutions).
  std::vector<std::size_t> chosen;
  std::size_t count = 0;
  for (auto&& sol : residual_search(pool, 0, ModularInteger(L), target, problem.min_size, chosen,
                                    problem.labels)) {
    co_yield std::move(sol);
    if (problem.max_solutions != 0 && ++count >= problem.max_solutions) co_return;
  }
}

}  // namespace pseudoprime
