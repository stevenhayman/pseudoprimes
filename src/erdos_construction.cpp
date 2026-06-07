#include "pseudoprime/erdos_construction.hpp"

#include <algorithm>

#include "pseudoprime/prime_set_strategy.hpp"

namespace pseudoprime {

std::vector<Integer> generate_prime_set(const Factorization& L, const PrimalityTest& primality,
                                        unsigned long order) {
  // Delegate to the default (divisor-enumeration) strategy and canonicalize the
  // order by prime, so results stay deterministic regardless of the strategy's
  // emission order. Other strategies (e.g. SievePrimeSetStrategy) can be used
  // directly via the PrimeSetStrategy interface.
  const DivisorPrimeSetStrategy strategy;
  std::vector<Integer> result;
  for (Integer&& p : strategy.generate(L, primality, order)) result.push_back(std::move(p));

  std::sort(result.begin(), result.end(),
            [](const Integer& a, const Integer& b) { return a < b; });
  return result;
}

std::generator<ConstructedCarmichael> construct_carmichael(const Factorization& L,
                                                           const SubsetProductSolver& solver,
                                                           std::size_t min_factors,
                                                           unsigned long order,
                                                           std::size_t max_solutions) {
  const std::vector<Integer> prime_set = generate_prime_set(L, LayeredPrimalityTest{}, order);

  const Integer Lval = L.value();
  SubsetProductProblem problem;
  problem.modulus = Lval;
  problem.target = 1;
  problem.min_size = min_factors;
  problem.max_solutions = max_solutions;
  problem.elements.reserve(prime_set.size());
  problem.labels.reserve(prime_set.size());
  for (const Integer& p : prime_set) {
    problem.elements.push_back(p % Lval);  // p mod L, the subset-product element
    problem.labels.push_back(p);
  }

  for (const SubsetProductSolution& sol : solver.solve(problem)) {
    ConstructedCarmichael cc;
    // The chosen primes are distinct, so the result is squarefree by construction.
    cc.factorization = Factorization::from_primes(sol.labels);
    cc.value = sol.product;
    co_yield cc;
  }
}

}  // namespace pseudoprime
