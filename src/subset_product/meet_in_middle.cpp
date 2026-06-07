#include "pseudoprime/subset_product.hpp"

#include <cstdint>
#include <generator>
#include <map>
#include <optional>
#include <stdexcept>
#include <vector>

#include "pseudoprime/modular_product.hpp"

namespace pseudoprime {

namespace {

// Reduce v into [0, m).
Integer mod_norm(const Integer& v, const Integer& m) {
  Integer r = v % m;
  if (r < 0) r += m;
  return r;
}

}  // namespace

std::generator<SubsetProductSolution> MeetInMiddleSolver::solve(
    const SubsetProductProblem& problem) const {
  const std::size_t n = problem.elements.size();
  const Integer& L = problem.modulus;
  const Integer target = mod_norm(problem.target, L);

  // Precondition: every element is a unit mod L, so any half-product is
  // invertible and the complementary lookup is well defined.
  for (const Integer& e : problem.elements) {
    if (gcd(e, L) != 1) {
      throw std::invalid_argument(
          "MeetInMiddleSolver requires every element to be coprime to the modulus");
    }
  }

  const std::size_t mid = n / 2;
  const std::size_t rest = n - mid;
  if (mid >= 63 || rest >= 63) {
    throw std::invalid_argument("MeetInMiddleSolver: too many elements for a 64-bit subset mask");
  }

  // Tabulate the first half: residue -> the half-subset products achieving it.
  std::map<Integer, std::vector<ModularProduct>> table;
  for (std::uint64_t mask = 0; mask < (std::uint64_t{1} << mid); ++mask) {
    ModularProduct mp(L);
    for (std::size_t i = 0; i < mid; ++i)
      if ((mask >> i) & 1u) mp.multiply(problem.elements[i], i);
    table[mp.residue()].push_back(std::move(mp));
  }

  // Sweep the second half and look up the complementary residue. `count` caps
  // the number of solutions yielded (problem.max_solutions; 0 == enumerate all).
  std::size_t count = 0;
  for (std::uint64_t mask = 0; mask < (std::uint64_t{1} << rest); ++mask) {
    ModularProduct mp_b(L);
    for (std::size_t j = 0; j < rest; ++j)
      if ((mask >> j) & 1u) mp_b.multiply(problem.elements[mid + j], mid + j);

    // Need a half-product whose residue is target * (mp_b residue)^{-1}.
    const std::optional<Integer> inv_b = mod_inverse(mp_b.residue(), L);
    if (!inv_b) continue;  // unreachable given the unit precondition, but stay safe
    const Integer need = mod_norm(target * *inv_b, L);

    const auto it = table.find(need);
    if (it == table.end()) continue;

    for (const ModularProduct& mp_a : it->second) {
      if (mp_a.size() + mp_b.size() < problem.min_size) continue;

      // First-half indices precede second-half indices, so history stays sorted.
      const ModularProduct full = mp_a.combine(mp_b);
      SubsetProductSolution sol;
      sol.indices = full.elements();
      sol.product = 1;
      for (std::size_t k : sol.indices) {
        sol.labels.push_back(problem.labels[k]);
        sol.product *= problem.labels[k];
      }
      co_yield sol;
      if (problem.max_solutions != 0 && ++count >= problem.max_solutions) co_return;
    }
  }
}

}  // namespace pseudoprime
