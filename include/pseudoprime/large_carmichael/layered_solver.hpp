#pragma once

#include <cstddef>
#include <generator>
#include <memory>

#include "pseudoprime/factorization.hpp"
#include "pseudoprime/large_carmichael/subgroup_chain.hpp"
#include "pseudoprime/subset_product.hpp"

namespace pseudoprime {

// A subset-product solver that scales to large element pools by first PAIRING
// the pool down through a subgroup chain of (Z/modulus Z)^x -- reusing the
// large-Carmichael pairing engine (pseudoprime/large_carmichael/pairing.hpp) --
// and then BRUTE-FORCING the small residual set of survivors.
//
// This is the answer to "subset-product enumeration does not scale": brute force
// is O(2^n) and meet-in-the-middle is O(2^(n/2)) in time and space, both hopeless
// once |P(L)| is in the thousands. Pairing combines elements whose coordinates
// cancel level by level, shrinking the pool toward the identity (the target,
// 1 mod L, for Carmichael construction); once it is at most `brute_threshold`,
// the residual is enumerated exactly.
//
// CONTRACT: this solver is HEURISTIC. Pairing drops unmatched leftovers ("chaff")
// at each level, so it finds SOME solutions, not provably all -- unlike
// BruteForceSolver / MeetInMiddleSolver, which are exhaustive. It is tuned for
// target == 1 (pairing drives toward the identity); other targets are still
// checked correctly in the residual search but benefit less from the pairing.
// With brute_threshold >= the pool size no pairing happens and the solver
// degenerates to an exhaustive brute force.
//
// The chain depends on the modulus's FACTORIZATION (never re-factored, per the
// project's policy), so the solver is constructed with the factored modulus and
// builds its chain eagerly. solve()'s problem.modulus must equal modulus.value().
class LayeredSubsetProductSolver : public SubsetProductSolver {
 public:
  // The default chain is the BALANCED congruence chain: for the rich moduli this
  // solver targets (|P(L)| >> 4^sqrt(log|G|)) its fewer levels deplete the pairing
  // pool far less than the prime-power chain, which is what lets the pairing reach
  // a solution. (The large-Carmichael construction keeps prime-power as its
  // default, where it empirically wins at smaller sizes.)
  explicit LayeredSubsetProductSolver(
      const Factorization& modulus,
      const SubgroupChainStrategy& chain_strategy = BalancedCongruenceChainStrategy{},
      std::size_t brute_threshold = 20, std::size_t max_pairing = 2);

  std::generator<SubsetProductSolution> solve(const SubsetProductProblem& problem) const override;

 private:
  Integer modulus_;                       // == the factored modulus's value
  std::unique_ptr<SubgroupChain> chain_;  // built once from the factorization
  std::size_t brute_threshold_;
  std::size_t max_pairing_;
};

}  // namespace pseudoprime
