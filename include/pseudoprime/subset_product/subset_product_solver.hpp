#pragma once

#include <cstddef>
#include <generator>
#include <vector>

#include "pseudoprime/factorization.hpp"

namespace pseudoprime {

// Find subsets of `elements` whose product is congruent to `target` modulo
// `modulus`, with at least `min_size` members. In the Erdos construction the
// elements are the residues p mod L of the candidate primes, `labels` are the
// primes themselves, `modulus` is L, and `target` is 1.
struct SubsetProductProblem {
  Integer modulus;                // L
  std::vector<Integer> elements;  // residues to multiply (units mod L)
  std::vector<Integer> labels;    // parallel array: the primes p_i
  Integer target = 1;             // desired residue of the product
  std::size_t min_size = 3;       // minimum number of chosen elements
  std::size_t max_solutions = 0;  // stop after this many solutions; 0 == enumerate all
};

struct SubsetProductSolution {
  std::vector<std::size_t> indices;  // chosen positions into elements/labels
  std::vector<Integer> labels;       // the chosen primes
  Integer product;                   // product of the chosen labels (a Carmichael number)
};

// Strategy interface so the search algorithm is pluggable (brute force and
// meet-in-the-middle now; Shallue pairing, back-propagated search later). Each
// solver yields solutions lazily, so callers choose find-first (break / take(1))
// or enumerate-all simply by how far they iterate.
class SubsetProductSolver {
 public:
  virtual ~SubsetProductSolver() = default;
  virtual std::generator<SubsetProductSolution> solve(
      const SubsetProductProblem& problem) const = 0;
};

}  // namespace pseudoprime
