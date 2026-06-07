#include "pseudoprime/large_carmichael/layered_solver.hpp"

#include <algorithm>
#include <gtest/gtest.h>
#include <stdexcept>
#include <utility>
#include <vector>

#include "pseudoprime/carmichael.hpp"
#include "pseudoprime/factorization.hpp"
#include "pseudoprime/subset_product.hpp"

namespace {

using pseudoprime::Factorization;
using pseudoprime::Integer;
using pseudoprime::PrimePower;

std::vector<Integer> products_of(const pseudoprime::SubsetProductSolver& solver,
                                 const pseudoprime::SubsetProductProblem& p) {
  std::vector<Integer> out;
  for (const pseudoprime::SubsetProductSolution& s : solver.solve(p)) out.push_back(s.product);
  std::sort(out.begin(), out.end());
  return out;
}

// L = 80 = 2^4 * 5, with P(L) = {3, 11, 17, 41}; the lone solution is 561.
Factorization L80() { return Factorization({{Integer(2), 4}, {Integer(5), 1}}); }
pseudoprime::SubsetProductProblem problem_L80() {
  pseudoprime::SubsetProductProblem p;
  p.modulus = Integer(80);
  p.elements = {Integer(3), Integer(11), Integer(17), Integer(41)};
  p.labels = p.elements;
  p.target = Integer(1);
  p.min_size = 3;
  return p;
}

// L = 120 = 2^3 * 3 * 5, with P(L) = {7,11,13,31,41,61}; solutions 41041,172081,852841.
Factorization L120() { return Factorization({{Integer(2), 3}, {Integer(3), 1}, {Integer(5), 1}}); }
pseudoprime::SubsetProductProblem problem_L120() {
  pseudoprime::SubsetProductProblem p;
  p.modulus = Integer(120);
  p.elements = {Integer(7), Integer(11), Integer(13), Integer(31), Integer(41), Integer(61)};
  p.labels = p.elements;
  p.target = Integer(1);
  p.min_size = 3;
  return p;
}

// With brute_threshold >= the pool size, no pairing happens and the layered
// solver degenerates to an exhaustive brute force -- so it must agree with it.
TEST(LayeredSolver, HighThresholdEqualsBruteForce) {
  for (auto [fac, prob] : {std::pair{L80(), problem_L80()}, std::pair{L120(), problem_L120()}}) {
    pseudoprime::LayeredSubsetProductSolver layered(fac, pseudoprime::PrimePowerChainStrategy{},
                                                    /*brute_threshold=*/1000);
    pseudoprime::BruteForceSolver bf;
    EXPECT_EQ(products_of(layered, prob), products_of(bf, prob));
  }
}

// The default chain (balanced) must construct and solve correctly. With a high
// threshold no pairing runs, so it yields the full known set regardless of chain.
TEST(LayeredSolver, DefaultChainConstructsAndSolves) {
  pseudoprime::LayeredSubsetProductSolver layered(L120());  // default chain + threshold
  pseudoprime::SubsetProductProblem p = problem_L120();
  p.max_solutions = 0;
  pseudoprime::BruteForceSolver bf;
  // Default brute_threshold is 20 >= 6 elements, so this is exhaustive here.
  EXPECT_EQ(products_of(layered, p), products_of(bf, p));
}

TEST(LayeredSolver, DegenerateFindsKnownSolutions) {
  pseudoprime::LayeredSubsetProductSolver layered(L120(), pseudoprime::PrimePowerChainStrategy{},
                                                  /*brute_threshold=*/1000);
  EXPECT_EQ(products_of(layered, problem_L120()),
            (std::vector<Integer>{41041, 172081, 852841}));
}

// Every solution the layered solver returns -- with pairing actually engaged
// (low threshold) -- must be a genuine Carmichael number: product 1 mod L,
// distinct strictly-increasing indices, and Korselt-valid. (Heuristic: the count
// is not asserted; pairing may drop chaff.)
TEST(LayeredSolver, ReturnedSolutionsAreValid) {
  const auto p = problem_L120();
  pseudoprime::LayeredSubsetProductSolver layered(L120(), pseudoprime::PrimePowerChainStrategy{},
                                                  /*brute_threshold=*/2);
  for (const pseudoprime::SubsetProductSolution& s : layered.solve(p)) {
    EXPECT_GE(s.indices.size(), p.min_size);
    EXPECT_TRUE(std::is_sorted(s.indices.begin(), s.indices.end()));
    EXPECT_EQ(std::adjacent_find(s.indices.begin(), s.indices.end()), s.indices.end());  // distinct

    Integer prod = 1;
    for (std::size_t idx : s.indices) prod *= p.labels[idx];
    EXPECT_EQ(prod % p.modulus, Integer(1));
    EXPECT_EQ(prod, s.product);
    EXPECT_TRUE(pseudoprime::is_carmichael(Factorization::from_primes(s.labels))) << s.product;
  }
}

TEST(LayeredSolver, HonorsMaxSolutions) {
  auto p = problem_L120();
  pseudoprime::LayeredSubsetProductSolver layered(L120(), pseudoprime::PrimePowerChainStrategy{},
                                                  /*brute_threshold=*/1000);
  p.max_solutions = 1;
  EXPECT_EQ(products_of(layered, p).size(), 1u);
  p.max_solutions = 2;
  EXPECT_EQ(products_of(layered, p).size(), 2u);
  p.max_solutions = 0;
  EXPECT_EQ(products_of(layered, p).size(), 3u);
}

TEST(LayeredSolver, ThrowsOnModulusMismatch) {
  pseudoprime::LayeredSubsetProductSolver layered(L80());  // built for modulus 80
  auto p = problem_L120();                                 // but asked to solve mod 120
  EXPECT_THROW(
      {
        for ([[maybe_unused]] const auto& s : layered.solve(p)) {
        }
      },
      std::invalid_argument);
}

}  // namespace
