#include "pseudoprime/subset_product.hpp"

#include <algorithm>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

namespace {

using pseudoprime::Integer;

std::vector<Integer> solve_all_with(const pseudoprime::SubsetProductSolver& solver,
                                    const pseudoprime::SubsetProductProblem& p) {
  std::vector<Integer> products;
  for (const pseudoprime::SubsetProductSolution& s : solver.solve(p)) products.push_back(s.product);
  std::sort(products.begin(), products.end());
  return products;
}

std::vector<Integer> solve_all(const pseudoprime::SubsetProductProblem& p) {
  pseudoprime::BruteForceSolver solver;
  return solve_all_with(solver, p);
}

TEST(BruteForceSolver, FindsSubsetWithTargetProduct) {
  // Residues of {3,11,17,41} mod 80; only {3,11,17} multiplies to 1 mod 80.
  pseudoprime::SubsetProductProblem p;
  p.modulus = Integer(80);
  p.elements = {Integer(3), Integer(11), Integer(17), Integer(41)};
  p.labels = {Integer(3), Integer(11), Integer(17), Integer(41)};
  p.target = Integer(1);
  p.min_size = 3;
  EXPECT_EQ(solve_all(p), (std::vector<Integer>{561}));
}

TEST(BruteForceSolver, RespectsMinSize) {
  // All elements are 1 mod 10; with min_size 3 only the full triple qualifies.
  pseudoprime::SubsetProductProblem p;
  p.modulus = Integer(10);
  p.elements = {Integer(1), Integer(1), Integer(1)};
  p.labels = {Integer(11), Integer(21), Integer(31)};
  p.target = Integer(1);
  p.min_size = 3;
  const auto all = solve_all(p);
  ASSERT_EQ(all.size(), 1u);
  EXPECT_EQ(all[0], Integer(11) * Integer(21) * Integer(31));
}

TEST(BruteForceSolver, StopsAfterFirstWhenCallerBreaks) {
  pseudoprime::SubsetProductProblem p;
  p.modulus = Integer(10);
  p.elements = {Integer(1), Integer(1), Integer(1), Integer(1)};
  p.labels = {Integer(11), Integer(21), Integer(31), Integer(41)};
  p.target = Integer(1);
  p.min_size = 3;

  pseudoprime::BruteForceSolver solver;
  int count = 0;
  for ([[maybe_unused]] const auto& s : solver.solve(p)) {
    ++count;
    break;  // lazy generation: nothing further is computed
  }
  EXPECT_EQ(count, 1);
}

TEST(BruteForceSolver, ReportsIndicesAndLabels) {
  pseudoprime::SubsetProductProblem p;
  p.modulus = Integer(80);
  p.elements = {Integer(3), Integer(11), Integer(17), Integer(41)};
  p.labels = {Integer(3), Integer(11), Integer(17), Integer(41)};
  p.target = Integer(1);
  p.min_size = 3;

  pseudoprime::BruteForceSolver solver;
  pseudoprime::SubsetProductSolution found;
  for (const auto& s : solver.solve(p)) {
    found = s;
    break;
  }
  EXPECT_EQ(found.indices, (std::vector<std::size_t>{0, 1, 2}));
  EXPECT_EQ(found.labels, (std::vector<Integer>{3, 11, 17}));
  EXPECT_EQ(found.product, Integer(561));
}

TEST(BruteForceSolver, NoSolutionWhenTooFewElements) {
  pseudoprime::SubsetProductProblem p;
  p.modulus = Integer(80);
  p.elements = {Integer(3), Integer(11)};
  p.labels = {Integer(3), Integer(11)};
  p.target = Integer(1);
  p.min_size = 3;  // cannot reach size 3
  EXPECT_TRUE(solve_all(p).empty());
}

// Residues of P(120) = {7,11,13,31,41,61}; product == 1 mod 120, size >= 3.
pseudoprime::SubsetProductProblem problem_L120() {
  pseudoprime::SubsetProductProblem p;
  p.modulus = Integer(120);
  p.elements = {Integer(7), Integer(11), Integer(13), Integer(31), Integer(41), Integer(61)};
  p.labels = p.elements;
  p.target = Integer(1);
  p.min_size = 3;
  return p;
}

TEST(MeetInMiddleSolver, FindsSubsetWithTargetProduct) {
  pseudoprime::SubsetProductProblem p;
  p.modulus = Integer(80);
  p.elements = {Integer(3), Integer(11), Integer(17), Integer(41)};
  p.labels = p.elements;
  p.target = Integer(1);
  p.min_size = 3;
  pseudoprime::MeetInMiddleSolver mitm;
  EXPECT_EQ(solve_all_with(mitm, p), (std::vector<Integer>{561}));
}

TEST(MeetInMiddleSolver, AgreesWithBruteForce) {
  const auto p = problem_L120();
  pseudoprime::BruteForceSolver bf;
  pseudoprime::MeetInMiddleSolver mitm;
  EXPECT_EQ(solve_all_with(mitm, p), solve_all_with(bf, p));
  EXPECT_EQ(solve_all_with(mitm, p), (std::vector<Integer>{41041, 172081, 852841}));
}

TEST(MeetInMiddleSolver, StopsAfterFirstWhenCallerBreaks) {
  const auto p = problem_L120();
  pseudoprime::MeetInMiddleSolver mitm;
  int count = 0;
  for ([[maybe_unused]] const auto& s : mitm.solve(p)) {
    ++count;
    break;
  }
  EXPECT_EQ(count, 1);
}

TEST(MeetInMiddleSolver, ThrowsOnNonUnitElement) {
  pseudoprime::SubsetProductProblem p;
  p.modulus = Integer(10);
  p.elements = {Integer(3), Integer(5), Integer(7)};  // 5 shares a factor with 10
  p.labels = p.elements;
  p.target = Integer(1);
  p.min_size = 2;
  pseudoprime::MeetInMiddleSolver mitm;
  EXPECT_THROW(solve_all_with(mitm, p), std::invalid_argument);
}

// max_solutions caps the number of solutions yielded (0 == all). L=120 has three.
TEST(MaxSolutions, CapsBruteForce) {
  auto p = problem_L120();
  pseudoprime::BruteForceSolver bf;

  p.max_solutions = 0;
  EXPECT_EQ(solve_all_with(bf, p).size(), 3u);
  p.max_solutions = 1;
  EXPECT_EQ(solve_all_with(bf, p).size(), 1u);
  p.max_solutions = 2;
  EXPECT_EQ(solve_all_with(bf, p).size(), 2u);
  p.max_solutions = 99;  // cap above the count yields all
  EXPECT_EQ(solve_all_with(bf, p).size(), 3u);
}

TEST(MaxSolutions, CapsMeetInMiddle) {
  auto p = problem_L120();
  pseudoprime::MeetInMiddleSolver mitm;

  p.max_solutions = 0;
  EXPECT_EQ(solve_all_with(mitm, p).size(), 3u);
  p.max_solutions = 1;
  EXPECT_EQ(solve_all_with(mitm, p).size(), 1u);
  p.max_solutions = 2;
  EXPECT_EQ(solve_all_with(mitm, p).size(), 2u);
}

}  // namespace
