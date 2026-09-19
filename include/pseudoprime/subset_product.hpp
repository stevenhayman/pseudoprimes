#pragma once
// Umbrella header for the subset-product solvers. Include this to get every
// solver, or include an individual header under pseudoprime/subset_product/ for
// just one.
//
// A third solver, `LayeredSubsetProductSolver`, lives in
// "pseudoprime/large_carmichael/layered_solver.hpp": it pairs the element pool
// down through a subgroup chain (reusing the large-Carmichael pairing engine)
// and then brute-forces the small residual. It is declared there rather than
// here because it depends on the subgroup-chain machinery; it is heuristic
// (finds some solutions, not provably all), unlike the two solvers included here.
#include "pseudoprime/subset_product/brute_force.hpp"
#include "pseudoprime/subset_product/meet_in_middle.hpp"
#include "pseudoprime/subset_product/subset_product_solver.hpp"
