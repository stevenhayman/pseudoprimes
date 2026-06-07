#pragma once

#include <cstddef>
#include <vector>

#include "pseudoprime/factorization.hpp"

namespace pseudoprime {

// This is the paper's omega function (Definition 1.2) -- a map from a unit a in
// (Z/LZ)^x to an index, NOT the standard number-theoretic omega(n) (= count of
// distinct prime factors). For L = prod Q_i (the Q_i = q_i^{h_i} are the
// prime-power components of L, ascending -- see Factorization::prime_powers),
// omega(a) is the LARGEST 1-based index i with a mod Q_i != 1, or 0 if a == 1
// mod Q_i for every i. So omega(a) <= m exactly when a lies in the subgroup
// (Z/(Q_1...Q_m)Z)^x.
//
// (The paper's Definition 1.3 omega-bar uses the SMALLEST such index and drives
// Algorithm 1 -- constructing many Carmichael numbers -- which is out of scope
// here. We implement only the max-index version for Algorithm 2.)
std::size_t weight(const Integer& a, const std::vector<Integer>& components);

}  // namespace pseudoprime
