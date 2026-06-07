#pragma once

#include "pseudoprime/factorization.hpp"

namespace pseudoprime {

// Multiplicative order of `a` modulo `n`: the least k > 0 with a^k == 1 (mod n).
//
// `group_order` is the factorization of a KNOWN multiple of the order -- in
// practice phi(n) or lambda(n), whose factorization the caller already has.
// The order is found by stripping prime factors from that multiple, so no
// factoring happens here.
//
// Precondition: gcd(a, n) == 1 and a^(group_order.value()) == 1 (mod n).
Integer multiplicative_order(const Integer& a, const Integer& n,
                             const Factorization& group_order);

}  // namespace pseudoprime
