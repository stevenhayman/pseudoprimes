#pragma once

#include "pseudoprime/integer.hpp"

namespace pseudoprime::detail {

// Does (p^r - 1) | L for every r with 2 <= r <= order? The r == 1 case,
// (p - 1) | L, is guaranteed by construction (p - 1 is a divisor of L), so it is
// not retested here. Returns true immediately for order <= 1. Shared by the
// divisor and sieve prime-set strategies.
bool satisfies_higher_order(const Integer& p, const Integer& L_value, unsigned long order);

}  // namespace pseudoprime::detail
