#pragma once

#include "pseudoprime/factorization.hpp"

namespace pseudoprime {

// Carmichael number of order m (Howe, "Higher-order Carmichael numbers"): a
// composite squarefree n such that for every prime divisor p of n and every r
// with 1 <= r <= order, there is an integer i >= 0 with n == p^i (mod p^r - 1).
// Order 1 is the ordinary Korselt criterion: squarefree, composite, and
// (p - 1) | (n - 1) for every prime p | n. (These force n odd with >= 3 prime
// factors.) An order-m Carmichael number is also a Carmichael number of every
// smaller order, down to order 1.
//
// This is the GENERAL notion: it admits non-rigid numbers (where some i != 0).
// Takes the known factorization so no factoring is performed.
bool is_carmichael(const Factorization& n, unsigned long order = 1);

// Rigid Carmichael number of order m (Howe, Section 5): the i = 0 special case
// of the above -- a composite squarefree n with n == 1 (mod p^r - 1) for every
// prime divisor p of n and every r with 1 <= r <= order. At order 1 this
// coincides with is_carmichael. Every element of the C(m,L) construction is
// rigid; not every order-m Carmichael number is.
bool is_rigid_carmichael(const Factorization& n, unsigned long order = 1);

}  // namespace pseudoprime
