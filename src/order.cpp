#include "pseudoprime/order.hpp"

namespace pseudoprime {

Integer multiplicative_order(const Integer& a, const Integer& n,
                             const Factorization& group_order) {
  // Start from the supplied multiple m of ord(a); for each prime q | m, remove
  // as many factors of q as still leave a^(ord) == 1 (mod n).
  Integer ord = group_order.value();
  for (const PrimePower& pe : group_order.factors()) {
    for (unsigned long k = 0; k < pe.exponent; ++k) {
      Integer candidate = ord / pe.prime;  // exact while q still divides ord
      if (mod_pow(a, candidate, n) == 1) {
        ord = candidate;
      } else {
        break;
      }
    }
  }
  return ord;
}

}  // namespace pseudoprime
