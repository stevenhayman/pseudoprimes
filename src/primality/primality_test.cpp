#include "pseudoprime/primality/primality_test.hpp"

#include "pseudoprime/primality/fellows_koblitz.hpp"
#include "pseudoprime/primality/lucas.hpp"

namespace pseudoprime {

bool FellowsKoblitzPrimalityTest::is_prime(const Integer& n,
                                           const Factorization& n_minus_one) const {
  return is_prime_fellows_koblitz(n, n_minus_one).verdict == Primality::Prime;
}

bool LayeredPrimalityTest::is_prime(const Integer& n, const Factorization& n_minus_one) const {
  if (n < 2) return false;
  if (n == 2 || n == 3) return true;
  if (mpz_even_p(n.get_mpz_t())) return false;

  // Base-2 Fermat: a definite composite proof for almost all composites (1 modexp).
  if (mod_pow(Integer(2), n - 1, n) != 1) return false;

  // Lucas: a definite primality proof for almost all primes (a few modexps).
  if (lucas_certifies_prime(n, n_minus_one, lucas_max_base_) == Primality::Prime) return true;

  // Rare residue (base-2 Fermat pseudoprimes; primes with a large least
  // primitive root): resolve definitively with the full test.
  return is_prime_fellows_koblitz(n, n_minus_one).verdict == Primality::Prime;
}

}  // namespace pseudoprime
