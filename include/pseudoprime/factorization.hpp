#pragma once

#include <cstddef>
#include <generator>
#include <utility>
#include <vector>

#include "pseudoprime/integer.hpp"

namespace pseudoprime {

// A prime raised to a positive power. `exponent` is always >= 1.
struct PrimePower {
  Integer prime;
  unsigned long exponent;
};

// The canonical representation of an integer whose factorization is KNOWN.
//
// Passing a Factorization (rather than a bare Integer) through a function
// signature is how we guarantee a number is never re-factored: anywhere the
// prime structure is needed, the caller hands over what it already knows.
//
// Invariant: `factors_` is sorted ascending by prime, primes are distinct, and
// every exponent is >= 1. An empty factorization represents the integer 1.
//
// NOTE: constructors trust the caller that the supplied bases are prime. This
// type records a factorization; it does not discover or verify one.
class Factorization {
 public:
  Factorization() = default;  // the empty product == 1

  explicit Factorization(std::vector<PrimePower> factors);

  static Factorization from_prime_power(const Integer& prime, unsigned long exponent = 1);

  // Build from prime factors with multiplicity (duplicates are merged).
  static Factorization from_primes(const std::vector<Integer>& primes);

  const std::vector<PrimePower>& factors() const { return factors_; }
  bool is_one() const { return factors_.empty(); }
  std::size_t num_distinct_primes() const { return factors_.size(); }

  Integer value() const;        // reconstruct n = product p^e
  Integer num_divisors() const;  // d(n) = product (e + 1); no enumeration
  bool is_squarefree() const;   // every exponent == 1
  std::vector<Integer> primes() const;                 // distinct primes
  std::vector<Integer> prime_powers() const;           // the components p^e, ascending by prime
  unsigned long multiplicity(const Integer& p) const;  // exponent of p, else 0

  Integer euler_phi() const;        // Euler totient phi(n)
  Integer carmichael_lambda() const;  // Carmichael function lambda(n)

  // All positive divisors of n, sorted ascending.
  std::vector<Integer> divisors() const;

  // All positive divisors, each carrying its own (sub-)factorization. This is
  // the key feeder for the Erdos construction: a divisor d of L yields the
  // candidate prime d+1 together with the full factorization of (d+1)-1 = d,
  // at no extra cost. Materializes every divisor; for large n prefer the
  // streaming divisors_lazy below.
  std::vector<Factorization> divisors_factored() const;

  // Lazy (streaming) form of divisors_factored: yields each factored divisor in
  // turn, holding only one at a time -- O(number of distinct primes) memory
  // rather than O(number of divisors). Order is unspecified.
  std::generator<Factorization> divisors_lazy() const;

 private:
  std::vector<PrimePower> factors_;
};

// gcd of two integers whose factorizations are known. Computed by intersecting
// the prime supports and taking the minimum exponent of each shared prime --
// no Euclidean algorithm needed (Bach & Shallit, computing with factored
// integers). The result is itself a Factorization.
Factorization gcd(const Factorization& a, const Factorization& b);

}  // namespace pseudoprime
