#include "pseudoprime/prime_set_strategy.hpp"

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

#include "detail.hpp"

namespace pseudoprime {

namespace {

// A prime dividing L, as a machine word, with its exponent in L. The exponent
// doubles as a cap during sieving: a candidate cannot contain this prime more
// times than L does and still divide L.
struct PrimeFactor {
  std::uint64_t prime;
  unsigned long exponent;
};

// Does the non-negative Integer value fit in a 64-bit unsigned word?
bool fits_in_uint64(const Integer& value) {
  return mpz_sgn(value.get_mpz_t()) >= 0 && mpz_sizeinbase(value.get_mpz_t(), 2) <= 64;
}

// Convert a non-negative Integer known to fit in 64 bits to std::uint64_t.
// (mpz_get_ui yields only an unsigned long, which is 32-bit on LLP64 Windows, so
// assemble the value from two limbs explicitly.)
std::uint64_t to_uint64(const Integer& value) {
  Integer remaining = value;
  const std::uint64_t low_bits = mpz_get_ui(remaining.get_mpz_t()) & 0xffffffffULL;
  remaining >>= 32;
  const std::uint64_t high_bits = mpz_get_ui(remaining.get_mpz_t()) & 0xffffffffULL;
  return (high_bits << 32) | low_bits;
}

// Build an Integer from a 64-bit unsigned word (unsigned long is 32-bit on LLP64
// Windows, so assemble from two halves).
Integer to_integer(std::uint64_t value) {
  return (Integer(static_cast<unsigned long>(value >> 32)) << 32) +
         Integer(static_cast<unsigned long>(value & 0xffffffffULL));
}

// Re-derive the factorization of a divisor of L using the prime factors of L. The
// divisor is known to be smooth over those primes within their exponents (the
// sieve guaranteed it), so the result is exact and complete and matches what
// divisors_lazy would have carried.
Factorization factorize_over(std::uint64_t divisor, const std::vector<PrimeFactor>& prime_factors) {
  std::vector<PrimePower> factors;
  for (const PrimeFactor& factor : prime_factors) {
    unsigned long exponent = 0;
    while (divisor % factor.prime == 0) {
      divisor /= factor.prime;
      ++exponent;
    }
    if (exponent > 0)
      factors.push_back({Integer(static_cast<unsigned long>(factor.prime)), exponent});
  }
  return Factorization(std::move(factors));  // ctor sorts; prime_factors is already ascending
}

// The prime factors of L (the primes dividing L) as machine words with their
// exponents. The primes of a highly composite L are tiny; require they fit in 64
// bits.
std::vector<PrimeFactor> prime_factors_of(const Factorization& L) {
  std::vector<PrimeFactor> prime_factors;
  prime_factors.reserve(L.factors().size());
  for (const PrimePower& factor : L.factors()) {
    if (!fits_in_uint64(factor.prime))
      throw std::invalid_argument("SievePrimeSetStrategy: prime factor of L exceeds 64 bits");
    prime_factors.push_back({to_uint64(factor.prime), factor.exponent});
  }
  return prime_factors;
}

// Resolve the inclusive sweep cap on the candidate value d = p - 1. cap == 0 means
// "complete": cap at the exact per-order ceiling ceil((L+1)^(1/order)), beyond
// which no order-m member can exist. The sieve runs in 64-bit words, so the
// effective cap must fit in 64 bits (for order 1 the complete ceiling is L + 1).
std::uint64_t resolve_cap(std::uint64_t cap, const Factorization& L, unsigned long order) {
  if (cap != 0) return cap;
  const Integer ceiling = provable_cap(L, order);
  if (!fits_in_uint64(ceiling))
    throw std::invalid_argument(
        "SievePrimeSetStrategy: ceil((L+1)^(1/order)) exceeds 64 bits; an explicit cap is "
        "required");
  return to_uint64(ceiling);
}

// Factorizing sieve over one window: set residual[i] to the part of d = window_low
// + i left after dividing out each prime factor of L up to its exponent. A slot
// that reaches 1 is exactly a divisor of L; a slot greater than 1 means d has a
// prime outside the factors of L, or one of them appeared more times than its
// exponent (a leftover power remains). Fills residual[0 .. window_high-window_low];
// residual must hold that many slots.
void sieve_window(std::vector<std::uint64_t>& residual, std::uint64_t window_low,
                  std::uint64_t window_high, const std::vector<PrimeFactor>& prime_factors) {
  const std::size_t slot_count = static_cast<std::size_t>(window_high - window_low + 1);
  for (std::size_t slot_index = 0; slot_index < slot_count; ++slot_index)
    residual[slot_index] = window_low + static_cast<std::uint64_t>(slot_index);

  for (const PrimeFactor& factor : prime_factors) {
    const std::uint64_t prime = factor.prime;
    const std::uint64_t first_multiple = ((window_low + prime - 1) / prime) * prime;
    for (std::uint64_t multiple = first_multiple; multiple >= window_low && multiple <= window_high;
         multiple += prime) {
      std::uint64_t& slot = residual[static_cast<std::size_t>(multiple - window_low)];
      unsigned long removed = 0;
      while (removed < factor.exponent && slot % prime == 0) {
        slot /= prime;
        ++removed;
      }
    }
  }
}

// Apply the exact order-m membership checks to a divisor of L and return the member
// prime p = divisor + 1, or nullopt if it fails one. Mirrors the per-candidate
// checks of DivisorPrimeSetStrategy: p odd, p does not divide L, p prime (certified
// from the free factorization of divisor = p - 1), and (p^r - 1) | L for
// r = 2..order.
std::optional<Integer> member_from_divisor(std::uint64_t divisor, const Integer& L_value,
                                           const std::vector<PrimeFactor>& prime_factors,
                                           const PrimalityTest& primality, unsigned long order) {
  if ((divisor & 1u) != 0) return std::nullopt;  // Carmichael factors are odd: p odd, so divisor even
  const Integer p = to_integer(divisor) + 1;
  if (L_value % p == 0) return std::nullopt;  // p must not divide L
  const Factorization divisor_factorization = factorize_over(divisor, prime_factors);
  if (!primality.is_prime(p, divisor_factorization)) return std::nullopt;
  if (order > 1 && !detail::satisfies_higher_order(p, L_value, order)) return std::nullopt;
  return p;
}

}  // namespace

Integer provable_cap(const Factorization& L, unsigned long order) {
  if (order == 0) order = 1;
  const Integer base = L.value() + 1;  // (L + 1)
  Integer root;
  // mpz_root sets root = floor(base^(1/order)) and returns non-zero iff exact.
  const int exact = mpz_root(root.get_mpz_t(), base.get_mpz_t(), order);
  if (exact == 0) root += 1;  // round up to the ceiling
  return root;                // ceil((L+1)^(1/order))
}

std::generator<Integer> SievePrimeSetStrategy::generate(
    const Factorization& L, const PrimalityTest& primality, unsigned long order) const {
  const Integer L_value = L.value();
  const std::vector<PrimeFactor> prime_factors = prime_factors_of(L);
  const std::uint64_t cap = resolve_cap(cap_, L, order);
  if (cap == 0) co_return;  // L == 1: no candidates

  const std::uint64_t window_width = static_cast<std::uint64_t>(window_);
  std::vector<std::uint64_t> residual(static_cast<std::size_t>(window_width));

  // Segmented sieve over d in [1, cap]: sieve each window to find the divisors of L
  // it contains, then confirm each as an order-m member exactly before yielding.
  for (std::uint64_t window_low = 1; window_low <= cap; window_low += window_width) {
    const std::uint64_t window_high =
        (cap - window_low < window_width - 1) ? cap : window_low + window_width - 1;  // inclusive
    sieve_window(residual, window_low, window_high, prime_factors);

    const std::size_t slot_count = static_cast<std::size_t>(window_high - window_low + 1);
    for (std::size_t slot_index = 0; slot_index < slot_count; ++slot_index) {
      if (residual[slot_index] != 1) continue;  // window_low + slot_index does not divide L
      const std::uint64_t divisor = window_low + static_cast<std::uint64_t>(slot_index);
      if (std::optional<Integer> p = member_from_divisor(divisor, L_value, prime_factors, primality, order))
        co_yield std::move(*p);
    }
  }
}

}  // namespace pseudoprime
