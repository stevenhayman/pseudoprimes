#include "pseudoprime/large_carmichael/construction.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

#include "pseudoprime/carmichael.hpp"
#include "pseudoprime/erdos_construction.hpp"
#include "pseudoprime/large_carmichael/pairing.hpp"
#include "pseudoprime/large_carmichael/weight.hpp"
#include "pseudoprime/modular_product.hpp"

namespace pseudoprime {

namespace {

// The pairing primitives now live in the shared pairing engine.
using detail::advance;
using detail::pair_pool;

// Natural log of a positive integer, via its mantissa/exponent (no overflow).
double ln_integer(const Integer& x) {
  signed long exp2 = 0;
  const double mant = mpz_get_d_2exp(&exp2, x.get_mpz_t());  // x = mant * 2^exp2
  return std::log(mant) + static_cast<double>(exp2) * std::log(2.0);
}

Integer norm(Integer v, const Integer& mod) {
  v %= mod;
  if (v < 0) v += mod;
  return v;
}

// The Erdos prime set P(L) with, for each prime, its residue mod L and its
// weight omega (Definition 1.2).
struct WeightedPrimeSet {
  std::vector<Integer> primes;
  std::vector<Integer> residues;  // residues[i] = primes[i] mod L
  std::vector<std::size_t> weights;  // weights[i] = weight(primes[i]); the paper's omega, tabulated
  std::size_t size() const { return primes.size(); }
};

WeightedPrimeSet build_weighted_prime_set(const Factorization& L,
                                          const std::vector<Integer>& components,
                                          const PrimalityTest& primality) {
  WeightedPrimeSet s;
  const Integer L_value = L.value();
  const std::vector<Integer> P = generate_prime_set(L, primality);
  s.primes.reserve(P.size());
  s.residues.reserve(P.size());
  s.weights.reserve(P.size());
  for (const Integer& p : P) {
    s.primes.push_back(p);
    s.residues.push_back(p % L_value);  // p mod L
    s.weights.push_back(weight(p, components));
  }
  return s;
}

// Choose m: the smallest j with N_j > (ln L) * 4^{sqrt(ln(Q_1...Q_j))}, where
// N_j = #{ primes with weight <= j }; falls back to r if no j qualifies.
std::size_t choose_modulus_count(const WeightedPrimeSet& s,
                                 const std::vector<Integer>& components) {
  const std::size_t r = components.size();
  std::vector<double> cumln(r + 1, 0.0);
  for (std::size_t j = 1; j <= r; ++j) cumln[j] = cumln[j - 1] + ln_integer(components[j - 1]);
  const double lnL = cumln[r];

  for (std::size_t j = 1; j <= r; ++j) {
    std::size_t Nj = 0;
    for (std::size_t w : s.weights)
      if (w <= j) ++Nj;
    const double threshold = lnL * std::pow(4.0, std::sqrt(cumln[j]));
    if (static_cast<double>(Nj) > threshold) return j;
  }
  return r;
}

// The distinguished element a0 = b^{-1} (b = product of all of P) after Phase-1
// reduction, together with the prime indices (T) folded into it.
struct DistinguishedElement {
  Integer residue;  // mod L
  std::vector<std::size_t> elements;
};

// Phase 1: form a0 = b^{-1} and cancel its components u = r .. m+1 by multiplying
// in matching primes, lowering omega(a0) to <= m. Marks the consumed primes in
// `used`. Returns nullopt if a needed matching prime is missing.
std::optional<DistinguishedElement> build_distinguished_element(
    const WeightedPrimeSet& s, const std::vector<Integer>& components, std::size_t m,
    const Integer& L_value, std::vector<bool>& used) {
  const std::size_t n = s.size();
  const std::size_t r = components.size();

  Integer b = 1;
  for (const Integer& rr : s.residues) b = (b * rr) % L_value;
  const std::optional<Integer> b_inverse = mod_inverse(b, L_value);
  if (!b_inverse) return std::nullopt;

  DistinguishedElement d;
  d.residue = *b_inverse;
  for (std::size_t u = r; u > m; --u) {
    if (weight(d.residue, components) != u) continue;  // top nontrivial component already below u
    const Integer& Qu = components[u - 1];

    // Need a prime p with p == a0^{-1} (mod Q_u).
    const std::optional<Integer> target = mod_inverse(norm(d.residue, Qu), Qu);
    if (!target) return std::nullopt;

    std::size_t found = n;
    for (std::size_t i = 0; i < n; ++i) {
      if (used[i] || s.weights[i] != u) continue;
      if (norm(s.residues[i], Qu) == *target) {
        found = i;
        break;
      }
    }
    if (found == n) return std::nullopt;
    d.residue = (d.residue * s.residues[found]) % L_value;
    used[found] = true;
    d.elements.push_back(found);
  }
  return d;
}

// The factorization of ghat = Q_1 * ... * Q_m (the first m prime-power components).
Factorization leading_factorization(const Factorization& L, std::size_t m) {
  std::vector<PrimePower> factors(L.factors().begin(), L.factors().begin() + static_cast<long>(m));
  return Factorization(std::move(factors));
}

// The Phase-2 pool P_m, as products mod ghat (the distinguished a0 is tracked
// separately by the caller, not placed in the pool).
std::vector<ModularProduct> build_pool(const WeightedPrimeSet& s, const std::vector<bool>& used,
                                       std::size_t m, const Integer& ghat_value) {
  std::vector<ModularProduct> pool;
  for (std::size_t i = 0; i < s.size(); ++i) {
    if (used[i] || s.weights[i] > m) continue;  // P_m
    pool.emplace_back(ModularInteger(norm(s.residues[i], ghat_value), ghat_value),
                      std::vector<std::size_t>{i});
  }
  return pool;
}

// Phase 2: drive the distinguished element to the identity of G by, at each
// level, advancing it (combined with pool partners) and then pairing the pool.
// `max_pairing` is the largest group size used when combining (2 = pairs only).
// Returns the identity product, or nullopt if it cannot be reached.
std::optional<ModularProduct> pair_to_identity(std::vector<ModularProduct> pool,
                                               ModularProduct distinguished,
                                               const SubgroupChain& chain,
                                               std::size_t max_pairing) {
  for (std::size_t level = 0; level < chain.levels(); ++level) {
    std::optional<ModularProduct> advanced =
        advance(std::move(distinguished), pool, chain, level, max_pairing);
    if (!advanced) return std::nullopt;
    distinguished = std::move(*advanced);
    pool = pair_pool(pool, chain, level, max_pairing);
  }
  if (distinguished.residue() != 1) return std::nullopt;
  return distinguished;
}

// The Carmichael number is the complement P \ elements (its product is 1 mod L).
// Verified with Korselt's criterion before being returned.
std::optional<Factorization> carmichael_from_complement(const WeightedPrimeSet& s,
                                                        const std::vector<std::size_t>& elements) {
  std::vector<bool> chosen(s.size(), false);
  for (std::size_t idx : elements) chosen[idx] = true;

  std::vector<Integer> factors;
  for (std::size_t i = 0; i < s.size(); ++i)
    if (!chosen[i]) factors.push_back(s.primes[i]);

  const Factorization result = Factorization::from_primes(factors);
  if (is_carmichael(result)) return result;
  return std::nullopt;
}

}  // namespace

std::optional<Factorization> construct_large_carmichael(
    const Factorization& L, const SubgroupChainStrategy& chain_strategy,
    const PrimalityTest& primality, std::size_t max_pairing) {
  const Integer L_value = L.value();
  const std::vector<Integer> components = L.prime_powers();  // Q_1..Q_r
  if (components.empty()) return std::nullopt;

  const WeightedPrimeSet s = build_weighted_prime_set(L, components, primality);
  if (s.size() == 0) return std::nullopt;

  const std::size_t m = choose_modulus_count(s, components);

  // Phase 1: build the distinguished element a0 (and mark the primes it absorbed).
  std::vector<bool> used(s.size(), false);
  const std::optional<DistinguishedElement> a0 =
      build_distinguished_element(s, components, m, L_value, used);
  if (!a0) return std::nullopt;

  // Phase 2: pair the pool to the identity of G = (Z/ghat Z)^x, with a0 tracked
  // separately and advanced first at each level.
  const Factorization ghat = leading_factorization(L, m);
  const Integer ghat_value = ghat.value();
  const std::unique_ptr<SubgroupChain> chain = chain_strategy.build(ghat);

  ModularProduct a0_product(ModularInteger(norm(a0->residue, ghat_value), ghat_value), a0->elements);
  std::vector<ModularProduct> pool = build_pool(s, used, m, ghat_value);

  const std::optional<ModularProduct> identity =
      pair_to_identity(std::move(pool), std::move(a0_product), *chain, max_pairing);
  if (!identity) return std::nullopt;

  return carmichael_from_complement(s, identity->elements());
}

}  // namespace pseudoprime
