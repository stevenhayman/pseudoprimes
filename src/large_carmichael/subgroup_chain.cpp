#include "pseudoprime/large_carmichael/subgroup_chain.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace pseudoprime {

CongruenceSubgroupChain::CongruenceSubgroupChain(std::vector<Integer> moduli)
    : moduli_(std::move(moduli)) {}

Integer CongruenceSubgroupChain::coordinate(const Integer& a, std::size_t i) const {
  Integer c = a % moduli_[i];
  if (c < 0) c += moduli_[i];
  return c;
}

Integer CongruenceSubgroupChain::inverse_key(const Integer& key, std::size_t i) const {
  const std::optional<Integer> inv = mod_inverse(key, moduli_[i]);
  if (!inv) {
    // Coset keys are units in G_i/G_{i+1}, so this is unreachable for valid input.
    throw std::invalid_argument("CongruenceSubgroupChain: key is not invertible modulo M_{i+1}");
  }
  return *inv;
}

Integer CongruenceSubgroupChain::combine_keys(const Integer& a, const Integer& b,
                                              std::size_t i) const {
  Integer c = (a * b) % moduli_[i];
  if (c < 0) c += moduli_[i];
  return c;
}

std::unique_ptr<SubgroupChain> PrimePowerChainStrategy::build(const Factorization& ghat) const {
  // Peel one prime factor at a time, ascending by prime: M runs 1 -> ... -> ghat.
  std::vector<Integer> moduli;
  Integer m = 1;
  for (const PrimePower& pe : ghat.factors()) {
    for (unsigned long k = 0; k < pe.exponent; ++k) {
      m *= pe.prime;
      moduli.push_back(m);
    }
  }
  return std::make_unique<CongruenceSubgroupChain>(std::move(moduli));
}

std::unique_ptr<SubgroupChain> BalancedCongruenceChainStrategy::build(
    const Factorization& ghat) const {
  // l = round(sqrt(log2 |G|)); aim for each factor group phi(M_{i+1})/phi(M_i) ~ 2^l.
  signed long exp2 = 0;
  const double mant = mpz_get_d_2exp(&exp2, ghat.euler_phi().get_mpz_t());  // |G| = mant*2^exp2
  const double log2_order = std::log2(mant) + static_cast<double>(exp2);
  const auto ell = std::max<long>(1, std::lround(std::sqrt(std::max(0.0, log2_order))));
  const double target = std::pow(2.0, static_cast<double>(ell));

  std::vector<Integer> moduli;
  Integer m = 1;
  double step_index = 1.0;  // cumulative phi-ratio since the last emitted modulus
  for (const PrimePower& pe : ghat.factors()) {
    const double q = mpz_get_d(pe.prime.get_mpz_t());
    for (unsigned long k = 0; k < pe.exponent; ++k) {
      // Index contribution of multiplying m by q: (q-1) the first time q enters
      // (a new prime, k==0), else q (raising an existing prime's power).
      step_index *= (k == 0) ? (q - 1.0) : q;
      m *= pe.prime;
      if (step_index >= target) {
        moduli.push_back(m);
        step_index = 1.0;
      }
    }
  }
  // Fold any remainder into the final level so the chain reaches ghat exactly.
  if (moduli.empty() || moduli.back() != ghat.value()) moduli.push_back(ghat.value());

  return std::make_unique<CongruenceSubgroupChain>(std::move(moduli));
}

}  // namespace pseudoprime
