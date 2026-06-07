#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "pseudoprime/factorization.hpp"

namespace pseudoprime {

// A descending chain of subgroups G = G_0 > G_1 > ... > G_t = {1} of the group
// G = (Z/ghat Z)^x, consumed by the large-Carmichael pairing engine. The engine never depends on
// how the chain is built -- only on these three operations -- so chain
// algorithms are interchangeable (strategy pattern).
class SubgroupChain {
 public:
  virtual ~SubgroupChain() = default;

  // Number of levels t (G_0/G_1, ..., G_{t-1}/G_t).
  virtual std::size_t levels() const = 0;

  // A key identifying a's coset in G_i / G_{i+1}. Two elements of G_i share a
  // key iff they lie in the same coset; elements already in G_{i+1} carry the
  // identity key (1) and are passed forward by the engine.
  virtual Integer coordinate(const Integer& a, std::size_t i) const = 0;

  // The key to pair `key` against at level i so that the product of the two
  // elements lands in G_{i+1}.
  virtual Integer inverse_key(const Integer& key, std::size_t i) const = 0;

  // The factor-group product of two coset keys at level i (the key of a*b in
  // G_i/G_{i+1}). Used by k-tuple pairing to find groups of elements whose
  // combined product lands in G_{i+1}.
  virtual Integer combine_keys(const Integer& a, const Integer& b, std::size_t i) const = 0;
};

// Factory for a SubgroupChain. The chain depends on ghat (the modulus of G),
// which is only known at run time once the construction has chosen it, so the
// chain is built lazily from ghat's factorization.
class SubgroupChainStrategy {
 public:
  virtual ~SubgroupChainStrategy() = default;
  virtual std::unique_ptr<SubgroupChain> build(const Factorization& ghat) const = 0;
};

// Congruence chain: G_i = { a : a == 1 (mod M_i) } for a chain of moduli
// 1 = M_0 | M_1 | ... | M_t = ghat. The coset key of a in G_i/G_{i+1} is a mod
// M_{i+1}; pairing uses the modular inverse. All residue arithmetic -- no
// discrete logs.
class CongruenceSubgroupChain : public SubgroupChain {
 public:
  // `moduli` is M_1, ..., M_t (ascending divisors of ghat, last == ghat). M_0 = 1.
  explicit CongruenceSubgroupChain(std::vector<Integer> moduli);

  std::size_t levels() const override { return moduli_.size(); }
  Integer coordinate(const Integer& a, std::size_t i) const override;
  Integer inverse_key(const Integer& key, std::size_t i) const override;
  Integer combine_keys(const Integer& a, const Integer& b, std::size_t i) const override;

 private:
  std::vector<Integer> moduli_;  // moduli_[i] == M_{i+1}
};

// Builds a CongruenceSubgroupChain by peeling ghat's prime-power structure one
// prime factor at a time (ascending by prime). Correct for arbitrary ghat; the
// factor groups are small (one prime each) but there are many levels, so the
// pairing pool depletes quickly.
class PrimePowerChainStrategy : public SubgroupChainStrategy {
 public:
  std::unique_ptr<SubgroupChain> build(const Factorization& ghat) const override;
};

// Builds a CongruenceSubgroupChain whose factor groups are approximately 2^l in
// size (l = round(sqrt(log2 |G|)), |G| = phi(ghat)), by accumulating ghat's
// prime-power factors into each modulus step until the step's index
// phi(M_{i+1})/phi(M_i) reaches ~2^l. This yields ~l levels instead of
// Omega(ghat), so the pairing pool depletes far less -- the balance point of the
// paper's Theorem 1.7. Still pure residue arithmetic (no discrete logs).
class BalancedCongruenceChainStrategy : public SubgroupChainStrategy {
 public:
  std::unique_ptr<SubgroupChain> build(const Factorization& ghat) const override;
};

}  // namespace pseudoprime
