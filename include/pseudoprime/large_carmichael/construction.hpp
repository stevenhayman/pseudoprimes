#pragma once

#include <optional>

#include "pseudoprime/factorization.hpp"
#include "pseudoprime/large_carmichael/subgroup_chain.hpp"
#include "pseudoprime/primality.hpp"

namespace pseudoprime {

// Construct a Carmichael number with many prime factors from the factored
// modulus L, via the "Large Carmichael subset-product" method (Algorithm 2 of
// "Constructing Carmichael Numbers Through Improved Subset Product Algorithms"),
// in the generalized prime-power-chain variation.
//
// Builds the Erdos prime set P(L), forms the distinguished element a0 = b^{-1}
// (b = product of all of P), drives a0 to the identity of (Z/LZ)^x by Phase-1
// component cancellation plus Phase-2 level-by-level pairing through the
// injected subgroup chain, and returns the factorization of the complementary
// subset P \ T, whose product is congruent to 1 mod L (hence a Carmichael
// number). The result is verified with is_carmichael before being returned.
//
// Returns nullopt if the (heuristic) pairing cannot drive a0 to the identity or
// the resulting complement is not a Carmichael number. The subgroup-chain
// algorithm is a pluggable strategy; the default is the prime-power congruence
// chain.
// `max_pairing` is the largest group size the pairing combines at once: 2 pairs
// classes c with c^{-1} (the default); 3 also combines triples c1*c2*c3 == 1,
// which raises the success rate when inverse-pairs are scarce, at extra cost.
std::optional<Factorization> construct_large_carmichael(
    const Factorization& L, const SubgroupChainStrategy& chain_strategy = PrimePowerChainStrategy{},
    const PrimalityTest& primality = LayeredPrimalityTest{}, std::size_t max_pairing = 2);

}  // namespace pseudoprime
