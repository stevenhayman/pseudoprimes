#pragma once

#include <cstddef>
#include <optional>
#include <vector>

#include "pseudoprime/large_carmichael/subgroup_chain.hpp"
#include "pseudoprime/modular_product.hpp"

// The subgroup-chain pairing engine shared by the large-Carmichael construction
// (drive one distinguished element to the identity) and the layered
// subset-product solver (pair a pool down before brute-forcing the residual).
//
// IMPORTANT: the pairing is order-sensitive. `pair_pool` reverses each coordinate
// bucket so that pop-from-back consumes front-first, reproducing the original
// pairing order; changing the consumption order changes which numbers are found.
namespace pseudoprime::detail {

// Advance a single element through one level by combining it with up to
// (max_pairing - 1) pool partners whose coordinates multiply to the inverse of
// the element's coordinate, so the product lands in the next subgroup. Consumes
// the partners from `pool`. Returns nullopt if no such group exists.
std::optional<ModularProduct> advance(ModularProduct element, std::vector<ModularProduct>& pool,
                                      const SubgroupChain& chain, std::size_t level,
                                      std::size_t max_pairing);

// Pair a pool through a single level: bucket by coordinate (the coset in
// G_level / G_{level+1}); the identity class carries forward, and the rest are
// combined into groups whose coordinates multiply to the identity -- pairs
// (c with c^{-1}), and, when max_pairing >= 3, triples (c1*c2*c3 == 1). Unmatched
// leftovers are dropped (chaff). Returns the survivors.
std::vector<ModularProduct> pair_pool(const std::vector<ModularProduct>& pool,
                                      const SubgroupChain& chain, std::size_t level,
                                      std::size_t max_pairing);

}  // namespace pseudoprime::detail
