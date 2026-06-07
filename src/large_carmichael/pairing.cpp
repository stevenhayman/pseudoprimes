#include "pseudoprime/large_carmichael/pairing.hpp"

#include <algorithm>
#include <cstddef>
#include <map>
#include <optional>
#include <vector>

namespace pseudoprime::detail {

std::optional<ModularProduct> advance(ModularProduct element, std::vector<ModularProduct>& pool,
                                      const SubgroupChain& chain, std::size_t level,
                                      std::size_t max_pairing) {
  const Integer c = chain.coordinate(element.residue(), level);
  if (c == 1) return element;  // already in the next subgroup
  const Integer need = chain.inverse_key(c, level);  // partner coordinates must multiply to this

  // Find a group of partners whose coordinates multiply to `need`, trying smaller
  // groups first (1 partner = a pair with `element`, 2 = a triple, ...). For each
  // size we take the lexicographically-smallest increasing-index combination,
  // generalizing the former pair/triple branches to any max_pairing. Pop order is
  // left-to-right; combine order is irrelevant to the result (residues commute and
  // the history is consumed as a set).
  const std::size_t n = pool.size();
  std::vector<Integer> coords(n);
  for (std::size_t k = 0; k < n; ++k) coords[k] = chain.coordinate(pool[k].residue(), level);

  std::vector<std::size_t> chosen;
  // DFS for `remaining` more partners (increasing index from `start`) whose
  // coordinates, times `acc`, equal `need`. Returns true on the first success,
  // which is the lexicographically-smallest combination.
  const auto search = [&](this auto&& self, std::size_t start, std::size_t remaining,
                          const Integer& acc) -> bool {
    if (remaining == 0) return acc == need;
    for (std::size_t i = start; i < n; ++i) {
      chosen.push_back(i);
      if (self(i + 1, remaining - 1, chain.combine_keys(acc, coords[i], level))) return true;
      chosen.pop_back();
    }
    return false;
  };

  for (std::size_t partners = 1; partners < max_pairing; ++partners) {
    chosen.clear();
    if (!search(0, partners, Integer(1))) continue;
    ModularProduct combined = element;
    for (std::size_t idx : chosen) combined = combined.combine(pool[idx]);
    // Erase consumed partners largest-index-first so earlier indices stay valid.
    for (std::size_t t = chosen.size(); t-- > 0;)
      pool.erase(pool.begin() + static_cast<long>(chosen[t]));
    return combined;
  }
  return std::nullopt;
}

std::vector<ModularProduct> pair_pool(const std::vector<ModularProduct>& pool,
                                      const SubgroupChain& chain, std::size_t level,
                                      std::size_t max_pairing) {
  const Integer one(1);
  std::map<Integer, std::vector<std::size_t>> bucket;
  for (std::size_t k = 0; k < pool.size(); ++k)
    bucket[chain.coordinate(pool[k].residue(), level)].push_back(k);

  std::vector<ModularProduct> survivors;
  if (const auto id = bucket.find(one); id != bucket.end()) {
    for (std::size_t idx : id->second) survivors.push_back(pool[idx]);
    id->second.clear();
  }

  // Consume from the front (reverse so the pop-from-back helpers do so), which
  // matches the original pairing order -- the heuristic is order-sensitive.
  for (auto& kv : bucket) std::reverse(kv.second.begin(), kv.second.end());

  const auto count = [&](const Integer& c) -> std::size_t {
    const auto it = bucket.find(c);
    return it == bucket.end() ? 0 : it->second.size();
  };
  const auto pop = [&](const Integer& c) -> std::size_t {
    auto& v = bucket.find(c)->second;
    const std::size_t idx = v.back();
    v.pop_back();
    return idx;
  };

  std::vector<Integer> keys;
  keys.reserve(bucket.size());
  for (const auto& kv : bucket) keys.push_back(kv.first);

  // Combine groups of k elements whose coordinates multiply to the identity, for
  // each group size k from 2 (pairs) up to `max_pairing`. A group is a canonical
  // NON-DECREASING coordinate tuple (c_1 <= ... <= c_k) of non-identity classes:
  // the first k-1 coordinates are chosen freely (non-decreasing, so each tuple is
  // enumerated once) and the last is forced to be the inverse of their product.
  // This is the common generalization of the pair (k=2) and triple (k=3) cases.
  std::vector<Integer> tuple;  // the coordinate tuple currently being assembled

  // Consume copies of `tuple` while the buckets hold enough: each distinct
  // coordinate needs at least its multiplicity in `tuple`. Pop left to right and
  // fold-combine, matching the original (order-sensitive) pairing order.
  const auto consume = [&] {
    const auto available = [&] {
      for (const Integer& c : tuple)
        if (count(c) < static_cast<std::size_t>(std::ranges::count(tuple, c))) return false;
      return true;
    };
    while (available()) {
      std::optional<ModularProduct> acc;
      for (const Integer& c : tuple) {
        ModularProduct e = pool[pop(c)];
        acc = acc ? acc->combine(e) : std::move(e);
      }
      survivors.push_back(std::move(*acc));
    }
  };

  // Pick the free coordinates c_1 <= ... <= c_{group_size-1} (recursing over
  // `keys` from `start` to stay non-decreasing), then force the last coordinate.
  const auto combine_groups = [&](this auto&& self, std::size_t group_size,
                                  std::size_t start) -> void {
    if (tuple.size() + 1 == group_size) {
      Integer prod = tuple.front();
      for (std::size_t t = 1; t < tuple.size(); ++t) prod = chain.combine_keys(prod, tuple[t], level);
      const Integer last = chain.inverse_key(prod, level);
      if (last != one && last >= tuple.back()) {  // canonical (non-decreasing) and non-identity
        tuple.push_back(last);
        consume();
        tuple.pop_back();
      }
      return;
    }
    for (std::size_t ki = start; ki < keys.size(); ++ki) {
      if (keys[ki] == one) continue;
      tuple.push_back(keys[ki]);
      self(group_size, ki);  // ki (not ki+1): repeats allowed, tuple stays non-decreasing
      tuple.pop_back();
    }
  };

  for (std::size_t k = 2; k <= max_pairing; ++k) combine_groups(k, 0);
  return survivors;
}

}  // namespace pseudoprime::detail
