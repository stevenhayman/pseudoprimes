#pragma once

#include <optional>

#include "pseudoprime/integer.hpp"

namespace pseudoprime {

// Outcome of a deterministic N-1 primality test.
//
// `Prime` and `Composite` are PROOFS (the N-1 tests are unconditional, with no
// probable-prime step). `Inconclusive` is an explicit "not proven either way"
// state: a test that searches witnesses up to a bound, or that is handed too
// small a factored part, may be unable to decide -- this is NOT evidence of
// compositeness. Callers must resolve an `Inconclusive` result by other means.
enum class Primality { Prime, Composite, Inconclusive };

// Result of a self-witnessing primality test. A `Composite` verdict may carry a
// nontrivial factor of n (1 < factor < n) discovered while testing witnesses.
// Finding one is not guaranteed, so `factor` is optional.
struct PrimalityCertificate {
  Primality verdict;
  std::optional<Integer> factor;
};

}  // namespace pseudoprime
