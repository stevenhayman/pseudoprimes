# pseudoprimes

A C++23 library for the construction of Carmichael numbers — ordinary ones and
**higher-order Carmichael numbers** — centered on the Erdős method and its generalizations.

## The Erdős construction

All Carmichael number construction in this library follows the Erdős method. Given a modulus `L`,
form the **prime set** `P(L) = { odd primes p : (p−1)|L, p∤L }`, then find subsets `S ⊆ P(L)`
with `|S| ≥ 3` whose product is `≡ 1 (mod L)`. By Korselt's criterion, each such product is a
Carmichael number. The algorithms differ only in how they solve this subset-product step; all
benefit from choosing `L` to be highly composite.

Every candidate prime `p = d + 1` comes from a divisor `d | L`, so `p − 1 = d` is already
factored — primality testing needs no additional factoring. The `Factorization` type is the
canonical handle for known factorizations; it threads through the entire API so nothing is
ever re-factored.

### Highly composite moduli

The prime set satisfies `|P(m,L)| ≤ d(L)` (the divisor count of `L`), so a richer `L` yields a
larger prime set and more Carmichael numbers. These generators stream highly composite `L` directly
as factorizations — born factored, ready to pass straight to `generate_prime_set` /
`construct_carmichael`:

- `superior_highly_composite(count)` — Ramanujan's **superior** highly composite numbers (OEIS
  A002201): `2, 6, 12, 60, 120, 360, 2520, …`. Each term differs from the previous by one prime
  factor; a clean "spine" for a growth-rate sweep.
- `candidate_highly_composite()` — every integer `2^a₁·3^a₂·…` with non-increasing exponents
  `a₁ ≥ a₂ ≥ … ≥ 1` (Ramanujan's necessary form), in increasing value order. A superset of the
  true HCNs; the flexible workhorse for exploration.
- `highly_composite()` — the true highly composite numbers (OEIS A002182):
  `1, 2, 4, 6, 12, 24, 36, 48, 60, …`. Obtained by filtering the candidate stream to strict
  records of `d(n)`, so no separate search is needed.

All three are lazy (`std::generator`). `factorization_from_exponents({a1, a2, …})` is the
low-level builder over the first few primes, and `Factorization::num_divisors()` returns
`d(n) = ∏(eᵢ + 1)` without enumerating divisors.

### Subset-product solvers

Three solvers implement the subset-product step, all behind the `SubsetProductSolver` strategy
interface and all honoring `SubsetProductProblem::max_solutions` (`0` = all, `1` = first,
`N` = stop after `N`):

- **`BruteForceSolver`** (`O(2ⁿ)`) — exhaustive search. Practical when `|P(L)|` is small.
- **`MeetInMiddleSolver`** (`O(2^(n/2))` time and space) — handles somewhat larger prime sets but
  becomes infeasible for `|P(L)|` in the thousands.
- **`LayeredSubsetProductSolver`** — heuristic: pairs elements of `P(L)` down a subgroup chain
  (see [below](#subgroup-chains)) until the pool is small enough to brute-force, then
  brute-forces the residual. The only solver implemented here that scales to `|P(L)|` in the hundreds of
  thousands. Finds *some* solutions rather than provably all.

#### The large-Carmichael algorithm

For very rich `L` where even enumerating solutions is infeasible, a specialized algorithm
(Alford, Grantham, Hayman, and Shallue, [arXiv:1203.6664](https://arxiv.org/abs/1203.6664))
targets a different goal: construct a *single* Carmichael number with many prime factors by
pairing elements of `P(L)` to the identity of `(ℤ/Lℤ)×` through a subgroup chain. The entry
point `construct_large_carmichael(L, chain, primality, max_pairing)` is a standalone path rather
than a pluggable `SubsetProductSolver`, since it produces one number rather than an enumeration.

This algorithm is heuristic — it succeeds for many `L` but not all (example: `L = 5040` →
an 18-prime Carmichael number). The `LayeredSubsetProductSolver` reuses the same pairing
engine, making it the natural choice when you want the `construct_carmichael` enumeration
interface over a huge `P(L)`: pair the pool down through a chain to reduce it, then brute-force
whatever remains.

##### Subgroup chains

Both `construct_large_carmichael` and `LayeredSubsetProductSolver` pair elements level by level
down a chain of subgroups of `(ℤ/Lℤ)×`. Two strategies are provided:

- **`PrimePowerChainStrategy`** — one prime-power factor group per level; default for
  `construct_large_carmichael`.
- **`BalancedCongruenceChainStrategy`** — factor groups of roughly equal size, fewer levels;
  default for `LayeredSubsetProductSolver`.

## Higher-order Carmichael numbers

Following Howe (*"Higher-order Carmichael numbers"*, [arXiv:math/9812089](https://arxiv.org/abs/math/9812089)), a **Carmichael number
of order `m`** is a composite squarefree `n` such that for every prime `p | n` and every `r` with
`1 ≤ r ≤ m`, there is an `i ≥ 0` with `n ≡ pⁱ (mod pʳ − 1)`. Order 1 is the ordinary notion. A
**rigid** Carmichael number of order `m` is the `i = 0` special case (`n ≡ 1 (mod pʳ − 1)` for
all such `p`, `r`).

- **Verification** — `is_carmichael(n, order = 1)` checks the general criterion (it admits
  non-rigid numbers); `is_rigid_carmichael(n, order = 1)` checks the rigid form. At order 1 the
  two coincide with Korselt. Both take a known factorization, so nothing is re-factored. The
  membership test `∃ i: n ≡ pⁱ` needs no discrete log — the reachable residues are exactly the
  `r` elements `{p⁰, …, p^{r−1}} (mod pʳ − 1)`.
- **Rigid construction** — the Erdős path generalizes directly: `generate_prime_set(L, …, order)`
  keeps the primes `p` with `(pʳ − 1) | L` for all `r ≤ m` (i.e. `P(m,L)`), and a subset
  product `≡ 1 (mod L)` over it yields a *rigid* Carmichael number of order `m` (Howe's
  `C(m,L)`). Feasible order-`m ≥ 2` examples require a richly composite `L` — the order-2 cases
  in the paper use `L` built from primes up to 31–37.

The §5 non-rigid construction and an order-`m` form of the large-Carmichael path are not
implemented.

## Primality testing

Every candidate prime is tested deterministically — no probable-prime tests. Two pluggable
`PrimalityTest` strategies are provided for use with `generate_prime_set` and
`construct_large_carmichael`:

- **`LayeredPrimalityTest`** (default) — base-2 Fermat (rejects most composites cheaply) →
  single-witness Lucas (proves most primes outright) → full Fellows–Koblitz fallback
  (`O((log n)⁶ / log log n)`, unconditional) for the rare remainder. Roughly 1000× faster per
  prime than running Fellows–Koblitz directly.
- **`FellowsKoblitzPrimalityTest`** — the full unconditional test directly, for when you want
  the highest rigor with no layering.

The underlying building blocks are also available as standalone functions:

- **`lucas_certifies_prime(n, n_minus_one)`** — proves primality via Lehmer's converse of
  Fermat; requires the complete factorization of n−1.
- **`bls_certifies_prime(n, f)`** — Brillhart–Lehmer–Selfridge N−1 test (BLS 1975): proves
  primality from a *partial* factorization — a fully-factored divisor `F` of n−1 with
  gcd(F, (n−1)/F) = 1. Succeeds as soon as F > n^(1/3) (Theorem 11) or F > √n (Pocklington
  branch), without needing the rest of n−1 factored. Returns `Prime`, `Composite`, or
  `Inconclusive` rather than a plain bool.

## Building

**Dependencies:** a C++23 compiler with `std::generator` (GCC 14+, MSVC 17.13+,
or LLVM 19+ / recent libc++), CMake ≥ 3.21, and GMP with `gmpxx`. GoogleTest is
used for tests and fetched automatically if not found locally.

This repo is developed on Windows with the MSYS2 MinGW64 toolchain (GMP from
pacman). `CMakePresets.json` also defines `linux-gcc` and `macos-clang` presets.

```pwsh
cmake --preset mingw
cmake --build build/mingw
```

GMP is located by `cmake/FindGMP.cmake` (GMP ships no CMake config). It uses
pkg-config as a hint and falls back to standard search paths; for a GMP in a
non-standard location, point it at the install prefix with `-DGMP_ROOT=...` (or
an environment variable of the same name).

## Using the library

The library has no command-line front end; you drive it from your own program by
linking the `pseudoprime` target. A modulus `L` is always supplied by its known
factorization, and the constructions yield results lazily, so you stop whenever
you like.

```cpp
#include <iostream>

#include "pseudoprime/carmichael.hpp"
#include "pseudoprime/erdos_construction.hpp"
#include "pseudoprime/factorization.hpp"
#include "pseudoprime/large_carmichael/construction.hpp"
#include "pseudoprime/large_carmichael/layered_solver.hpp"
#include "pseudoprime/subset_product.hpp"

using namespace pseudoprime;

int main() {
  // L by its factorization. Here L = 2^3 * 3 * 5 = 120; use a richly composite L
  // (e.g. 2^12*3^6*5^3*7^2*11^2*13*17*19*23*29*31*37*41*43*47) for the layered
  // and large paths below.
  const Factorization L({{Integer(2), 3}, {Integer(3), 1}, {Integer(5), 1}});

  // (1) Erdos construction with BruteForceSolver — enumerate every Carmichael
  //     number whose prime factors are a subset of P(L) with product == 1 (mod L).
  //     Solvers are pluggable: BruteForceSolver / MeetInMiddleSolver (exhaustive)
  //     or LayeredSubsetProductSolver (heuristic, scales to a huge P(L)).
  BruteForceSolver brute;
  for (const ConstructedCarmichael& cc : construct_carmichael(L, brute)) {
    std::cout << cc.value << "  (" << cc.factorization.num_distinct_primes() << " primes)\n";
  }

  // (2) Erdos construction with LayeredSubsetProductSolver — pairs P(L) down a
  //     subgroup chain, then brute-forces the residual. The trailing args are
  //     min_factors, order, max_solutions (0 = all, 1 = first, N = stop after N).
  LayeredSubsetProductSolver layered(L);  // balanced chain by default
  for (const ConstructedCarmichael& cc :
       construct_carmichael(L, layered, /*min_factors=*/3, /*order=*/1, /*max_solutions=*/1)) {
    std::cout << "layered: " << cc.value << "\n";
  }

  // (3) Erdos construction, higher-order (rigid) — pass order = m. Verify with
  //     is_carmichael(n, m) (general) or is_rigid_carmichael(n, m).
  for (const ConstructedCarmichael& cc : construct_carmichael(L, brute, 3, /*order=*/2)) {
    std::cout << "order-2: " << cc.value
              << (is_rigid_carmichael(cc.factorization, 2) ? "  [rigid]\n" : "  [BAD]\n");
  }

  // (4) Large-Carmichael algorithm — one number with many prime factors (arXiv:1203.6664).
  if (std::optional<Factorization> n = construct_large_carmichael(L)) {
    std::cout << "large: " << n->value() << "  (" << n->num_distinct_primes() << " primes)\n";
  }
  return 0;
}
```

Build it by adding this repository to your own CMake project and linking the
`pseudoprime` target (it carries its include paths and GMP linkage transitively):

```cmake
cmake_minimum_required(VERSION 3.21)
project(my_driver CXX)
add_subdirectory(path/to/pseudoprimes pseudoprimes-build)
add_executable(driver driver.cpp)
target_compile_features(driver PRIVATE cxx_std_23)
target_link_libraries(driver PRIVATE pseudoprime)
```

Or, once the library is built, compile directly against the static lib:

```pwsh
g++ -std=c++23 -Iinclude driver.cpp build/mingw/libpseudoprime.a -lgmpxx -lgmp -o driver
```

Entry points, by capability:

| Capability | API |
| --- | --- |
| Enumerate Carmichael numbers from `L` | `construct_carmichael(L, solver, min_factors, order, max_solutions)` |
| Subset-product solvers | `BruteForceSolver`, `MeetInMiddleSolver`, `LayeredSubsetProductSolver` |
| Large-Carmichael algorithm (one number, many prime factors) | `construct_large_carmichael(L, chain, primality, max_pairing)` |
| Subgroup chain strategies | `PrimePowerChainStrategy`, `BalancedCongruenceChainStrategy` |
| Verify (order `m`, general / rigid) | `is_carmichael(n, m)`, `is_rigid_carmichael(n, m)` |
| Build `P(m,L)` directly | `generate_prime_set(L, primality, order)` |
| Generate rich `L` (SHCN / HCN) | `superior_highly_composite()`, `candidate_highly_composite()`, `highly_composite()` |
| Primality strategies | `LayeredPrimalityTest`, `FellowsKoblitzPrimalityTest` |
| Primality building blocks | `lucas_certifies_prime()`, `bls_certifies_prime()` |

## Running the tests

Tests use GoogleTest and run under CTest:

```pwsh
ctest --test-dir build/mingw --output-on-failure
```

To skip the tests entirely (and the GoogleTest fetch), configure with
`-DPP_BUILD_TESTS=OFF`.
