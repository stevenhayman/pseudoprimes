# Project guidelines for coding agents

`pseudoprimes` is a C++23 library for constructing Carmichael numbers, ordinary and
higher-order, via the Erdős method.

## Relationship to the papers

Two papers supply the mathematics and the notation:

- Howe, *"Higher-order Carmichael numbers"* (arXiv:math/9812089) — the order-`m` and rigid
  definitions, and the `C(m,L)` construction.
- Alford, Grantham, Hayman, and Shallue, *"Constructing Carmichael Numbers through Improved
  Subset-Product Algorithms"* (arXiv:1203.6664) — the large-Carmichael pairing algorithm.

The library follows these papers for definitions, notation, and the two algorithms named above.
It is not a transcription of them:

- Data structures, solver design, and search order are this codebase's own.
- Coverage is partial. Howe's Section 5 non-rigid construction and an order-`m` form of the
  large-Carmichael path are not implemented.
- Several components come from neither paper: `BruteForceSolver`, `MeetInMiddleSolver`, the
  highly-composite modulus generators, and `LayeredSubsetProductSolver`, which reuses the
  large-Carmichael pairing engine for a goal the paper does not address.

Do not assume a function reproduces a published algorithm step for step. The header comment is
the authority on what a component actually guarantees — in particular which solvers are
exhaustive (`BruteForceSolver`, `MeetInMiddleSolver`) and which are heuristic
(`LayeredSubsetProductSolver`, `construct_large_carmichael`). Read it before changing behavior,
and keep it accurate when you do.

## Naming

- `snake_case` for locals, free functions, and variables; `PascalCase` for types;
  trailing underscore for private members (e.g. `factors_`).
- No camelCase locals: write `L_value` / `ghat_value`, not `Lval` / `ghatVal`.
- Prefer full words for ordinary program identifiers (`n_minus_one`, not `nm1`; `exponent`,
  not `exp`).
- **Paper notation wins.** Identifiers carrying the papers' symbols keep their shape, capitals
  and subscript suffixes included: `m`, `r`, `b`, `u`, `ghat`, `a0`, `L_value`, `G_i`, `U_k`,
  `Q_u`, `N_j`. Idiomatic short locals (`pe`, `i`, `j`) are fine. Do not expand or "correct"
  these; they keep the code readable alongside the papers.

## Comments

Comments in `.hpp`/`.cpp` are ASCII-only: write `phi`, `lambda`, `(Z/L Z)^x`, `ghat`, `==` for
congruence, `>` for "is a supergroup of". Unicode math belongs in the Markdown files
(`README.md`, `apps/*.md`), never in source.

## Core invariants — preserve these

- **Born factored.** A number whose factorization matters is passed as a `Factorization`, never
  re-factored. Thread the known factorization through APIs instead of recovering it.
- **Pluggable strategies.** New algorithms slot in behind the existing strategy interfaces:
  `PrimalityTest`, `SubsetProductSolver`, `PrimeSetStrategy`, `SubgroupChainStrategy`.
- **Lazy results.** Enumerations return `std::generator<...>` so callers pay only for what they
  consume; keep new enumerations lazy.
- **Deterministic primality only.** No probable-prime tests — a `true`/`false` primality answer
  must be a proof. The pseudoprimality checkers in `bpsw.hpp` and `lucas_sequence.hpp` verify
  congruences of proven composites (like `is_carmichael`); they live outside `primality/` on
  purpose, are never `PrimalityTest`s, and nothing in `primality/` may call them.

## Formatting

Indentation, line width, and whitespace are in `.editorconfig`. There is no `.clang-format`;
for anything `.editorconfig` does not cover, match the surrounding code.

## Build & test

```pwsh
cmake --preset mingw
cmake --build build/mingw
ctest --test-dir build/mingw --output-on-failure
```

Configure with `-DPP_BUILD_TESTS=OFF` to skip tests (and the GoogleTest fetch). This repo
develops on Windows with the MSYS2 MinGW64 toolchain (GMP from pacman); `CMakePresets.json`
also defines `linux-gcc` and `macos-clang` presets.
