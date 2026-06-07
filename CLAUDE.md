# Project guidelines for Claude

`pseudoprimes` is a C++23 library for constructing Carmichael numbers (ordinary
and higher-order). It is mathematical software: the code mirrors the algorithms
in Howe, *"Higher-order Carmichael numbers"*, and *"Constructing Carmichael
Numbers through Improved Subset-Product Algorithms"*. Keep edits consistent with
the conventions below — they are the house style this codebase already follows.

## Naming

- `snake_case` for locals, free functions, and variables; `PascalCase` for
  types; trailing-underscore for private members (e.g. `factors_`).
- No camelCase locals (do not write `Lval`/`ghatVal`; write `L_value` /
  `ghat_value`).
- Prefer full words over abbreviations when the name is a normal program
  identifier (`n_minus_one`, not `nm1`; `exponent`, not `exp`).
- **Paper-notation exception:** short single-symbol names that map directly to
  the source papers are *encouraged*, not flagged — `m`, `r`, `b`, `u`, `Q_u`,
  `N_j`, `ghat`, `a0`, and idiomatic loop/short locals (`pe`, `i`, `j`). Do not
  expand these; they keep the code readable alongside the papers.

## Comments

- Code comments are **ASCII-only**: write `phi`, `lambda`, `(Z/L Z)^x`, `ghat`,
  `==` for congruence, `>` for "is a supergroup of". Do not use Unicode math
  symbols in `.hpp`/`.cpp` files.
- Unicode math is fine in `README.md` prose (and only there).

## Core invariants — preserve these

- **Born factored.** A number whose factorization matters is passed as a
  `Factorization`, never re-factored. Thread the known factorization through APIs
  instead of recovering it.
- **Pluggable strategies.** New algorithms slot in behind the existing strategy
  interfaces: `PrimalityTest`, `SubsetProductSolver`, `PrimeSetStrategy`,
  `SubgroupChainStrategy`.
- **Lazy results.** Enumerations return `std::generator<...>` so callers pay only
  for what they consume; keep new enumerations lazy.
- **Deterministic primality only.** No probable-prime tests — a `true`/`false`
  primality answer must be a proof.

## Formatting

- 2-space indentation, ~100-column lines. There is no `.clang-format`; match the
  surrounding code.

## Build & test

```pwsh
cmake --preset mingw
cmake --build build/mingw
ctest --test-dir build/mingw --output-on-failure
```

Configure with `-DPP_BUILD_TESTS=OFF` to skip tests (and the GoogleTest fetch).
This repo develops on Windows with the MSYS2 MinGW64 toolchain (GMP from pacman);
`CMakePresets.json` also defines `linux-gcc` and `macos-clang` presets.
