# Agent Notes for PALP

This file contains project-specific context and conventions for agents working
on the PALP codebase.

## Project overview

PALP (a Package for Analyzing Lattice Polytopes) is a command-line package for
toric geometry and reflexive polytope classification. This repository is a fork
undergoing a C → C++20 modernization. The goal is idiomatic, RAII-based,
safe C++20 while preserving byte-for-byte identical CLI output and behavior.

- Original authors: Maximilian Kreuzer and Harald Skarke (GPLv3).
- Upstream: https://gitlab.com/stringstuwien/PALP
- PALP website: http://hep.itp.tuwien.ac.at/~kreuzer/CY/CYpalp.html
- PALP online documentation: http://palp.itp.tuwien.ac.at/wiki/index.php/PALP_online_documentation

## Build system

CMake only. No legacy Makefiles remain.

```bash
cmake -S . -B build -D CMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Key CMake options:

- `POLY_Dmax` (default 6): compile-time maximum polytope dimension. Larger
  values increase memory use because many arrays are sized from it.
- `PALP_TEST_DIMENSIONS` (default `"4;5;6;11"`): dimensions for which
  dimension-suffixed executables and CTest tests are built.
- `PALP_TESTS_LONG` (default OFF): enable long-running tests.

Sanitizer builds are documented in `PLAN.md`.

## Test policy

**Every change must keep the 196-test CTest suite passing** (Release + ASan).
The suite is in `tests/*.sh` and uses `tests/lib/run-test.sh`. Each test script
is run once per configured dimension (`PALP_TEST_DIMENSIONS`) with `DIM=<dim>`
set in the environment.

Do not modify expected test outputs unless the output change is intentional
and equivalent behavior is preserved. The modernization targets byte-for-byte
output compatibility.

## Code style and conventions

- **Language**: C++20, `-Wall -Wextra`, no extensions.
- **Formatting**: run `clang-format -i -style=LLVM` on changed files before
  every commit.
- **Naming**: LLVM style.
  - Functions/variables: `camelCase`
  - Structs/classes/aliases: `PascalCase`
  - Macros and compile-time dimension limits: `UPPER_CASE`
  - Member fields: `camelCase` (no leading underscore)
  - Parameters: `camelCase` (drop the legacy `_P`, `_V`, `_E` prefixes)
  - Keep single-letter math variables (`d`, `n`, `v`) where they mirror paper
    notation.
- **Headers**: public headers live under `include/palp/`, use `#pragma once`.
  Legacy private headers are in `src/`.
- **Memory**: prefer RAII. Avoid `malloc`/`free`/`new`/`delete` in new code.
  Use `std::vector`, `std::array`, `std::unique_ptr`, `std::string`.
- **Strings/I/O**: `printf`/`fprintf` calls that produce CLI output must stay
  byte-for-byte compatible. Internal/diagnostic prints may use `<format>`.
- **Error handling**: new library code should avoid `exit(1)`. Prefer
  `throw std::runtime_error` or return codes; `exit()` is confined to driver
  `main()` functions. There is a `[[noreturn]] void Die(const char*)` helper in
  `cws.cpp` that may serve as a model for a future `palp::die()` helper.
- **C headers**: use `<cXXX>` equivalents. POSIX headers (`<unistd.h>`,
  `<sys/*.h>`, `<fcntl.h>`) are isolated in `include/palp/PosixProcess.h`.
- **No `NULL`**: use `nullptr`.
- **No C-style casts**: use `static_cast`/`const_cast`/`reinterpret_cast`, or
  eliminate the cast by changing types.
- **No `typedef struct`**: use plain `struct Name { ... };` or `using` aliases.
- **No function-local `static` mutable state**: move to anonymous-namespace
  variables or thread through call sites.
- **`POLY_Dmax` stays compile-time** (macro). Do not introduce dynamic
  dimensioning in the migration phase.

## Important types and constants

- `Long` and `LLong` are `using` aliases for `long` and `long long` in
  `include/palp/palp_types.h`. Both are expected to be 64-bit
  (`static_assert`s enforce this).
- Dimension limits are macros driven by `POLY_Dmax`: `POINT_Nmax`,
  `VERT_Nmax`, `FACE_Nmax`, `SYM_Nmax`, `EQUA_Nmax`, `AMBI_Dmax`.
- `PalpContext` (`include/palp/Global.h`) is the global I/O context used by
  legacy code after Phase 5.11 removed `inFILE`/`outFILE` globals.

## Source layout

| Directory | Contents |
| --------- | -------- |
| `src/` | All implementation files (`*.cpp`) and legacy private headers |
| `include/palp/` | Public C++ headers (`Global.h`, `Subpoly.h`, `Mori.h`, `Nef.h`, `LG.h`, `Rat.h`, `palp_types.h`, `PosixProcess.h`) |
| `tests/` | Shell test scripts and reference inputs |
| `build/` | Default CMake build directory (not tracked) |

## Executables

The build produces both unsuffixed executables (`poly.x`, `class.x`, ...) and
dimension-suffixed copies (`poly-4d.x`, `class-4d.x`, ...). For the 4D
reflexive classification, use `class-4d.x` and `cws-4d.x`, which are compiled
with `POLY_Dmax=4`.

## Planning and issue tracking

- `PLAN.md` contains the phased modernization plan. Completed phases are
  marked `[x]`; open phases are `[ ]`.
- `ISSUES.md` tracks audit findings and bugs. Each issue has a severity,
  status, file/line, and description. Closed issues keep their numbers.
- Prefer one concern or one file per commit. Update `ISSUES.md` and `PLAN.md`
  when closing issues or completing phases.

## Common gotchas

- Output compatibility is the top constraint. Do not change `printf`/`fprintf`
  format strings or ordering unless you are certain the test outputs still
  match.
- `class.x` for the 4D reflexive classification must be run with
  `POLY_Dmax=4` (`class-4d.x`).
- `nef.x` may require `POLY_Dmax >= dim(N) + codim - 1`.
- `mori.x -MD` may require `POLY_Dmax >= n_points - dim(N) - 1`.
- The codebase still contains many `exit(1)` calls in library functions.
  Converting them is planned (Phase 14) but is a large refactor; do not mix
  it casually with unrelated fixes.
- Several `goto` statements and C-style casts remain; address them in their
  dedicated phases (13.1–13.2).

## When in doubt

- Re-run `cmake --build build && ctest --test-dir build` after every edit.
- Keep changes minimal and focused.
- Update this file if you change build steps, conventions, or project
  structure.
