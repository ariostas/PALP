# Agent Guide

## Project Shape

PALP is a C codebase with five main CLI programs:

- `poly.x`: polytope analysis.
- `class.x`: classification/database workflows.
- `cws.x`: combined weight systems.
- `nef.x`: nef partitions.
- `mori.x`: triangulations, Mori cone, and Singular integration.

Core shared files are `Global.h`, `Coord.c`, `Rat.c`, `Vertex.c`, and `Polynf.c`. Specialized modules include `LG.c`, `Subpoly.c`, `Subadd.c`, `Subdb.c`, `E_Poly.c`, `Nefpart.c`, `MoriCone.c`, and `SingularInput.c`.

## Build and Test

Preferred build:

```sh
make -j2
```

Full regression suite:

```sh
make check
```

`mori.x` tests require `Singular` on `PATH`. In this workspace it was installed into the `vibe` conda environment:

```sh
conda install -n vibe -c conda-forge singular
```

Long tests:

```sh
make checklong
```

CMake build:

```sh
cmake -S . -B build -D CMAKE_BUILD_TYPE=Release
cmake --build build
```

C++ header smoke test:

```sh
c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc
```

The maintained make-based build is `GNUmakefile`. The top-level `Makefile` is only a compatibility wrapper that delegates to `GNUmakefile`.

## Important Build Details

- `POLY_Dmax` controls many compile-time array sizes and must be consistent across all object files.
- `GNUmakefile` handles per-dimension builds such as `poly-4d.x`, `poly-5d.x`, `poly-6d.x`, and `poly-11d.x`.
- CMake removes `-DNDEBUG` because this code relies on assertions for runtime correctness. Do not add `NDEBUG` casually.
- Object files and `.x` executables are generated in the repository root.

## Coding Constraints

- Preserve exact CLI output unless the task explicitly changes behavior.
- Treat `assert()` as part of current behavior until explicit runtime checks replace it.
- Avoid broad refactors when fixing a bug. Add a focused test first when practical.
- Be careful with fixed-size arrays from `Global.h`: `POINT_Nmax`, `VERT_Nmax`, `FACE_Nmax`, `EQUA_Nmax`, `FIB_Nmax`, and `AMBI_Dmax`.
- Many functions use global `FILE *inFILE` and `FILE *outFILE`; changing this requires coordinated CLI testing.
- The code uses `Long` and `LLong` macros for numeric assumptions. Do not change widths without checking historical output.

## Test Protocol

Tests are shell scripts under `tests/`. Each script sources `tests/lib/run-test.sh`, sets:

- `COMMAND`
- `DESCRIPTION`
- `EXPECTED`

and then calls `run_test`.

`make check` runs scripts across multiple `DIM` values. Individual tests can usually be run directly, for example:

```sh
DIM=6 tests/2.1-polytope-input.sh
```

## Current High-Risk Areas

See `ISSUES.md` before making modernization changes. In particular:

- `LG.c` has a use-after-free candidate.
- `MoriCone.c` mixes 64-bit incidence storage with 32-bit shifts.
- `Polynf.c`, `E_Poly.c`, and `Subdb.c` have matrix-shape and bounds warnings.
- Parser code often uses fixed arrays and assertion-only validation.

## Suggested Workflow

1. Read the relevant header first, especially `Global.h`.
2. Build once with `make -j2`.
3. Make the smallest behavior-preserving change.
4. Run the most targeted test script.
5. Run `make check` before considering the change complete.
6. Record any new suspected bug in `ISSUES.md` instead of mixing investigation with unrelated refactoring.
