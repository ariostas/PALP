# Modernization Plan

Goal: move PALP from C toward maintainable modern C++ while preserving the historical outputs, especially for 4D reflexive polytopes. Each step should be small enough to verify with the existing shell tests plus targeted new tests.

## Phase 0: Preserve Current Behavior

1. Establish a clean baseline.
   - Run `make cleanall && make -jN`.
   - Run `make check`; run `make checklong` separately when time permits.
   - Capture compiler version, flags, and test output.

2. Make tests easier to run in CI.
   - Add CMake/CTest entries that call the existing `tests/*.sh` scripts.
   - Keep `GNUmakefile check` working.
   - Verification: `ctest` and `make check` run the same test scripts.

3. Add sanitizer build modes.
   - Add CMake presets or documented commands for AddressSanitizer and UndefinedBehaviorSanitizer.
   - Do not change production flags yet.
   - Verification: default build still matches baseline; sanitizer failures are filed in `ISSUES.md`.

## Phase 1: Make C Safer Before C++

4. Replace known undefined behavior without changing interfaces.
   - Fix the `LG.c` use-after-free.
   - Fix 64-bit incidence shifts.
   - Fix unchecked `fgets()` in `SingularInput.c`.
   - Verification: targeted regression tests plus full `make check`.

5. Convert assertion-only user checks into explicit runtime checks.
   - Start with parsers and fixed-size array boundaries.
   - Keep assertions for internal invariants.
   - Verification: add tests for invalid input and capacity errors.

6. Normalize matrix and array shape contracts.
   - Introduce named C typedefs or small wrapper structs for `PairMat`, normal-form matrices, and incidence arrays.
   - Eliminate call sites that pass smaller arrays to larger formal types.
   - Verification: compiler warnings for `-Warray-parameter` and `-Wstringop-overflow` decrease.

7. Remove stale build paths.
   - Either delete/update `Makefile` or make it delegate to `GNUmakefile`.
   - Make CMake and GNUmakefile source lists match exactly.
   - Verification: both build systems produce all five main executables.

## Phase 2: Prepare for C++ Compilation

8. Add include guards and C++ compatibility shims.
   - Guard every header.
   - Remove duplicate local declarations where headers should be authoritative.
   - Add `extern "C"` only if C and C++ objects must coexist temporarily.
   - Verification: C build unchanged; headers can be included from a trivial `.cc` file.

9. Introduce fixed-width type aliases.
   - Replace macro aliases like `Long` and `LLong` with a central typedef/using layer.
   - Document required ranges for classified 4D workflows.
   - Verification: binary output remains unchanged on baseline tests.

10. Compile selected modules as C++ without semantic changes.
    - Start with leaf-like modules: `Rat.c`, then `Coord.c`, then `Vertex.c`.
    - Rename one file at a time only after it compiles cleanly as C++.
    - Verification: full test suite after each file.

## Phase 3: Introduce Modern C++ Types

11. Wrap rational numbers.
    - Replace free-function `Rat` operations with a small value type.
    - Add explicit construction and checked division.
    - Verification: unit tests for arithmetic and all existing CLI tests.

12. Wrap bounded arrays.
    - Use `std::array` for compile-time bounded vectors/matrices.
    - Provide `.size()`-checked helpers before changing algorithms.
    - Verification: no output changes; sanitizer builds improve.

13. Replace manual allocation with RAII.
    - Use `std::vector` for variable-size work buffers.
    - Use automatic storage or `std::unique_ptr` where fixed historical layout is still required.
    - Verification: run sanitizers and compare CLI outputs.

14. Introduce a runtime context.
    - Replace global `inFILE`/`outFILE` and scattered option globals with a context object.
    - Migrate one executable at a time.
    - Verification: CLI behavior and prompts remain compatible.

## Phase 4: Modularize Algorithms

15. Separate parsing, core algorithms, and presentation.
    - Give `poly`, `cws`, `nef`, `mori`, and `class` thin CLI front ends.
    - Move reusable computation into library modules.
    - Verification: add direct library tests for normal form, vertex/facet computation, and Hodge output.

16. Make dimensions explicit.
    - Preserve `POLY_Dmax` builds initially.
    - Later evaluate template-based dimensions or runtime-sized containers for selected algorithms.
    - Verification: existing 4d/5d/6d/11d test matrix remains green.

17. Add golden datasets for historically critical outputs.
    - Keep small inputs in git.
    - For larger 4D reflexive checks, store hashes or sampled canonical outputs.
    - Verification: modernization changes must reproduce golden outputs exactly unless a bug fix intentionally changes them.

## Phase 5: Finish C++ Cleanup

18. Replace macros with scoped constants/functions.
    - Convert `min`, `max`, capacity macros where feasible.
    - Remove macro side effects from headers.
    - Verification: warnings stay low and tests pass.

19. Replace `exit()` from library code with errors.
    - Use result types or exceptions at the C++ library boundary.
    - Keep CLI behavior by translating errors into the same messages and exit codes.
    - Verification: invalid-input tests preserve expected diagnostics.

20. Enforce formatting and static checks.
    - Add `clang-format` after large mechanical moves are complete.
    - Add warning gates gradually, starting with warnings that are already fixed.
    - Verification: CI rejects new warnings in migrated modules.

## Verification Rule

After each step, record:

- Build command and compiler version.
- `make check` or `ctest` result.
- Any sanitizer result if applicable.
- Whether CLI output changed intentionally or unintentionally.
