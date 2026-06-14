# Potential Logic Issues and Bug Candidates

This file records issues found during an initial code-reading and build-warning pass. These are not all proven defects; each item is a candidate for focused reproduction, minimization, and regression tests before changing behavior.

## 1. `LG.c` returns through a freed `EqList`

- Location: `LG.c:15-20`
- Evidence: `Is_Gen_CY()` allocates `E`, uses `E->ne`, calls `free(E)`, and then returns `i==E->ne`.
- Risk: This is a real use-after-free. It may usually work because the freed block is not overwritten immediately, but the result is undefined behavior and could change under C++, sanitizers, different allocators, or small refactors.
- Suggested check: Save `E->ne` before `free(E)`, add a test covering `Is_Gen_CY()`, and run under AddressSanitizer.
- Status: Fixed in Phase 1 item 4 by saving `E->ne` before `free(E)`; verified with focused `poly -l` tests and `make check`.

## 2. Incidence bit operations use 32-bit shifts in a 64-bit type

- Location: `MoriCone.c:97-101`, `MoriCone.c:1734-1737`
- Evidence: `makeN()` correctly uses `((Inci64) 1) << N`, but `putN()`, `setN()`, and the origin-entry clearing code use `1 << N`.
- Risk: For `N >= 32`, shifting a signed `int` is undefined or truncating before assignment to `Inci64`. `VERT_Nmax` is 64, and Mori triangulation code explicitly handles `p <= 64`, so this can corrupt incidences for larger vertex counts.
- Suggested check: Change all incidence shifts to `(Inci64)1 << N`, add tests for bit positions 31, 32, and 63, and audit all `1 <<` expressions touching `Inci64`.
- Status: Fixed in Phase 1 item 4 by casting incidence masks to `Inci64` before shifting; verified with focused Mori tests and `make check`.

## 3. Runtime correctness depends on assertions being enabled

- Location: `CMakeLists.txt:14-15`, many `assert()` checks across core algorithms
- Evidence: CMake explicitly strips `-DNDEBUG` because compiling with it breaks the code. Bounds checks and even important input validations are often only `assert()` calls, for example `Coord.c:178-179`, `MoriCone.c:1528-1532`, and many fixed-array writes.
- Risk: Release builds or downstream embedders that define `NDEBUG` can silently remove checks and turn malformed or simply large inputs into memory corruption.
- Suggested check: Introduce explicit checked preconditions for user input, capacity limits, and allocation results. Reserve `assert()` for internal invariants after the explicit checks exist.

## 4. Pairing-matrix buffers are passed with incompatible dimensions

- Location: `Global.h:149`, `Global.h:426`, `Polynf.c:1475-1478`, `Polynf.c:1729-1752`, `Polynf.c:1811`, `E_Poly.c:339`, `Subdb.c:1718`
- Evidence: `PairMat` is `Long[EQUA_Nmax][VERT_Nmax]`. For `POLY_Dmax > 4`, `EQUA_Nmax` defaults to 1280, but several call sites allocate `Long X[VERT_Nmax][VERT_Nmax]` and pass it to `Make_VEPM()`. GCC warns that `Make_VEPM` may access 655360 bytes in a 32768-byte region.
- Risk: If an `EqList` at one of these call sites can have more than `VERT_Nmax` equations, `Make_VEPM()` can write past the local matrix. This is especially relevant outside the fully classified 4D reflexive case.
- Suggested check: For each call site, prove and assert `E.ne <= VERT_Nmax`, or allocate a true `PairMat`. Then encode the matrix shape in a wrapper type.
- Status: Fixed for the known `Make_VEPM()` call sites in Phase 1 item 6 by using `PairMat` and by making the `Complete_Poly()` prototype match its implementation. Related `Make_ANF()`, `Read_HyperSurf()`, and `SimplexVolume()` shape warnings were also normalized.

## 5. Reflexive fibration recursion appears to overrun fixed arrays

- Location: `Polynf.c:1442-1472`
- Evidence: GCC repeatedly reports out-of-bounds accesses to `s`, `T`, and `G` in `Fiber_Rec_New_Point()` when `make check` builds dimensions 4, 5, and 6. The recursion indexes `s[r]`, `T[r]`, and `G[r]` while `r` can reach the requested fiber dimension.
- Risk: This may be either a false positive from aggressive inlining or a real off-by-one when `fdim >= POLY_Dmax` or when recursion reaches `r == fdim`. If real, fibration output can be corrupted.
- Suggested check: Add explicit guards on `fdim`, `r`, and array bounds. Create a minimal test for `poly.x` fibration options at dimensions 4, 5, and 6 under AddressSanitizer.

## 6. Combined-weight-system and weight parsers write before all bounds are checked

- Location: `Coord.c:139-179`, `Coord.c:285-321`, `LG.c:32-74`, `LG.c:191-204`
- Evidence: Input numbers are collected into fixed arrays and then mapped into fixed-size `CWS`/`Weight` fields. Several important capacity checks are assertions, and some dimensions are checked only after partially filling arrays.
- Risk: Malformed or unusually large CWS input can trigger assertion aborts in normal builds and memory corruption if assertions are disabled. Failed `fscanf()` calls can also leave local integers uninitialized.
- Suggested check: Replace the scanf loops with a checked tokenizer that validates token count, conversion success, sign, and destination capacity before writing into fixed arrays.
- Status: Fixed for the listed parser paths in Phase 1 item 5 by adding checked integer reads, explicit weight/count diagnostics, `/Z` capacity checks, and invalid-input regressions. A broader tokenizer abstraction can still be introduced later during parser modernization.

## 7. Rational arithmetic assumes nonzero denominators and can divide by zero

- Location: `Rat.c:7-13`, `Rat.c:23-51`, `Rat.c:81-90`, `Rat.c:171-177`, `Rat.c:187-215`, `Rat.c:245-254`
- Evidence: `rR(a,b)` calls `Fgcd(a,b)` and divides by `g` without checking `b != 0`. `Fgcd()` itself starts with `a %= b`, which is undefined when `b == 0`. Similar assumptions exist in the `LRat` path.
- Risk: Unexpected zero denominators from matrix reduction or input parsing can crash or corrupt computations instead of producing a controlled error.
- Suggested check: Add denominator preconditions to rational constructors/division, use `NNgcd()` where zeros are valid, and add tests for zero-denominator rejection.

## 8. `DYNComplete_Poly()` can read past ordered facets if rank is deficient

- Location: `E_Poly.c:60-148`
- Evidence: The loop `while (rank<n)` increments `i` and indexes `_E->e[OrdFac[i]]` without checking `i < _E->ne`. GCC also warns that `BasFac` may be used uninitialized later.
- Risk: Non-full-dimensional, malformed, or numerically degenerate equation sets can read past `OrdFac` and leave basis data uninitialized.
- Suggested check: Add an explicit failure path when no full-rank facet basis exists, and build a focused malformed-input test.

## 9. Temporary hypersurface input file is global and fragile

- Location: `SingularInput.c:49-87`, `SingularInput.c:172-179`
- Evidence: `Read_HyperSurf()` uses the fixed filename `HEInput.txt` in the working directory to persist hypersurface input across calls.
- Risk: Parallel runs, interrupted runs, or two `mori.x -H` invocations in the same directory can read stale or cross-contaminated hypersurface data.
- Suggested check: Store this state in memory on `MORI_Flags` or use a per-process temporary file that is removed deterministically.

## 10. `Read_HyperSurf()` can parse uninitialized stack content on EOF or read failure

- Location: `SingularInput.c:49-87`
- Evidence: `fgets()` return values are not checked before `strtok(string, ...)`.
- Risk: Empty stdin or an unreadable cached file can leave `string` uninitialized and parse arbitrary stack bytes.
- Suggested check: Check `fgets()` and fail with a clear diagnostic if no line was read.
- Status: Fixed in Phase 1 item 4 by checking `fgets()` in both input paths before calling `strtok()`; verified through the full `make check` suite.

## 11. `cws.c` weight generation triggers array-bound warnings

- Location: `cws.c:511`, `cws.c:601-620`
- Evidence: GCC reports that recursive `createweights()` can access `X->points[npoints]` above the `points[W_Nmax][W_Nmax]` bounds.
- Risk: This may be a false positive caused by recursion and compile-time dimension constants, but it is in classification-related code. If real, generated weight systems can be missed or corrupted.
- Suggested check: Add `assert(0 <= npoints && npoints < W_Nmax)` before writes and a separate explicit stopping condition before `testweisys()` is called with `npoints > W_Nmax`.

## 12. Legacy `Makefile` is stale and likely broken

- Location: `Makefile:14-15`, `Makefile:48`
- Evidence: It refers to `Moricone.c` instead of `MoriCone.c` and has a malformed `mori` rule with dependencies and the link command on the same line. The maintained build appears to be `GNUmakefile` and CMake.
- Risk: Users or agents invoking `make -f Makefile` can get a broken or incomplete build and then debug the wrong thing.
- Suggested check: Remove the stale file, redirect it to `GNUmakefile`, or update it to match the maintained source list.
- Status: Fixed in Phase 1 item 7 by replacing `Makefile` with a compatibility wrapper that delegates to `GNUmakefile`; verified with wrapper builds, CMake, and `make -f Makefile check`.

## 13. Global `FILE *` input/output state limits reentrancy

- Location: `Global.h:104`, definitions in `poly.c`, `class.c`, `nef.c`
- Evidence: Most parser and printer functions implicitly use global `inFILE` and `outFILE`.
- Risk: This makes library bindings, nested calls, tests, and future C++ wrappers fragile. It also makes error handling via `exit()` harder to replace incrementally.
- Suggested check: Introduce an execution/context object that owns input, output, options, and diagnostics, then migrate one CLI at a time.

## 14. AddressSanitizer build changes `nef -v` point statistics output

- Location: `nef.c:305-311`, exposed by `tests/6.4.18-nef-v.sh`
- Evidence: `ctest --preset asan` passed 193/196 tests but failed `6.4.18-nef-v.sh` for dimensions 5, 6, and 11. The expected point-statistics lines are `35# 1` and `100# 1`; the ASan build instead prints every index from `0` through `511` with values around `-1094795586`, except the two expected bins are offset by one. The same tests pass in the release and UBSan presets.
- Risk: This looks like an initialization or lifetime bug in the `Pstat` path that only becomes visible under ASan-instrumented stack/heap layout. It may make the `nef -v` summary depend on compiler instrumentation even when no sanitizer trap is emitted.
- Suggested check: Audit `Pstat` allocation and initialization before `Print_Pstat()`, add a targeted test for `nef -v` point counts under ASan, and check whether `Pstat::P` is fully zeroed before `Print_VP()` increments bins.

## 15. `MakeMobius()` stores `int **` into under-aligned `int` storage

- Location: `LG.c:753-755`
- Evidence: `ctest --preset ubsan` passed 193/196 tests but failed `3.2.11-poly-l.sh` for dimensions 5, 6, and 11. UBSan reports `runtime error: store to misaligned address ... for type 'int *', which requires 8 byte alignment` at `LG.c:754` in `MakeMobius()`, called from `LGO_VaHo()` and then `poly.c:212`. The failing command is the `/Z3` `poly -fl` example from the PALP manual.
- Risk: `MakeMobius()` allocates one raw block sized as `int` data plus `int *` data, then casts the middle of the `int` array to `int **`. On 64-bit systems this can place pointer slots at only 4-byte alignment. That is undefined behavior and can break under stricter architectures, sanitizers, or C++ allocators.
- Suggested check: Split the allocation into separately aligned `int *d` and `int **mt` allocations, or allocate a struct/byte buffer with explicit alignment. Preserve the existing triangular indexing and rerun `ctest --preset ubsan`.

## 16. `poly -A` asserts when `POLY_Dmax` equals the input dimension

- Location: `Polynf.c:2503-2504`, reached from `poly.c:247-250`
- Evidence: `echo '5 1 1 1 1 1' | ./poly-4d.x -fA` aborts with `Make_ANF: Assertion 'P->n<POLY_Dmax' failed`. The same smoke test succeeds for `poly-5d.x`, `poly-6d.x`, and `poly-11d.x`.
- Risk: The affine normal form path temporarily adds one coordinate row, so it needs `POLY_Dmax > P->n`. As written, this is a user-facing capacity condition enforced only by `assert()`, producing an abort instead of a clear diagnostic.
- Suggested check: Replace the assertion with an explicit `POLY_Dmax` capacity check, add an invalid-input regression for `poly-4d.x -A`, and verify `poly -fA` output remains unchanged for dimensions with enough headroom.

## Build observations from this pass

- `make -j2` completed and produced the default `.x` executables.
- `make check` initially failed at the Mori/Singular tests because `Singular` was not installed. After installing `singular` 4.4.1 into the `vibe` conda environment, `make check` completed successfully.
- Notable warning categories: use-after-free, possible array bounds violations, matrix-shape mismatches, uninitialized values, misleading indentation, and `abs(long)` truncation.
