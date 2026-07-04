# PALP Potential Issues / Bugs

Found during code inspection for the C → C++17 migration. Each item has a file,
line number, severity, and description. Fixes are applied in Phase 5 (after
migration), one regression-tested fix at a time.

Line numbers refer to the original C source as of the start of the
modernization effort. They will shift as migration proceeds; the function
names and descriptions should remain locatable.

---

## Critical (data correctness risk)

### 1. Side effect inside `assert()`
- **File**: `MoriCone.c`
- **Line**: 408
- **Severity**: Critical
- **Description**: `assert(++m < binco)` — the increment of `m` is a side
  effect inside `assert()`. If `NDEBUG` is defined (which CMake currently
  strips, but this is a latent bug), `m` is never incremented, causing an
  infinite loop or wrong results. Must move the increment out of the assert:
  `++m; assert(m < binco);`

### 2. Misaligned store in `MakeMobius`
- **File**: `LG.c`
- **Line**: 753–754
- **Severity**: Critical
- **Description**: `MakeMobius` allocates a single buffer of `int`s and then
  casts `&d[(n*(n+3))/2]` to `int **`. When `Pint` is `int` (the default for
  `POLY_Dmax < 7`), the `int **` array is stored at an offset that may not be
  aligned to `sizeof(int *)`. UBSAN reports:
  `runtime error: store to misaligned address ... for type 'int *', which requires 8 byte alignment`.
  This is undefined behavior and can cause crashes or wrong Möbius function
  values on strict-alignment architectures. The data structure should be split
  into two separately allocated arrays (`d[]` and `mt[]`) rather than packing
  them into one buffer.

### 3. Memory leaks on error paths in `Find_Equations` / `IP_Check` / `Ref_Check`
- **File**: `Vertex.cpp`
- **Lines**: 984–1125
- **Severity**: Critical
- **Status**: Fixed
- **Description**: `CEq`, `CEq_I`, and `F_I` were originally raw `malloc`ed
  buffers. They are now `std::unique_ptr<CEqList>` and `std::array<INCI,...>`,
  so memory is released automatically on every exit path including the
  `exit(1)` calls in subroutines. Remaining bounds checks that previously used
  `assert` were also converted to explicit `fputs`/`exit(1)` errors.

### 4. `realloc` loses old pointer on failure
- **File**: `E_Poly.cpp`
- **Line**: ~60
- **Severity**: High
- **Status**: Fixed
- **Description**: `DYNadd_for_completion` previously used `realloc` directly on
  `_CP->L`. During migration the buffer was converted to `std::vector`, so resize
  failures are handled by the standard library and the old pointer is no longer
  lost. A capacity check (`_CP->np >= _CP->L.size()`) triggers `resize()` before
  the new element is written.

### 5. Partial array initialization
- **File**: `E_Poly.cpp`
- **Line**: ~722, ~1463
- **Severity**: Medium
- **Status**: Fixed
- **Description**: `int h[POLY_Dmax][POLY_Dmax] = {}` now used consistently for
  full zero-initialization in both `Make_Mirror` and `Compute_E_Poly`.

### 6. `assert` guarding division by zero
- **File**: `Polynf.cpp`
- **Line**: ~1528, ~1985, ~3720
- **Severity**: High
- **Status**: Fixed
- **Description**: `assert(g > 0)` was followed by division by `g`. Converted
  to explicit `if (g <= 0) { fputs(...); exit(1); }` checks in
  `Aux_Vol_Barycent`, `Add_Square_To_Rel`, and `Normalize_QuotientZ`.

### 7. Dead function: `FileRW()`
- **File**: `cws.c`
- **Line**: 1713
- **Severity**: Medium
- **Description**: `FileRW()` is declared and defined but never called. The
  function takes a `FILE*` by value (so any opened file is never returned to
  the caller). Dead code that should be removed.

---

## Error handling (exit code / robustness)

### 8. `exit(0)` used for all error conditions
- **File**: ALL files (~200+ instances)
- **Lines**: everywhere
- **Severity**: High
- **Description**: `exit(0)` (success exit code) is used for ALL error
  conditions across the entire codebase. Exit code 0 means success; errors
  should return non-zero. Tests check stdout so this is currently masked,
  but scripts wrapping PALP cannot detect failures. All instances should be
  `exit(1)` or a proper exception.

### 9. Memory leaks in `mori.c` main loop
- **File**: `mori.cpp`
- **Lines**: 89–103
- **Severity**: High
- **Status**: Fixed
- **Description**: `CW`, `E`, `DE`, `_P`, and `_DP` were raw heap allocations.
  During migration they were converted to `std::unique_ptr`, so memory is
  released automatically when the pointers go out of scope. `PM` and `DPM`
  are stack-allocated `PairMat` arrays (typedef'd as `Long[EQUA_Nmax][VERT_Nmax]`)
  and do not leak.

### 10. Huge struct never freed in `RgcWeights` — RESOLVED by migration
- **File**: `cws.cpp`
- **Line**: 530
- **Severity**: High
- **Description**: `RgcWeights` allocated `RgcClassData *X` which contains
  `Equation wli[WDIM]` where `WDIM = 800000`. This was a massive allocation that
  was never freed. During migration it was replaced by
  `std::unique_ptr<RgcClassData>`, so the memory is now released automatically.
  No separate Phase 5 fix is needed.

### 11. Missing NULL check on `malloc`
- **File**: `SingularInput.c`
- **Line**: 173
- **Severity**: Medium
- **Status**: Fixed
- **Description**: `he = malloc(divclassnr * sizeof(int))` — no NULL check
  on `malloc` return. Subsequent use of `he` could segfault. During migration
  this buffer was replaced by `std::vector<int> he(divclassnr)`, so allocation
  failure throws `std::bad_alloc` instead of producing a NULL pointer.

### 12. Command injection via `system()`
- **File**: `SingularInput.cpp`
- **Line**: ~99–601
- **Severity**: Medium
- **Status**: Fixed
- **Description**: `system(SingularCall)` used a command string built from
  `TMPDIR` and the temp file path, allowing shell metacharacters. Replaced with
  `fork`/`dup2`/`execvp` of `Singular` with explicit argv (`-q`, script file),
  and added a check for `mkstemp` failure.

---

## Memory safety

### 13. Memory leak on `exit(0)` in `Eval_Poly_NF`
- **File**: `Polynf.c`
- **Line**: 262
- **Severity**: Medium
- **Description**: `Eval_Poly_NF` mallocs `CL`; if called functions `exit(0)`
  (which many do), `CL` is leaked. This is a systemic issue: any `exit(0)`
  inside a function called from a malloc-owner leaks that allocation.

### 14. Memory leaks on error paths in `Fano5d`
- **File**: `Polynf.cpp`
- **Lines**: 4826–5096
- **Severity**: Medium
- **Status**: Fixed
- **Description**: `Fano5d` allocated `Q` and `F` with raw `malloc` and `M`/`G`
  with `Init_Matrix`, but freed them only on the success path. Replaced `Q` and
  `F` with `std::unique_ptr` and let RAII free them on every return path.
  `M`/`G` are still stack-allocated `Matrix` structures whose internal buffers
  are freed by `Free_Matrix` on the success path.

### 15. Static mutable state in `Subdb.c`
- **File**: `Subdb.c`
- **Line**: 1596
- **Severity**: Low
- **Description**: `static unsigned char uc[NUC_Nmax]` and `static int ms3` —
  persistent mutable state, not thread-safe, not reentrant.

### 16. Static state in `LG.c`
- **File**: `LG.c`
- **Lines**: 280, 491
- **Severity**: Low
- **Description**: `static int MaxPoNum` and `static int M` — persistent state
  across calls. Function behavior depends on previous calls.

### 17. Static counters in `Polynf.c`
- **File**: `Polynf.c`
- **Lines**: 1168, 1947, 2635
- **Severity**: Low
- **Description**: `static int` counters in `ConifoldSing`, `GL_Lattice_Basis`,
  and `Fano5d` — not thread-safe, state persists unexpectedly across calls.

### 18. Global mutable pointer in `Subadd.c`
- **File**: `Subadd.c`
- **Line**: 96
- **Severity**: Low
- **Description**: `NF_List *AuxNFLptr = NULL` — global mutable pointer used as
  a "dirty trick" for statistics. Creates hidden coupling between functions.

---

## Integer / type issues

### 19. Integer overflow in volume computation
- **File**: `Polynf.cpp`
- **Lines**: ~1438, ~1450
- **Severity**: Medium
- **Status**: Fixed
- **Description**: `I *= VM[i][i]` — product of diagonal entries could overflow
  `Long`. Added overflow checks before each multiplication and hard-error
  `exit(1)` if the intermediate product would exceed `std::numeric_limits<Long>::max()`.

### 20. `lcm` macro overflow
- **File**: `lgotwist.cpp`
- **Line**: ~56
- **Severity**: Medium
- **Status**: Fixed
- **Description**: `Lcm(a,b)` computed `a * (b / gcd(a,b))`, risking overflow
  even when the result would fit. Rewrote to divide first: `(a / gcd) * b`,
  and handle `gcd == 0`.

### 21. `printf` format truncation
- **File**: `Rat.c`
- **Line**: 59
- **Severity**: Low
- **Status**: Fixed
- **Description**: `Rpr` uses `%d` format with `(int) c.N` — truncation if
  `Long` exceeds `int` range. Changed to `%ld` with `(long) c.N`/`(long) c.D`.

### 22. Platform-dependent integer sizes
- **File**: `Global.h`
- **Lines**: 12–13
- **Severity**: Low
- **Description**: `Long` defined as `long`, `LLong` as `long long` — sizes
  are platform-dependent (32 vs 64 bit for `long`). Should use `<cstdint>`
  types (`int64_t`, etc.) in C++ version for consistency.

### 23. Unusual declaration order
- **File**: `nef.c`
- **Line**: 485
- **Severity**: Low
- **Status**: Fixed
- **Description**: `int long nl` — unusual declaration order (legal C but
  confusing). Changed to `long int nl`.

---

## Input validation / buffer safety

### 24. Off-by-one in `char c[999]` buffer
- **File**: `LG.c`
- **Lines**: 34, 50–54
- **Severity**: Medium
- **Description**: `char c[999]` fixed buffer. The loop at line 50 reads
  `c[n] = fgetc(inFILE)` before checking `n == 999` at line 54. If exactly 999
  characters are read, the write to `c[999]` is out of bounds (off-by-one).
  The check should be before the write, or the buffer should be `[1000]`.

### 25. No bounds check on skeleton input
- **File**: `lgotwist.c`
- **Line**: 167
- **Severity**: Medium
- **Description**: `s->p[s->N++] = i - '0'` — no bounds check against `NM`
  (9 or 13). Buffer overflow if input has more elements than expected.

### 26. `sprintf` buffer overflow risk
- **File**: `cws.c`
- **Line**: 1882
- **Severity**: Low
- **Description**: `char command[100]` used with `sprintf` — potential buffer
  overflow if the formatted string exceeds 99 characters. Should use
  `snprintf`.

### 27. Tight `sprintf` buffer
- **File**: `Subdb.c`
- **Line**: 1468
- **Severity**: Low
- **Description**: `char com[35]` with `sprintf` — if the format produces more
  than 34 characters, buffer overflow. Should use `snprintf`.

### 28. Unchecked `fscanf` return values
- **File**: ALL files
- **Lines**: many
- **Severity**: Medium
- **Description**: `fscanf` return values are never checked anywhere in the
  codebase. Malformed input causes undefined behavior (variables left
  uninitialized, loops reading garbage). All `fscanf` calls should check
  return values and handle failures gracefully.

---

## Code smells / maintainability

### 29. Local `min`/`max` macro redefinition
- **File**: `E_Poly.c`
- **Lines**: 6–7
- **Severity**: Low
- **Description**: `#define min(a,b)` and `#define max(a,b)` redefined locally
  — conflicts with system headers and other PALP files that define the same
  macros. Should use `std::min`/`std::max` in C++.

### 30. `min`/`max` macros in header
- **File**: `LG.h`
- **Lines**: 5–6
- **Severity**: Low
- **Description**: `#define min(a,b)` / `#define max(a,b)` in a header file
  — pollutes all including files. Should be removed; use `std::min`/`std::max`.

### 31. `min`/`max` macros in `Subpoly.h`
- **File**: `Subpoly.h`
- **Lines**: 220–221
- **Severity**: Low
- **Description**: Same `#define min(a,b)` / `#define max(a,b)` redefinition.
  Should be removed.

### 32. Macro redefinitions in `lgotwist.c`
- **File**: `lgotwist.c`
- **Lines**: 48–54
- **Severity**: Low
- **Description**: `#define mod`, `#define abs`, `#define min`, `#define max`
  — all conflict with standard library definitions. Should use standard
  library functions.

### 33. `register` keyword usage
- **File**: `LG.c`, `lgotwist.c`, `Rat.c`
- **Lines**: various (LG.c:70–103, 622; Rat.c throughout; lgotwist.c:70–103)
- **Severity**: Low
- **Description**: `register` keyword — obsolete in C++17 (ignored or error).
  Must be removed during migration.

### 34. Debug print in production code
- **File**: `Polynf.c`
- **Line**: 1497
- **Severity**: Low
- **Description**: `puts("PM");` — debug print left in production code. Should
  be removed or gated behind a debug flag.

### 35. Debug `printf` in production
- **File**: `MoriCone.cpp`
- **Line**: 949
- **Severity**: Low
- **Status**: Fixed (as part of TRACE_TRIANGULATION cleanup)
- **Description**: `printf` debug statement left in production code. Removed with the
  `#if TRACE_TRIANGULATION` cleanup in `MoriCone.cpp`.

### 36. Dead code blocks
- **File**: `MoriCone.c`
- **Lines**: 685, 860
- **Severity**: Low
- **Description**: `#ifdef OLD_code` and
  `#ifdef FIRST_TRY__TOO_COMPLICATED_BUT_MIGHT_BE_VIABLE` — dead code blocks
  that should be removed.

### 37. Duplicated rational arithmetic
- **File**: `lgotwist.c`
- **Lines**: 46–185
- **Severity**: Medium
- **Description**: The entire `Rat.h`/`Rat.c` rational arithmetic is duplicated
  inline in `lgotwist.c` — code duplication with divergent bug potential.
  Should use the shared `Rat.cpp` functions.

### 38. `#define TEST` / `#undef TEST` toggling mid-file
- **File**: `Polynf.c`
- **Lines**: 706–707, 1378–1681
- **Severity**: Low
- **Description**: `#define TEST` / `#undef TEST` / `#define TEST` toggling
  mid-file — very confusing, makes it hard to know which code paths are
  active. Should use `constexpr bool` flags at the top of the file.

### 39. Interactive `scanf` in library code
- **File**: `Nefpart.c`
- **Line**: 705
- **Severity**: Low
- **Description**: `scanf("%c", &c)` interactive debug pause in library code —
  breaks batch processing. Should be removed or gated behind `#ifndef NDEBUG`.

### 40. German comments
- **File**: `Polynf.c`
- **Lines**: 2754–3222
- **Severity**: Low
- **Description**: German comments in the `Make_Fano5d` section — mixed language
  in the codebase. Should be translated to English for consistency.

---

## Global state

### 41. Global `inFILE` / `outFILE`
- **File**: ALL driver files (`poly.c`, `cws.c`, `class.c`, `nef.c`, `mori.c`)
- **Severity**: High
- **Description**: `FILE *inFILE, *outFILE` globals — hidden dependencies
  across all library functions, prevents thread safety, functions mutate
  without restoring. Should be replaced with a `PalpContext` struct passed
  explicitly to functions that need I/O.

### 42. `inFILE`/`outFILE` mutated without restore
- **File**: `cws.c`
- **Line**: 5
- **Severity**: High
- **Description**: `inFILE`, `outFILE` are set globally and never restored. The
  filter mode (`inFILE = NULL` → `stdin`) vs file mode creates fragile state
  that persists across calls to `Read_CWS_PP` and similar functions.

### 43. `static EqList E` shared between functions
- **File**: `Subdb.c`
- **Lines**: 1029, 707
- **Severity**: Medium
- **Description**: `static EqList E` — shared mutable state between
  `DB_to_Hodge` and `PH_Sublat_Polys`. Functions interfere with each other's
  state.

### 44. `static int InputOK` state leak
- **File**: `Coord.c`
- **Lines**: 132, 228, 277
- **Severity**: Medium
- **Description**: `static int InputOK` — persistent across calls in
  `ReadCwsPp` / `Read_PP` / `Read_CWS`. Controls whether the `-h` help message
  is printed. State leaks between calls; the second call to `Read_*` behaves
  differently from the first.

---

## Algorithmic concerns (verify mathematical correctness)

### 45. `Inci64` limited to 64 points
- **File**: `MoriCone.cpp`
- **Lines**: ~614, ~668
- **Severity**: High
- **Status**: Fixed
- **Description**: `p <= 64` asserts in `Triang_from_SR` and `StanleyReisner`
  converted to hard errors (`fputs(..., stderr); exit(1);`) so the limit remains
  enforced even if asserts are disabled.

### 46. `binco` cap at 2999
- **File**: `MoriCone.c`
- **Lines**: 403, 421
- **Severity**: Medium
- **Description**: `binco` is capped at 2999. If the actual binomial
  coefficient exceeds 2999, `assert(++m < binco)` (see issue #1) fails. This
  limits the maximum number of simplices that can be handled.

### 47. `volatile` used to prevent optimization
- **File**: `Polynf.c`
- **Lines**: 348, 470
- **Severity**: Low
- **Description**: `volatile` used to prevent compiler optimization — suggests
  compiler-specific behavior issues. In C++, `volatile` has different
  semantics than in C. Should be replaced with `std::atomic` or removed with
  proper memory barriers, or the root cause of the optimization issue should
  be addressed.

### 48. Euler characteristic check via `assert`
- **File**: `Vertex.cpp`
- **Line**: ~260
- **Severity**: Medium
- **Status**: Fixed
- **Description**: `assert(M == 2*(d%2))` in face-incidence construction converted
  to a hard error: prints the F-vector, polytope, and face info, then `exit(1)`.