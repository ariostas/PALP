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
- **File**: `Vertex.c`
- **Lines**: 633–654, 667–683, 696–710
- **Severity**: Critical
- **Description**: These three functions `malloc` `CEq`, `CEq_I`, `F_I`. On
  certain paths (e.g., `GLZ_Start_Simplex` returning nonzero), they `free` all
  three and return. But on the `exit(0)` path inside `Finish_Find_Equations`
  and subroutines, the memory is leaked. The three functions also have
  inconsistent error handling: `Find_Equations` returns 0 on simplex-codim > 0
  but `IP_Check` returns 0 and leaks in an analogous path.

### 4. `realloc` loses old pointer on failure
- **File**: `E_Poly.c`
- **Line**: 53
- **Severity**: High
- **Description**: `DYNadd_for_completion` assigns `realloc` result directly
  to `_CP->L` without saving the old pointer. If `realloc` returns NULL, the
  old memory is leaked. Standard C bug. Fix: save old pointer, check for NULL,
  free old on failure.

### 5. Partial array initialization
- **File**: `E_Poly.c`
- **Line**: 673
- **Severity**: Medium
- **Description**: `int h[POLY_Dmax][POLY_Dmax] = {{0},{0}}` — only initializes
  the first two rows explicitly. While C99 does zero-fill the rest, this is
  misleading and fragile on some compilers. Should be `= {}` or `= {{0}}` for
  full zeroing.

### 6. `assert` guarding division by zero
- **File**: `Polynf.c`
- **Line**: 894
- **Severity**: High
- **Description**: `assert(g > 0)` is followed by division by `g`. If assert
  is disabled (`NDEBUG`), division by zero occurs. Should be an explicit
  `if (g <= 0) { /* error */ }` check, not an assert.

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
- **File**: `mori.c`
- **Lines**: 78–88
- **Severity**: High
- **Description**: All allocations (`CW`, `E`, `DE`, `_P`, `_DP`, `PM`, `DPM`)
  are never freed — memory leak on every input iteration of the `while` loop.
  The program returns 0 at the end without freeing anything.

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
- **Description**: `he = malloc(divclassnr * sizeof(int))` — no NULL check
  on `malloc` return. Subsequent use of `he` could segfault.

### 12. Command injection via `system()`
- **File**: `SingularInput.c`
- **Line**: 584
- **Severity**: Medium
- **Description**: `system(SingularCall)` — the command string incorporates
  `getenv("TMPDIR")` which could contain shell metacharacters. Also, the
  return value only prints a message; the error is not propagated. Should
  sanitize the temp directory path or use `execve`-style calls.

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
- **File**: `Polynf.c`
- **Lines**: 2632–2735
- **Severity**: Medium
- **Description**: `Fano5d` allocates Q, F, M, G — but frees them only on the
  success path (line 2735). All error/early-return paths leak all four
  allocations.

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
- **File**: `Polynf.c`
- **Lines**: 855, 861
- **Severity**: Medium
- **Description**: `I *= VM[i][i]` — product of diagonal entries can overflow
  `Long` for large polytopes. No overflow check. Could produce wrong simplex
  volumes, affecting fibration analysis.

### 20. `lcm` macro overflow
- **File**: `lgotwist.c`
- **Line**: 52
- **Severity**: Medium
- **Description**: `#define lcm(a,b) ((a)*(b)/gcd((a),(b)))` — the product
  `a*b` is computed before the division by `gcd`, risking overflow even when
  the result would fit. Should divide first: `((a)/gcd(a,b))*(b)`.

### 21. `printf` format truncation
- **File**: `Rat.c`
- **Line**: 59
- **Severity**: Low
- **Description**: `Rpr` uses `%d` format with `(int) c.N` — truncation if
  `Long` exceeds `int` range. Should use `%ld` or a C++ stream.

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
- **Description**: `int long nl` — unusual declaration order (legal C but
  confusing). Should be `long int nl` or just `long nl`.

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
- **File**: `MoriCone.c`
- **Line**: 949
- **Severity**: Low
- **Description**: `printf` debug statement left in production code.

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
- **File**: `MoriCone.c`
- **Lines**: 402, 418
- **Severity**: High
- **Description**: `p <= 64` assert — `Inci64` (unsigned long long) limited to
  64 points. If a polytope has more than 64 points, silent failure or wrong
  triangulation results. The assert is the only protection; with `NDEBUG`, it
  disappears entirely.

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
- **File**: `Vertex.c`
- **Line**: 154
- **Severity**: Medium
- **Description**: `assert(M == 2*(d%2))` — if this assert is disabled, incorrect
  face counts go undetected, producing wrong Hodge numbers silently. This is
  a mathematical correctness check that should be a hard error, not an assert.