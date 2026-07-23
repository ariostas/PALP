# PALP Potential Issues / Bugs

Found during code inspection for the C → C++20 migration. Each item has a file,
line number, severity, and description. Fixes are applied in the appropriate
phase (see `PLAN.md`), one regression-tested fix at a time.

Line numbers refer to the current source files (they shift as migration
proceeds; function names and descriptions should remain locatable). Items
#1–#48 were found during the first C→C++17 audit; items #49+ were found during
the second C++20-target audit.

---

## Fixed (summary)

The following issues were resolved during migration/Phase 5. See the git log
and the detailed commit messages for the full context of each fix.

- **Critical correctness**
  - `#1` `assert(++m < binco)` side effect in `MoriCone.cpp`.
  - `#2` Misaligned store in `LG.cpp` `MakeMobius` (split packed buffer).
  - `#3` Memory leaks on error paths in `Find_Equations` / `IP_Check` / `Ref_Check` (`Vertex.cpp`).
  - `#6` `assert` guarding division by zero in `Polynf.cpp`.
  - `#45` `Inci64` limited to 64 points converted to hard error in `MoriCone.cpp`.
  - `#48` Euler characteristic check via `assert` in `Vertex.cpp`.

- **Error handling / robustness**
  - `#8` `exit(0)` on error paths replaced with `exit(1)`.
  - `#9` Memory leaks in `mori.cpp` main loop.
  - `#10` Huge struct never freed in `RgcWeights` (`cws.cpp`) — resolved by RAII.
  - `#11` Missing NULL check on `malloc` in `SingularInput.cpp` — replaced by `std::vector`.
  - `#12` Command injection via `system()` in `SingularInput.cpp` — replaced by `fork/execvp`.

- **Memory safety**
  - `#4` `realloc` losing old pointer in `E_Poly.cpp`.
  - `#5` Partial array initialization in `E_Poly.cpp`.
  - `#14` Memory leaks on error paths in `Fano5d` (`Polynf.cpp`).

- **Integer / type issues**
  - `#19` Integer overflow in volume computation (`Polynf.cpp`).
  - `#20` `lcm` macro overflow in `lgotwist.cpp`.
  - `#21` `printf` format truncation in `Rat.cpp`.
  - `#22` Platform-dependent integer sizes — `static_assert`s added in `palp_types.h`.
  - `#23` Unusual `int long nl` declaration in `nef.cpp`.

- **Input validation / buffer safety**
  - `#24` Off-by-one in `char c[999]` buffer in `LG.cpp`.
  - `#25` No bounds check on skeleton input in `lgotwist.cpp`.
  - `#26` `sprintf` buffer overflow risk in `cws.cpp`.
  - `#27` Tight `sprintf` buffer in `Subdb.cpp`.
  - `#28` Unchecked `fscanf` return values (all non-`nef.cpp` `Read_WPCICY` calls).

- **Code smells / maintainability**
  - `#29`–`#32` Local `min`/`max`/`mod`/`abs` macro redefinitions in `E_Poly.cpp`,
    `LG.h`, `Subpoly.h`, `lgotwist.c` — replaced with inline helpers or standard
    functions.
  - `#33` `register` keyword usage removed across `LG.cpp`, `lgotwist.cpp`, `Rat.cpp`.
  - `#34` Debug print in `Polynf.cpp` removed.
  - `#35` Debug `printf` in `MoriCone.cpp` removed with `TRACE_TRIANGULATION` cleanup.
  - `#36` Dead `#ifdef OLD_code` blocks in `MoriCone.cpp` removed.
  - `#37` Duplicated rational arithmetic in `lgotwist.cpp` removed; the file now
    uses `Rat`/`Rpr` from `Rat.cpp` and the global `PalpContext` for I/O.
  - `#38` `#define TEST` toggling in `Polynf.cpp` converted to `constexpr` flags.
  - `#39` Interactive `scanf` debug pause in `Nefpart.cpp` removed.
  - `#40` German comments in `Polynf.cpp` — still present, low priority.

- **Global state**
  - `#43` `static EqList E` in `Subdb.cpp` — still present, low priority.
  - `#44` `static int InputOK` in `Coord.cpp` — still present, low priority.

- **Algorithmic concerns**
  - `#46` `binco` cap at 2999 in `MoriCone.cpp` — architectural limit, unchanged.
  - `#47` `volatile` usage in `Polynf.cpp` — still present, low priority.

---

## Open / remaining

### 41. Global `inFILE` / `outFILE` — Fixed in Phase 5.11
- **File**: all driver files (`poly.cpp`, `cws.cpp`, `class.cpp`, `nef.cpp`, `mori.cpp`)
- **Severity**: High
- **Status**: Closed
- **Description**: `FILE *inFILE, *outFILE` globals created hidden dependencies
  and prevented thread safety. They were replaced with explicit `FILE *in` and
  `FILE *out` parameters threaded through all I/O helpers. Each driver `main()`
  now owns a local input/output file handle. This unblocks `lgotwist.cpp`
  sharing `Rat.cpp` in a follow-up step.

### 42. `inFILE`/`outFILE` mutated without restore — Fixed in Phase 5.11
- **File**: `cws.cpp`
- **Severity**: High
- **Status**: Closed
- **Description**: `inFILE`/`outFILE` were mutated globally inside `cws.cpp`
  helpers and never restored. The migration replaced the global state with
  explicit local `FILE *` parameters, eliminating the fragile filter-mode vs
  file-mode global state.

### 15. Static mutable state in `Subdb.cpp`
- **File**: `Subdb.cpp`
- **Line**: 1596 (original)
- **Severity**: Low
- **Status**: Open
- **Description**: `static unsigned char uc[NUC_Nmax]` and `static int ms3` —
  persistent mutable state, not thread-safe, not reentrant.

### 16. Static state in `LG.cpp`
- **File**: `LG.cpp`
- **Lines**: 280, 491 (original)
- **Severity**: Low
- **Status**: Closed
- **Description**: `static int MaxPoNum` in `TEST_WeightMakePoints` and the
  `static int M` progress counter in `Add_Mono_2_Poly` were converted to
  translation-unit anonymous-namespace variables. The persistent single-threaded
  state is now explicit; full reentrance would require passing counters through
  call sites, which is out of scope for this cleanup.

### 17. Static counters in `Polynf.cpp`
- **File**: `Polynf.cpp`
- **Lines**: 1168, 1947, 2635 (original)
- **Severity**: Low
- **Status**: Open
- **Description**: `static int` counters in `ConifoldSing`, `GL_Lattice_Basis`,
  and `Fano5d` — not thread-safe, state persists across calls.

### 18. Global mutable pointer in `Subadd.cpp`
- **File**: `Subadd.cpp`
- **Line**: 96 (original)
- **Severity**: Low
- **Status**: Closed
- **Description**: `NF_List *AuxNFLptr = NULL` was referenced in the original
  C source but was already removed during the migration. No further action was
  required; closing as resolved.

### 40. German comments in `Polynf.cpp`
- **File**: `Polynf.cpp`
- **Lines**: 2754–3222 (original)
- **Severity**: Low
- **Status**: Open
- **Description**: German comments in the `Make_Fano5d` section — mixed language.

### 46. `binco` cap at 2999
- **File**: `MoriCone.cpp`
- **Lines**: 403, 421 (original)
- **Severity**: Medium
- **Status**: Open (architectural limit)
- **Description**: `binco` is capped at 2999. This limits the maximum number of
  simplices that can be handled. Raising it would require a larger review of
  dependent buffers.

### 47. `volatile` used to prevent optimization
- **File**: `Subdb.cpp`
- **Lines**: 1640 (original)
- **Severity**: Low
- **Status**: Closed
- **Description**: A single `volatile unsigned int *_nf` pointer was used in
  `Subtract_Aux_from_DB` to force the compiler to re-read `FIo.NFnum[v][nu]`
  through a pointer. With no concurrent access in this single-threaded code,
  the `volatile` qualifier was unnecessary; removed and replaced with a plain
  `unsigned int *`.

---

## Second-audit findings (C++20 target)

### 49. Wrong index in `QComplete_Poly`
- **File**: `src/Vertex.cpp`
- **Line**: 1591
- **Severity**: High
- **Status**: Closed
- **Description**: `EyD[j] = _E->e[i].c * Den;` used `i` (the loop variable
  of the preceding `for` loop, which holds `n` at this point) instead of `j`
  (the current outer-loop index). Was fixed during an earlier migration pass;
  the current code reads `_E->e[j].c * Den`.

### 50. Self-swap no-op in `Calc_VaHo`
- **File**: `src/LG.cpp`
- **Line**: 2427
- **Severity**: High
- **Status**: Closed
- **Description**: `if (w[j] < w[i]) swap(&w[j], &w[j]);` swapped `w[j]`
  with itself. Was fixed during an earlier migration pass; the current code
  reads `swap(&w[i], &w[j])`.

### 51. Wrong loop index in `Read_WZeight`
- **File**: `src/LG.cpp`
- **Lines**: 250–251
- **Severity**: High
- **Status**: Closed
- **Description**: `for (k = 0; k < a; k++) Za[j] /= g;` divided `Za[j]`
  (where `j == d`) repeatedly instead of `Za[k]`. Was fixed during an earlier
  migration pass; the current code reads `Za[k] /= g`.

### 52. Division-by-zero check on wrong element in `Compute_X0`
- **File**: `src/Coord.cpp`
- **Line**: 977
- **Severity**: Medium
- **Status**: Closed
- **Description**: `if (_C->W[j][0]) { if (_C->d[j] % _C->W[j][N]) return 0; }`
  checked `W[j][0]` for zero but divided by `W[j][N]`. Was fixed during an
  earlier migration pass; the current code checks the same element it divides
  by in both the base case (`W[j][0]`) and the recursive case (`W[j][N]`).

### 53. XOR instead of multiplication in `prideco`
- **File**: `src/lgotwist.cpp`
- **Line**: 212
- **Severity**: High
- **Status**: Closed
- **Description**: `(p = prime[n++]) ^ 2` used bitwise XOR (`^`) instead of
  `p * p`. Was fixed during an earlier migration pass; the current code reads
  `(p = prime[n++]) * p`.

### 54. Wrong indices in `Make_211_CWS` comparison
- **File**: `src/cws.cpp`
- **Line**: 2206
- **Severity**: Medium
- **Status**: Closed
- **Description**: `if ((W[0].w[0] != W[0].w[1]) && (W[1].w[0] != W[2].w[0]))`
  compared `W[1].w[0]` with `W[2].w[0]` instead of `W[1].w[0]` with `W[1].w[1]`.
  Was fixed during an earlier migration pass; the current code reads
  `(W[1].w[0] != W[1].w[1])`.

### 55. `Fgcd` / `LFgcd` divide by zero
- **File**: `src/Rat.cpp`
- **Lines**: 105–108, 355
- **Severity**: High
- **Status**: Open
- **Description**: `Fgcd(Long a, Long b)` does `a %= b` without checking
  `b == 0`. If `b` is zero, this is undefined behavior (integer division by
  zero). Same bug in `LFgcd` (line 355). The other gcd variants
  (`NNgcd:116`, `LNNgcd:366`) correctly guard with `if (!b) return a;`.

### 56. `lcm` / `Flcm` overflow risk
- **File**: `src/cws.cpp`
- **Lines**: 16, 377
- **Severity**: Medium
- **Status**: Open
- **Description**: `inline Long lcm(Long a, Long b) { return (a * b) /
  NNgcd(a, b); }` computes `a * b` before dividing, which can overflow even
  when `lcm(a, b)` fits in `Long`. Should be `(a / NNgcd(a, b)) * b`. Same
  in `Flcm` (line 377). The `Lcm` in `lgotwist.cpp:55` already does this
  correctly.

### 57. UB shift in `Inci64` operations
- **File**: `src/MoriCone.cpp`
- **Line**: 3055
- **Severity**: High
- **Status**: Open
- **Description**: `T.I[i] &= ~(1 << (T.v - 1));` uses `int` literal `1`.
  When `T.v - 1 >= 32` (i.e., `VERT_Nmax > 32`, which holds for
  `POLY_Dmax >= 5` where `VERT_Nmax == 64`), this is undefined behavior
  (shift by >= width of int). Should be `Inci64(1) << (T.v - 1)`.

### 58. Use of out-of-scope loop variable in `MakeRefWeights`
- **File**: `src/LG.cpp`
- **Line**: 986
- **Severity**: Medium
- **Status**: Open
- **Description**: `puts(i ? "" : " 0");` uses `i` after the `for` loop that
  declared it has ended. In C++ this is ill-formed; in practice `i` holds
  `P->n`, so the expression is always truthy and the function always prints
  `""`, never `" 0"` for empty polynomials. Likely intended to print `" 0"`
  when `P->n == 0`.

### 59. EOF never detected in `lgotwist.cpp` skeleton reader
- **File**: `src/lgotwist.cpp`
- **Lines**: 1192–1197
- **Severity**: Medium
- **Status**: Open
- **Description**: `char c; while ('\n' != (c = fgetc(palpContext.in)))
  if (c == EOF) {...}` — `c` is `char`, but `fgetc` returns `int` and
  `EOF` is `-1`. If `char` is unsigned, `c == EOF` is never true; if
  signed, `c` becomes `-1` (not `'\n'` or `' '`), so the loop reads past
  EOF. Should use `int c`.

### 60. `static Along totNF` leaks state across databases
- **File**: `src/Subdb.cpp`
- **Line**: 3257
- **Severity**: Medium
- **Status**: Open
- **Description**: `static Along totNF;` in `Read_H_ucNF_from_DB` persists
  across calls with different `DB` arguments. Calling with DB1 then DB2
  compares DB2's totals against DB1's leftover `totNF`, producing wrong
  results. Should be an anonymous-namespace variable or threaded through
  the call site.

### 61. `int` overflow in `Ind * Ind * Ind`
- **File**: `src/Polynf.cpp`
- **Line**: 2447
- **Severity**: Medium
- **Status**: Open
- **Description**: `int Ind = Divisibility_Index(P, V), I3 = Ind * Ind * Ind;`
  — `Ind * Ind * Ind` overflows `int` if `Ind` is large (e.g., `Ind = 1626`
  gives `I3 ≈ 4.3e9 > INT_MAX`). Should use `Long`.

### 62. Missing EOF checks on `fgetc` used as array index
- **File**: `src/Subdb.cpp`
- **Lines**: 1172, 1177, 2228, 2233, 3920, 3926
- **Severity**: Medium
- **Status**: Open
- **Description**: `v = fgetc(F);` and `nu = fgetc(F);` are used as array
  indices (`L.nNUC[v]`, `DB->Fv[v]`, etc.) without checking for EOF. If
  `fgetc` returns `-1`, the subsequent array access is out-of-bounds.

### 63. Syntax error in dead `TRIANG_CHECKSUM` block
- **File**: `src/MoriCone.cpp`
- **Line**: 849
- **Severity**: Low
- **Status**: Open
- **Description**: `sumS + > T->I[i];` inside `#ifdef TRIANG_CHECKSUM`
  (never defined) is a syntax error. The block is dead but would fail to
  compile if the macro were ever enabled. Remove the block.

### 64. Dead static counter in `GL_Lattice_Basis`
- **File**: `src/Polynf.cpp`
- **Line**: 3865
- **Severity**: Low
- **Status**: Open
- **Description**: `static int x;` is incremented (`x++` at line 3866) but
  never read. Dead code; remove.

### 65. Static mutable state in `Make_New_CEqs`
- **File**: `src/Vertex.cpp`
- **Lines**: 821–822
- **Severity**: Low
- **Status**: Open
- **Description**: `static CEqList Bad_C; static INCI Bad_C_I[CEQ_Nmax];` —
  persistent mutable state, not thread-safe, not reentrant. Move to
  anonymous namespace.

### 66. Static `PolyPointList` in `FE_Close_the_Hole`
- **File**: `src/Subpoly.cpp`
- **Line**: 148
- **Severity**: Low
- **Status**: Open
- **Description**: `static PolyPointList P;` — persistent mutable state,
  not thread-safe. Move to anonymous namespace or allocate locally.

### 67. Static `Inci64 *I` in `TriList_to_MoriList`
- **File**: `src/MoriCone.cpp`
- **Line**: 2984
- **Severity**: Low
- **Status**: Open
- **Description**: `static Inci64 *I = NULL;` owns a `malloc`'d buffer
  cached across calls. Not reentrant; prevents RAII conversion. Move to
  `std::vector<Inci64>` owned by the caller.

### 68. Large commented-out code blocks
- **File**: multiple (`MoriCone.cpp`, `Polynf.cpp`, `Subadd.cpp`,
  `Subdb.cpp`, `lgotwist.cpp`, `cws.cpp`)
- **Severity**: Low
- **Status**: Open
- **Description**: ~80+ lines of commented-out diagnostic code in
  `MoriCone.cpp` (lines 844–853, 1318–1329, 1664–1692, 3062–3071,
  3083–3091, 3097–3105), plus smaller blocks in `Polynf.cpp:2692–2703`,
  `Subadd.cpp:2399–2402`, `Subdb.cpp:801–804`,
  `lgotwist.cpp:1313–1322`, `cws.cpp:2453–2458`. Remove.

### 69. Dead configuration flags `LLong_EEV` / `TEST_EEV`
- **File**: `src/Vertex.cpp`
- **Lines**: 323–326
- **Severity**: Low
- **Status**: Open
- **Description**: `constexpr bool LLong_EEV = true; constexpr bool TEST_EEV
  = false;` are declared but never referenced. `EEV_To_Equation` always
  uses the `LLong` path unconditionally. Remove.

### 70. ~1,000+ `exit(1)` calls in library functions
- **File**: all `src/*.cpp` except driver `main()` functions
- **Severity**: Medium (architectural)
- **Status**: Open
- **Description**: `exit(1)` is called in over 1,000 places in non-`main`
  functions, making PALP unusable as an embeddable library. The largest
  offenders: `Subdb.cpp` (227), `Polynf.cpp` (169), `MoriCone.cpp` (121),
  `LG.cpp` (110), `Subadd.cpp` (92), `cws.cpp` (68), `Coord.cpp` (56).
  Should be converted to `throw std::runtime_error` (caught in `main()`),
  via a `palp::die()` helper.

### 71. `fgetc` return stored in `char` (EOF handling)
- **File**: `src/Coord.cpp` (lines 21–29, 122, 219, 406, 508),
  `src/cws.cpp` (line 599), `src/lgotwist.cpp` (lines 1192–1197),
  `src/MoriCone.cpp` (lines 2517–2520, 2541–2551)
- **Severity**: Medium
- **Status**: Open
- **Description**: `char c = fgetc(in);` stores the `int` return of
  `fgetc` in a `char`, losing the EOF sentinel. If `char` is unsigned,
  EOF becomes 255 and is indistinguishable from a valid byte; if signed,
  EOF becomes -1 but is not explicitly checked. Should use `int c`.

### 72. VLA usage (non-standard C++)
- **File**: `src/MoriCone.cpp` (lines 2004, 2006),
  `src/SingularInput.cpp` (line 195)
- **Severity**: Medium
- **Status**: Open
- **Description**: `Inci64 T[naT]`, `int nt[ANfan][ANtri]`, and
  `int DegreeVec[dim]` are variable-length arrays, a non-standard C++
  extension. Replace with `std::vector`.

### 73. `qsort` with C-style comparator
- **File**: `src/Vertex.cpp`
- **Line**: 289
- **Severity**: Low
- **Status**: Open
- **Description**: `int diff(const void *a, const void *b)` with
  `*((int *)a) - *((int *)b)` is a `qsort` comparator using C-style casts.
  Replace with `std::sort` and a lambda.

### 74. 212 C-style casts
- **File**: all `src/*.cpp`
- **Severity**: Low
- **Status**: Open
- **Description**: 212 C-style casts remain (`(int)`, `(char *)`, `(Long)`,
  etc.). The bulk are `(int)` casts in `printf` format args that vanish with
  C++20 `<format>`. Replace with `static_cast`/`const_cast`/
  `reinterpret_cast` or eliminate by changing types.

### 75. 91 `typedef` statements (should be `using`)
- **File**: `include/palp/*.h`, `src/*.cpp`
- **Severity**: Low
- **Status**: Open
- **Description**: 91 `typedef` statements remain, including 50+
  `typedef struct { ... } Name;` (should be `struct Name { ... };`) and
  array typedefs (should be `using Name = std::array<T, N>;`).

### 76. Duplicate type definitions
- **File**: multiple
- **Severity**: Low
- **Status**: Open
- **Description**: `Matrix` (Mori.h:157 and Polynf.cpp:1706 — identical
  bodies, different TUs), `AmbiPointList` (LG.h:51 forward decl vs
  Nef.h:66 full def), `symlist` (LG.cpp:2546 struct vs lgotwist.cpp:96
  `long[NS]` — incompatible meanings), `VPerm`/`VPermList`
  (Polynf.cpp:1373 and poly.cpp:85), `CWLatticeBasis` (Coord.cpp:9 and
  nef.cpp:24), `subl_int` (Subpoly.h:66 vs Subpoly.cpp:8). Consolidate
  into shared headers.

### 77. `lgotwist.cpp` file-scope mutable globals
- **File**: `src/lgotwist.cpp`
- **Lines**: 32, 105–129, 143, 175, 621–622, 1147
- **Severity**: Medium
- **Status**: Open
- **Description**: ~20 file-scope mutable globals (`n`, `d`, `det`, `wei`,
  `prdet`, `lcmd`, `specnum`, etc.) constitute the entire application
  state as plain globals, preventing reentrance and thread safety. At
  minimum, `mask[]` (line 32) and `prime[]` (line 175) should be
  `constexpr`; the rest should be consolidated into the `LgoTwistState`
  struct.

### 78. Last function-like macro `BZangle`
- **File**: `src/MoriCone.cpp`
- **Line**: 962
- **Severity**: Low
- **Status**: Open
- **Description**: `#define BZangle(a, b) (ConeAngle(B.x[Z[a]], B.x[Z[b]]))`
  is the last remaining function-like macro (the others were converted in
  Phase 4). It captures `B` and `Z` implicitly. Should become an inline
  function taking `Matrix &B, const int *Z, int a, int b`.

### 79. `COEFF_Nmax` macro references local variables
- **File**: `src/LG.cpp`
- **Lines**: 26, 2338
- **Severity**: Low
- **Status**: Open
- **Description**: `#define COEFF_Nmax (d * D + 2 * N)` is an object-like
  macro that expands to an expression referencing local variables `d`,
  `D`, `N` from the enclosing function. It is defined twice (lines 26 and
  2338) in different functions with the same formula. Should be a
  `const int` local or a small inline helper.

### 80. Magic numbers in physics formulas
- **File**: `src/Polynf.cpp` (lines 641–647, 1432, 2285, 2428),
  `src/E_Poly.cpp` (line 634)
- **Severity**: Low
- **Status**: Open
- **Description**: Hardcoded coefficients like `48 + 6 * (ho[1] - ho[2] +
  ho[3])` (Euler characteristic), `24 * (ho[1] - ho[2] + ho[3] - ho[4])`,
  `44 * h[0][0] + 4 * h[1][1] - 2 * h[1][2] + ...` are from physics
  formulas but uncommented. Add named `constexpr` constants or comments
  citing the source.

### 81. Unchecked `fputc(fgetc(F), Fv)` in `Subdb.cpp`
- **File**: `src/Subdb.cpp`
- **Lines**: 168, 186, 2788, 2790–2797, 3041–3051
- **Severity**: Medium
- **Status**: Open
- **Description**: `fputc(fgetc(F), Fv);` copies bytes without checking
  for EOF. If `fgetc` returns EOF, `fputc(EOF, ...)` writes 255 (or -1
  cast to unsigned char), corrupting the output file.

### 82. Missing `return` after `exit()` in stubs
- **File**: `src/MoriCone.cpp` (lines 108–109, 316–317)
- **Severity**: Low
- **Status**: Open
- **Description**: `int INCI_LmR(...) { puts("Implement INCI_LmR");
  exit(1); }` ends without `return` after `exit(1)`. Technically UB if
  `exit` didn't terminate (it does, but compilers warn). Add `return 0;`
  or mark `[[noreturn]]`.
