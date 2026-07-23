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
- **Status**: Closed
- **Description**: `Fgcd(Long a, Long b)` did `a %= b` without checking
  `b == 0`. Was fixed during an earlier migration pass; both `Fgcd` and `LFgcd`
  now guard with `if (!b) return a;` at the top, matching `NNgcd`/`LNNgcd`.

### 56. `lcm` / `Flcm` overflow risk
- **File**: `src/cws.cpp`
- **Lines**: 16, 377
- **Severity**: Medium
- **Status**: Closed
- **Description**: `lcm` and `Flcm` computed `(a * b) / NNgcd(a, b)`, which
  can overflow even when the result fits in `Long`. Was fixed during an earlier
  migration pass; both now compute `(a / g) * b`.

### 57. UB shift in `Inci64` operations
- **File**: `src/MoriCone.cpp`
- **Line**: 3055
- **Severity**: High
- **Status**: Closed
- **Description**: `T.I[i] &= ~(1 << (T.v - 1));` used an `int` literal `1`.
  Was fixed during an earlier migration pass; the current code reads
  `~(Inci64(1) << (T.v - 1))`.

### 58. Use of out-of-scope loop variable in `MakeRefWeights`
- **File**: `src/LG.cpp`
- **Line**: 986
- **Severity**: Medium
- **Status**: Closed
- **Description**: `puts(i ? "" : " 0");` used `i` after the `for` loop that
  declared it had ended. Was fixed during an earlier migration pass; the
  current code prints based on `P->n` directly:
  `if (P->n) puts(""); else puts(" 0");`.

### 59. EOF never detected in `lgotwist.cpp` skeleton reader
- **File**: `src/lgotwist.cpp`
- **Lines**: 1192–1197
- **Severity**: Medium
- **Status**: Closed
- **Description**: `char c; while ('\n' != (c = fgetc(palpContext.in)))
  if (c == EOF) {...}` stored `fgetc`'s `int` return in a `char`. Was fixed
  during an earlier migration pass; the current `ReadEOL()` uses `int c` and
  checks `c == EOF`.

### 60. `static Along totNF` leaks state across databases
- **File**: `src/Subdb.cpp`
- **Line**: 3257
- **Severity**: Medium
- **Status**: Closed
- **Description**: `static Along totNF;` in `Read_H_ucNF_from_DB` persisted
  across calls with different `DB` arguments. Was fixed during an earlier
  migration pass; the counter is now `DB->readHucNF_TotNF`, stored per database.

### 61. `int` overflow in `Ind * Ind * Ind`
- **File**: `src/Polynf.cpp`
- **Line**: 2447
- **Severity**: Medium
- **Status**: Closed
- **Description**: `int Ind = Divisibility_Index(P, V), I3 = Ind * Ind * Ind;`
  overflowed `int` for large `Ind`. Was fixed during an earlier migration
  pass; the current code computes `Long I3 = static_cast<Long>(Ind) * Ind * Ind;`.

### 62. Missing EOF checks on `fgetc` used as array index
- **File**: `src/Subdb.cpp`
- **Lines**: 1172, 1177, 2228, 2233, 3920, 3926
- **Severity**: Medium
- **Status**: Closed
- **Description**: `v = fgetc(F);` and `nu = fgetc(F);` were used as array
  indices without checking for EOF. Was fixed during an earlier migration
  pass; all four read sites now use `int ch = fgetc(F);` and check
  `if (ch == EOF)` before assigning to `v` or `nu`.

### 63. Syntax error in dead `TRIANG_CHECKSUM` block
- **File**: `src/MoriCone.cpp`
- **Line**: 849
- **Severity**: Low
- **Status**: Closed
- **Description**: `sumS + > T->I[i];` inside `#ifdef TRIANG_CHECKSUM` was a
  syntax error. Was fixed during an earlier migration pass; the block has been
  removed.

### 64. Dead static counter in `GL_Lattice_Basis`
- **File**: `src/Polynf.cpp`
- **Line**: 3865
- **Severity**: Low
- **Status**: Closed
- **Description**: `static int x;` was incremented (`x++`) but never read.
  Was fixed during an earlier migration pass; the counter has been removed.

### 65. Static mutable state in `Make_New_CEqs`
- **File**: `src/Vertex.cpp`
- **Lines**: 821–822
- **Severity**: Low
- **Status**: Closed
- **Description**: `static CEqList Bad_C; static INCI Bad_C_I[CEQ_Nmax];` were
  persistent mutable state. Was fixed during an earlier migration pass; they
  have been replaced with local/RAII storage.

### 66. Static `PolyPointList` in `FE_Close_the_Hole`
- **File**: `src/Subpoly.cpp`
- **Line**: 148
- **Severity**: Low
- **Status**: Closed
- **Description**: `static PolyPointList P;` was persistent mutable state.
  Was fixed during an earlier migration pass; it has been replaced with local
  allocation.

### 67. Static `Inci64 *I` in `TriList_to_MoriList`
- **File**: `src/MoriCone.cpp`
- **Line**: 2984
- **Severity**: Low
- **Status**: Closed
- **Description**: `static Inci64 *I = NULL;` owned a `malloc`'d buffer
  cached across calls. Was fixed during an earlier migration pass; it has been
  replaced with a `std::vector<Inci64>` owned by the caller.

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
- **Status**: Closed
- **Description**: `constexpr bool LLong_EEV = true; constexpr bool TEST_EEV
  = false;` were declared but never referenced. Removed; `EEV_To_Equation`
  continues to use the `LLong` path unconditionally.

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
- **Status**: Closed
- **Description**: `char c = fgetc(in);` stored the `int` return of `fgetc`
  in a `char`, losing the EOF sentinel. Fixed by changing all sites in
  `Coord.cpp`, `cws.cpp`, and `MoriCone.cpp` to `int c` and rewriting
  whitespace-skip loops to stop on EOF. `lgotwist.cpp` was already fixed in
  #59.

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
- **Status**: Closed
- **Description**: `#define BZangle(a, b) (ConeAngle(B.x[Z[a]], B.x[Z[b]]))`
  was the last function-like macro and implicitly captured `B` and `Z`.
  Replaced with `inline int BZangle(Matrix &B, const int *Z, int a, int b)`
  and updated all call sites.

### 79. `COEFF_Nmax` macro references local variables
- **File**: `src/LG.cpp`
- **Lines**: 26, 2338
- **Severity**: Low
- **Status**: Closed
- **Description**: `#define COEFF_Nmax (d * D + 2 * N)` was an object-like
  macro referencing local variables `d`, `D`, `N`. Replaced with
  `const int coeffNmax = d * D + 2 * N;` in `Calc_VaHo`; the duplicate
  top-level macro was also removed.

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
- **Status**: Closed
- **Description**: `fputc(fgetc(F), Fv);` copied bytes without checking
  for EOF. Replaced with a namespace helper `copyByteChecked(src, dst)`
  that aborts on EOF or write error. All four sites in `Polyi_2_DBo` and
  `Make_Hodge_db` were updated.

### 82. Missing `return` after `exit()` in stubs
- **File**: `src/MoriCone.cpp` (lines 108–109, 316–317)
- **Severity**: Low
- **Status**: Closed
- **Description**: The `INCI_LmR` stub referenced in the issue no longer
  exists; the only remaining stub in that area (`Print_Inci64_list`) is
  `void`-returning, so no return-after-`exit` warning applies. Closing as
  resolved by earlier cleanup.

---

## Complete 4D enumeration audit

The following findings come from tracing the documented 4D classification
workflow through seed generation, recursive enumeration, binary-list handling,
database merging, and completeness validation. Scale-dependent findings are
marked as such; the deterministic binary-reader regression was also reproduced
with a locally built `class-4d.x` and a newly generated two-dimensional list.

### 83. C++ evaluation order corrupts every binary list header read
- **File**: `src/Subadd.cpp:300`; `src/Subdb.cpp:116,1044,2096,3788`
- **Severity**: Critical
- **Status**: Open
- **Description**: Expressions of the form
  `NFnum[v][nu = fgetc(F)] = fgetUI(F)` rely on the left side being evaluated
  first. Since C++17, the assignment right operand is sequenced first, so
  `fgetUI()` consumes the `nu` byte and the first three count bytes, then
  `fgetc()` consumes the last count byte as the array index. All five readers
  therefore decode `[nu][uint count]` incorrectly. A list freshly written by
  `class-4d.x` failed its own `-C`, `-M`, and `-b` readers; for example, an
  encoded count of 1 was read as 16,777,216. This breaks recovery, exclusion
  lists, conversion, consistency checking, and database creation. Read `nu`
  and the count in separate checked statements before indexing.

### 84. Documented 4D command converts seeds instead of enumerating descendants
- **File**: `README.md:149-155`; `src/class.cpp:335-337,442-445,487-488`;
  `src/Subpoly.cpp:1041-1090,1093-1147`
- **Severity**: Critical
- **Status**: Open
- **Description**: The documented command
  `class-4d.x -a -pi cws4d.txt -po class4d.bin` treats `cws4d.txt` as an
  existing binary exclusion list because `-pi` is binary input. In addition,
  `-a` calls `Ascii_to_Binary`, which adds only each input polytope's normal
  form and never calls `Start_Make_All_Subpolys`. Thus the command either fails
  while parsing ASCII as binary or produces only seed polytopes, not their
  reflexive subpolytopes. The complete workflow must invoke normal
  classification with the ASCII seed stream as positional/stdin input.

### 85. Documented seed stream omits all simple 4D weight systems
- **File**: `README.md:131-151`; `src/cws.cpp:812-835,1083-1189,1609-1631,1744-1802`
- **Severity**: Critical
- **Status**: Open
- **Description**: The README writes `cws-4d.x -w4` to `weights.txt` but never
  consumes that file. Bare `-c4` emits combined systems constructed from the
  built-in lower-dimensional weight tables; it does not also emit the simple
  five-weight 4D systems generated by `-w4`. Classifying only `cws4d.txt`
  therefore omits whole descendant families. The production seed stream must
  contain both sets. The summary line emitted by `-w4` must also be removed or
  placed last because the current parser treats it as end of input (issue
  #86).

### 86. ASCII parsing silently truncates the seed stream
- **File**: `src/Coord.cpp:112-123,193-231,287-290,375-382`
- **Severity**: High
- **Status**: Open
- **Description**: `ReadCwsPp` skips only literal spaces before a record and
  returns the same value for clean EOF and an unexpected separator. After one
  valid record, a blank line, comment-only line, or leading non-space
  whitespace ends classification successfully and leaves all later seeds
  unprocessed. A valid final CWS or point-list record without a terminating
  newline is also discarded at EOF. Both behaviors were reproduced locally.
  Replace the character-at-a-time parser with a line-based parser that
  distinguishes record, clean EOF, blank/comment line, and malformed input.

### 87. Conversion and checking modes ignore streams opened by `class.cpp`
- **File**: `src/class.cpp:423-453`; `src/Subpoly.cpp:1041-1060`;
  `src/Subadd.cpp:2640-2661`; `include/palp/Subpoly.h:181-202`
- **Severity**: High
- **Status**: Open
- **Description**: `Ascii_to_Binary` and `Gen_Ascii_to_Binary` receive no
  input stream and call `Read_CWS_PP` with its default `stdin`, ignoring the
  positional input file already opened by `main`. Conversely, `Bin_2_ascii`
  and `Check_NF_Order` are called without the opened positional output stream
  and always write to `stdout`. A requested file can therefore be ignored or
  left empty, which is especially misleading for missing-mirror output. Thread
  explicit `FILE *in`/`FILE *out` arguments through every mode.

### 88. Missing recovery file leaves recursion depth uninitialized
- **File**: `src/Subadd.cpp:132-150,430-464`; `src/Subpoly.cpp:656-658,1012-1025`
- **Severity**: High
- **Status**: Open
- **Description**: With `-r`, `Init_NF_List` initializes `rd` only if the
  `.aux` file opens successfully. A missing file merely prints a message, so
  `Start_Make_All_Subpolys` copies an indeterminate `rd` into `rf`; recovery
  logic can then use uninitialized branch entries and skip recursive branches.
  Initialize `rd` unconditionally and make a requested but missing checkpoint
  a hard error.

### 89. Recovery checkpoints are not bound to the run they resume
- **File**: `README.md:190-192`; `src/Subadd.cpp:575-584,586-733`;
  `src/Subpoly.cpp:1012-1025,1093-1147`
- **Severity**: High
- **Status**: Open
- **Description**: Checkpoints contain list and branch state but no input hash,
  current seed identity, option set, format version, or `POLY_Dmax`. Supplying
  a reordered/changed seed stream, stale checkpoint, or differently built
  executable can apply saved branch indices to another polytope and silently
  omit branches. The README recovery command also omits the original ASCII
  input. Version checkpoints and verify the build, options, input identity,
  seed index, and seed normal form before resuming.

### 90. A rejected seed does not fail classification
- **File**: `src/Subpoly.cpp:647-663,1115-1147`
- **Severity**: Medium
- **Status**: Open
- **Description**: If `IP_Check` rejects a seed,
  `Start_Make_All_Subpolys` prints `IP_check negative!` and returns. The outer
  loop continues, writes a final list without that seed's descendants, and
  exits successfully. Propagate failure or record rejected seed identities and
  refuse to report successful completion when any input was skipped.

### 91. Missing-mirror output is mixed with diagnostics
- **File**: `README.md:157-167`; `src/Subdb.cpp:980-1065,1135-1289`
- **Severity**: High
- **Status**: Open
- **Description**: `-M` writes consistency banners, counts, and statistics to
  the same `stdout` stream as missing normal forms. Consequently the documented
  redirection cannot produce an empty `missing_mirrors.txt`, even when no
  mirror is missing, and downstream scripts cannot distinguish data from
  diagnostics. Send diagnostics to `stderr`, write records only to the selected
  output, and return a machine-readable missing count/status.

### 92. Mirror closure is documented as a completeness proof
- **File**: `README.md:162-169`
- **Severity**: Critical
- **Status**: Open
- **Description**: An empty missing-mirror set proves only mirror closure; any
  mirror-closed proper subset passes. Presenting emptiness as proof that no
  polytope was missed can certify a severely incomplete enumeration. Full 4D
  validation also needs an independently validated complete seed set, the exact
  total of 473,800,776 distinct normal forms, integrity checks, and preferably
  a canonical database manifest/hash.

### 93. Consistency checking ignores honest-list order and checksum failures
- **File**: `src/Subdb.cpp:921-945,1250-1289`
- **Severity**: High
- **Status**: Open
- **Description**: `Check_NF_Order` discards the return from
  `Check_hnf_order`, then unconditionally prints `order o.k.`. A mirror-flag
  checksum mismatch is only printed, and `-M` returns successfully immediately
  afterward. Deduplication and database lookup assume sorted lists, so accepting
  an unsorted/corrupt list can suppress valid new records. Treat either failure
  as fatal and print success only after all checks pass.

### 94. Recommended database command cannot create or complete a database
- **File**: `README.md:171-180`; `src/class.cpp:478-481`;
  `src/Subdb.cpp:431-450`
- **Severity**: High
- **Status**: Open
- **Description**: The documented
  `class-4d.x -di class4d_db -do class4d_db_complete` dispatches to
  `Add_Polya_2_DBi`, which immediately requires `-pa`; it neither copies nor
  completes a database. Document the actual binary-list-to-database and merge
  stages, including required inputs, and add round-trip validation.

### 95. `Add_Polya_2_DBi` rejects every usable input pair
- **File**: `src/Subdb.cpp:433,477-518`
- **Severity**: High
- **Status**: Open
- **Description**: While reading a database, `tNF` is never incremented, so
  every nonempty database fails the total-count check. If that is fixed, the
  auxiliary recursion-byte test rejects every nonempty file, including the
  required zero byte, because it tests for EOF rather than `rd_byte == 0`.
  Accumulate counts and read/validate the recursion byte once before calling
  `Read_Bin_Info`.

### 96. Database merge mangles output names and can dismantle its input
- **File**: `src/Subdb.cpp:443-455,640-685,807-827,868-870`
- **Severity**: High
- **Status**: Open
- **Description**: `Ofn` is a raw database prefix, but output paths use
  `replace(size - 4, 4, extension)`, replacing the last four characters of the
  prefix instead of appending `.vNN`, `.sl`, or `.info`. In in-place mode a
  source shard is renamed to a backup and then reopened under its now-missing
  original name; other source components are removed before a complete
  replacement database exists. Construct paths by appending extensions and
  always stage a complete out-of-place database before switching prefixes.

### 97. Binary payload readers turn truncation into valid-looking `0xff` data
- **File**: `src/Subadd.cpp:361-364,377-386,1999-2009`;
  `src/Subdb.cpp:313-321,408-410,2119,2431,3811`
- **Severity**: High
- **Status**: Open
- **Description**: Several payload loops assign `fgetc()` directly to
  `unsigned char`. Ordinary EOF becomes `0xff`, counters continue to advance,
  and `ferror()` remains false, so truncated normal forms can pass byte-count
  checks and poison sorted-list/database lookups. This is distinct from closed
  issues #62 and #81, which covered selected metadata and copy sites. Use a
  checked byte-read helper or exact-size `fread` for every record.

### 98. Large database-block arithmetic overflows before widening
- **File**: `src/Subadd.cpp:300-320`; `src/Subdb.cpp:112-130,175-180,226-321,795`;
  `include/palp/Subpoly.h:81-90,101-109`
- **Severity**: High (scale-dependent)
- **Status**: Open
- **Description**: Products such as `count * nu` are evaluated in 32-bit
  types before being added to `Along`; database conversion additionally stores
  block byte counts in `int`, and RAM sample positions use `int`/`UPint`.
  Multi-gigabyte `(vertex-count, encoded-length)` blocks can therefore wrap,
  underallocate buffers, skip copy loops, or corrupt the lookup index. The
  exact trigger depends on the full 4D block distribution. Use checked
  `uint64_t`/`size_t` arithmetic throughout and verify metadata against actual
  file sizes.

### 99. Fixed sublattice buffers can abort or overflow full-scale merges
- **File**: `include/palp/Subpoly.h:48-73`; `src/Subadd.cpp:114-125,1065-1091,2090-2122,2182-2237`;
  `src/Subdb.cpp:531-609,1365-1403,1554`
- **Severity**: High (scale-dependent)
- **Status**: Open
- **Description**: Classification hard-stops at `SL_Nmax == 65,536` for a 4D
  build. Merge/subtract paths use fixed `SLp[SL_Nmax]` arrays and byte pools
  sized by the average-size assumption `CperR_MAX == 32`, but increment/write
  union entries without checking either capacity. Large sublattice unions can
  thus abort enumeration or write out of bounds. Replace count and byte pools
  with vectors sized from checked file metadata and grow them on insertion.

### 100. Enumeration outputs are not committed transactionally
- **File**: `src/Subadd.cpp:744-800`; `src/Subdb.cpp:64-213`
- **Severity**: High
- **Status**: Open
- **Description**: Final lists and checkpoints are opened directly with
  `"wb"`, immediately destroying any previous good result. Database metadata
  is published before all shards. Most write paths check `ferror()` before
  `fclose()` but ignore close-time failures, so interruption, ENOSPC, or delayed
  I/O failure can leave a plausible-looking incomplete result. Write sibling
  temporary files, check write/flush/close, validate counts and sizes, publish
  database metadata last, and atomically rename only after success.

### 101. Merge/subtract output may truncate one of its inputs
- **File**: `src/Subadd.cpp:2090-2121`; `src/Subdb.cpp:1365-1403`
- **Severity**: High
- **Status**: Open
- **Description**: These operations open source files and then open the output
  with `"wb"` without rejecting the same path/inode. If output aliases an
  input, the source enumeration is truncated before it is read and can be lost
  even though the operation later fails. Reject aliasing by inode or always
  stage output to a distinct temporary file.

### 102. `Close_the_Hole` checks capacity after an out-of-bounds access
- **File**: `src/Subpoly.cpp:248-260`
- **Severity**: Critical (conditional memory-safety defect)
- **Status**: Open
- **Description**: When `_CEq->ne == CEQ_Nmax`, the function writes
  `CEq_INCI[_CEq->ne]` and passes it to `Irrel` before testing the bound. If the
  incidence is irrelevant, the check is skipped entirely. Reaching the fixed
  candidate-equation capacity therefore causes undefined behavior that can
  corrupt pruning state and omit a recursive branch. Check capacity before
  every access at `_CEq->ne`; use a dynamically sized container if the bound is
  not mathematically guaranteed.

### 103. Production 4D targets depend on the test configuration
- **File**: `CMakeLists.txt:12-14,99-149`; `README.md:127-129,184-186`
- **Severity**: Low
- **Status**: Open
- **Description**: Dimension-suffixed executables, including the recommended
  `class-4d.x`, are created only when `BUILD_TESTING` is enabled and only for
  `PALP_TEST_DIMENSIONS`. A tests-disabled build leaves only `class.x`, whose
  default is `POLY_Dmax=6`, despite the workflow requiring a 4D build. Separate
  production build dimensions from the test matrix, or document configuring
  the unsuffixed target with `-DPOLY_Dmax=4`; report the compiled limit at
  runtime.

### 104. No test exercises `class.x` or the classification file formats
- **File**: `tests/*.sh`; `CMakeLists.txt:163-175`
- **Severity**: Medium
- **Status**: Open
- **Description**: None of the test scripts invokes `class.x` or
  `class-4d.x`. There is no seed-to-classification smoke test, binary
  write/read round trip, recovery test, database conversion/merge test,
  malformed/truncated input test, or known-count completeness check. This gap
  allowed issue #83 to make a newly written list unreadable while all 196 tests
  remained green. Add small deterministic end-to-end fixtures to CI and a
  separately gated full 473,800,776-count/manifest validation.
