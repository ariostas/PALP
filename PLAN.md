#PALP Modernization Plan : C → C++ 20

## Goal

Migrate PALP from C to modern C++20 while preserving byte-for-byte identical CLI
output and behavior. The target is idiomatic, RAII-based, safe C++20 using
`std::array`/`std::vector`/`std::span`/`std::unique_ptr`, LLVM naming
conventions, and `<cXXX>` C++ headers exclusively. Bugs found during the audit
are tracked in `ISSUES.md` and fixed incrementally.

## Principles

- **Behavior preservation first**: every step must keep the full 196-test
  CTest suite passing (Release + ASan).
- **One file or one concern per step**: convert a single source file or a
  single cross-cutting concern at a time. Commit after every verified step.
- **Run `clang-format -i -style=LLVM` before every commit.**
- **`POLY_Dmax` stays compile-time** (macro); no dynamic dimensioning in the
  migration phase.
- **No new features**: modernize language constructs only, not algorithms.
- **Exact output compatibility**: every `tests/*.sh` output must match the
  pre-migration baseline byte-for-byte.
- **No `exit()` in library code**: new error paths use exceptions
  (`throw std::runtime_error`) or return codes; `exit()` is confined to
  driver `main()` functions.
- **No manual memory management**: `malloc`/`free`/`new`/`delete` are replaced
  by RAII (`std::vector`, `std::unique_ptr`, `std::string`, stack objects).
- **No legacy C headers**: `<stdio.h>` → `<cstdio>`, etc. POSIX-only headers
  (`<unistd.h>`, `<sys/*.h>`) are wrapped behind a small facade.
- **No C-style casts**: use `static_cast`/`reinterpret_cast`/`const_cast`, or
  better, eliminate the cast by changing types.
- **No `NULL`**: use `nullptr`.
- **No `typedef`**: use `using` aliases; drop the redundant `typedef struct`
  tag and use plain `struct Name { ... };`.
- **No function-local `static` mutable state**: move to anonymous-namespace
  variables (already done for `LG.cpp`, `Subadd.cpp`) or thread through call
  sites.
- **No `goto`**: refactor to structured control flow.
- **Naming**: LLVM style — `camelCase` for functions/methods/variables,
  `PascalCase` for classes/structs/aliases, `UPPER_CASE` for macros and
  `constexpr` compile-time dimension limits that drive `#if`/array bounds.
  Member fields are `camelCase` (no leading underscore); parameters are
  `camelCase` (no `_` prefix).

## Test harness

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build

# Sanitizer builds
cmake -S . -B build/asan -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-O1 -g -fno-omit-frame-pointer -fsanitize=address" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address"
cmake --build build/asan && ctest --test-dir build/asan

cmake -S . -B build/ubsan -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-O1 -g -fno-omit-frame-pointer -fsanitize=undefined -fno-sanitize-recover=undefined" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=undefined"
cmake --build build/ubsan && ctest --test-dir build/ubsan
```

---

## Completed migration work

### Phase 0 — Infrastructure
- [x] `CMakeLists.txt`: CXX-only, C++17, `-Wall -Wextra`, `palp_sources()` helper.
  Headers moved to `include/palp/` with `#pragma once`; include paths updated.
  `palp_types.h` shim added and integrated via `Global.h`.

### Phase 1 — Foundation libraries (C→C++17 file-by-file)
- [x] `Rat.cpp`, `Vertex.cpp`, `Coord.cpp`, `Polynf.cpp`, `LG.cpp`.

### Phase 2 — Driver programs
- [x] `poly.cpp`, `cws.cpp`, `class.cpp`, `nef.cpp`, `mori.cpp`.

### Phase 3 — Secondary libraries
- [x] `E_Poly.cpp`, `Nefpart.cpp`, `MoriCone.cpp`, `SingularInput.cpp`,
  `Subdb.cpp`, `Subpoly.cpp`, `Subadd.cpp`.

### Phase 4 — First modernization pass
- [x] Removed `.c` files from CMake, deleted legacy Makefiles, moved sources
  to `src/`, removed `extern "C"` shims, consolidated `palp::min`/`palp::max`.
- [x] `#define` → `constexpr`/`using` sweep where safe (kept `POLY_Dmax` and
  array-size-driving macros as macros).
- [x] Fixed-size local buffers → `std::array`/`std::vector` where applicable.
- [x] `lgotwist.cpp` shares `Rat.cpp` helpers and `PalpContext`.
- [x] Removed `-DNDEBUG` workaround; Release + ASAN + UBSAN pass 196/196 for
  `POLY_Dmax` = 4, 5, 6, 11.
- [x] Converted safe local `#define` constants/helpers to `constexpr` /
  inline templates (`Polynf.cpp`, `Subadd.cpp`, `MoriCone.cpp`).
- [x] Replaced function-local `static` mutable state in `LG.cpp` and
  `Subadd.cpp` with anonymous-namespace variables.
- [x] Removed unnecessary `volatile` qualifier (`Subdb.cpp`, issue #47).

### Phase 5 — First bug-fix pass
- [x] **5.1** `assert(++m < binco)` side effect moved out of `assert()`.
- [x] **5.2** `exit(0)` on error paths → `exit(1)`.
- [x] **5.3/5.3a** Memory leaks in `Fano5d` and `mori.cpp` fixed via RAII.
- [x] **5.4/5.5** `realloc`/partial-init issues in `E_Poly.cpp` fixed via
  `std::vector`.
- [x] **5.6** `fscanf` return-value checks added (except `nef.cpp`
  `Read_WPCICY`).
- [x] **5.7** `system()` injection in `SingularInput.cpp` → `fork/execvp`.
- [x] **5.8** Buffer-overflow risks fixed (`LG.cpp`, `lgotwist.cpp`,
  `cws.cpp`, `Subdb.cpp`).
- [x] **5.9** Integer-overflow risks fixed (`Polynf.cpp` volume, `lgotwist.cpp`
  `Lcm`).
- [x] **5.10/5.13/5.14** Critical `assert()`s → explicit errors.
- [x] **5.11** Global `inFILE`/`outFILE` → explicit `FILE *in`/`FILE *out`
  parameters + `PalpContext` (5.11.1–5.11.5).
- [x] **5.12** Dead code removed (`FileRW`, `OLD_code`, debug prints,
  `TEST` toggles).

---

## Remaining work (C++20 target)

A full audit (see `AUDIT.md` for the raw findings) identified the following
categories of remaining work. Items are ordered so that each step is
independently testable and keeps the suite green.

### Phase 6 — C++20 baseline and headers

- [x] **6.1** Bump `CMAKE_CXX_STANDARD` 17 → 20 in `CMakeLists.txt:17`.
  Verify all 196 tests still pass on the compilers of interest.
- [x] **6.2** Replace legacy C headers in `include/palp/Global.h` (lines 4–8)
  and `include/palp/Subpoly.h` (line 3) with their `<cXXX>` equivalents:
  `<cassert>`, `<cstdio>`, `<cstdlib>`, `<cstring>`, `<ctime>`, `<climits>`.
- [x] **6.3** Wrap POSIX-only headers in `src/SingularInput.cpp` (lines 14–17:
  `<fcntl.h>`, `<sys/types.h>`, `<sys/wait.h>`, `<unistd.h>`) behind a small
  `palp/PosixProcess.h` facade; POSIX APIs have no standard C++ equivalent.
- [x] **6.4** Mechanically replace `NULL` → `nullptr` across `src/` (~214
  occurrences). `time(NULL)` is left as `time(nullptr)` (it remains a C-style
  seed and is acceptable with `<ctime>`).
- [x] **6.5** Convert `#ifdef __cplusplus` `Long`/`LLong` macros in
  `Global.h:32-33` to be C++-only (the C fallback is dead); keep the `using`
  aliases in `palp_types.h:14-15` as the single source of truth.

### Phase 7 — Critical bug fixes (from audit)

These are confirmed logic bugs that produce wrong results. Fix one per commit,
each verified by the full test suite.

- [x] **7.1** `src/Vertex.cpp:1591` — `EyD[j] = _E->e[i].c * Den;` uses wrong
  index `i` (should be `j`). (ISSUES.md #49)
- [x] **7.2** `src/LG.cpp:2427` — `swap(&w[j], &w[j]);` swaps element with
  itself; sort is a no-op. Should be `swap(&w[i], &w[j])`. (ISSUES.md #50)
- [x] **7.3** `src/LG.cpp:250-251` — `Za[j] /= g;` uses `j` instead of `k` in
  the loop. (ISSUES.md #51)
- [x] **7.4** `src/Coord.cpp:977` — checks `W[j][0]` but divides by
  `W[j][N]`; potential division by zero. (ISSUES.md #52)
- [x] **7.5** `src/lgotwist.cpp:212` — `(p = prime[n++]) ^ 2` uses XOR instead
  of `p*p`; prime decomposition broken. (ISSUES.md #53)
- [x] **7.6** `src/cws.cpp:2206` — `W[1].w[0] != W[2].w[0]` compares wrong
  indices; should be `W[1].w[0] != W[1].w[1]`. (ISSUES.md #54)
- [x] **7.7** `src/Rat.cpp:105,355` — `Fgcd`/`LFgcd` divide `a %= b` without
  checking `b == 0`; UB. (ISSUES.md #55)
- [x] **7.8** `src/cws.cpp:16,377` — `lcm`/`Flcm` compute `a*b` before gcd;
  overflow risk. Use `(a/g)*b`. (ISSUES.md #56)
- [x] **7.9** `src/MoriCone.cpp:3055` — `1 << (T.v - 1)` uses `int` literal;
  UB shift when `T.v > 32`. Use `Inci64(1) << ...`. (ISSUES.md #57)
- [x] **7.10** `src/LG.cpp:986` — `puts(i ? "" : " 0");` uses out-of-scope
  loop variable `i`. (ISSUES.md #58)
- [x] **7.11** `src/lgotwist.cpp:1192-1197` — `if (c == EOF)` where `c` is
  `char`; EOF never detected if `char` is unsigned. (ISSUES.md #59)
- [x] **7.12** `src/Subdb.cpp:3257` — `static Along totNF;` in
  `Read_H_ucNF_from_DB` leaks state across different `DB` arguments.
  (ISSUES.md #60)
- [x] **7.13** `src/Polynf.cpp:2447` — `Ind * Ind * Ind` overflows `int`;
  use `Long`. (ISSUES.md #61)
- [x] **7.14** `src/Subdb.cpp:1172,1177,2228,2233,3920,3926` — `fgetc` return
  used as array index without EOF check. (ISSUES.md #62)
- [x] **7.15** `src/MoriCone.cpp:849` — `sumS + > T->I[i];` syntax error in
  dead `#ifdef TRIANG_CHECKSUM` block; remove the block. (ISSUES.md #63)
- [x] **7.16** `src/Polynf.cpp:3865` — `static int x;` incremented but never
  read; dead code. (ISSUES.md #64)
- [x] **7.17** `src/cws.cpp:2206` duplicate check; see 7.6.
- [x] **7.18** Remaining `static` mutable state in `Subdb.cpp` (#15),
  `Polynf.cpp` (#17), `Coord.cpp` (#44), `Vertex.cpp` (#65),
  `Subpoly.cpp` (#66), `MoriCone.cpp` (#67) → anonymous namespace or
  thread-through.

### Phase 8 — RAII: eliminate `malloc`/`free`/`new`/`delete`

~82 `malloc`/`free` pairs and 1 raw `new`/`delete` pair remain. Convert
file-by-file to `std::vector`, `std::unique_ptr`, or stack objects. Each
conversion is one commit.

- [ ] **8.1** `src/MoriCone.cpp` — `Inci64 *A/B/IV/SRG/I` allocations
  (lines 699, 771, 2635, 2909–2911, 2987) → `std::vector<Inci64>` or
  `std::unique_ptr<Inci64[]>`. Already has TODO comments.
- [ ] **8.2** `src/Polynf.cpp` — `FaceInfo *`, `PolyPointList *`, `FibW *`,
  `ek3fli *`, `Long **root` allocations → `std::unique_ptr` / stack.
- [ ] **8.3** `src/E_Poly.cpp` — `Poset_Element *`, `SPoly *`, `Interval *`,
  `BPoly *` calloc/free (lines 1295–1342) → `std::vector`.
- [ ] **8.4** `src/Subadd.cpp` — `PEnt *`, `PPEnt *`, `int *SLp`,
  `unsigned char *NewNF`, `Base_List *`, `NF_List *`, filename `char *`
  (lines 116–119, 351, 437, 457, 756, 758, 1570, 1635, 2123, 2643) →
  `std::vector` / `std::unique_ptr` / `std::string`.
- [ ] **8.5** `src/Subdb.cpp` — `new DataBase` / `delete DB` (lines 3163,
  3253) → `std::unique_ptr<DataBase>` or stack ownership.
- [ ] **8.6** Remaining filename `char *` buffers across `Subdb.cpp` (already
  half-migrated to `std::vector<char>`) → `std::string` + `+=` / `+` /
  `std::format`.

### Phase 9 — `std::array` / `std::vector` / `std::span` for arrays

- [ ] **9.1** Convert struct member C-arrays with compile-time-known sizes to
  `std::array`. Start with the small, leaf structs (`Equation`, `PEnt`,
  `PPEnt`, `BaHo`, `VaHo`, `Step`, `V_Flag`, `M_Rank`, `Subset`, `DxD`,
  `VPerm`). Large structs (`PolyPointList`, `FibW`, `FaceInfo`, `PartList`)
  are deferred until RAII is complete and the impact on stack usage is
  assessed.
- [ ] **9.2** Convert function signatures taking decayed C-array parameters
  (`Long M[][VERT_Nmax]`, `Long *V`, `int *d`) to `std::span<T>` (C++20) or
  reference-to-`std::array` where the size is known. Replace `int *d`/`int *v`
  output parameters with `int& d`/`int& v`.
- [ ] **9.3** Convert VLAs (non-standard C++) to `std::vector`:
  `src/MoriCone.cpp:2004` `Inci64 T[naT]`, `:2006` `int nt[ANfan][ANtri]`,
  `src/SingularInput.cpp:195` `int DegreeVec[dim]`.
- [ ] **9.4** Convert `qsort` + comparator in `src/Vertex.cpp:289` to
  `std::sort` with a lambda.

### Phase 10 — `typedef` → `using` / `struct Name {}`

- [ ] **10.1** Convert all 91 `typedef struct { ... } Name;` to
  `struct Name { ... };` and `typedef T Name[N];` to
  `using Name = std::array<T, N>;` (or `using Name = T[N];` where the
  typedef is used as a function-parameter decay helper). File-by-file.
- [ ] **10.2** Resolve duplicate type definitions: `Matrix` (Mori.h:157 and
  Polynf.cpp:1706 → one shared header), `AmbiPointList` (LG.h:51 forward
  decl vs Nef.h:66 full def → consolidate), `symlist` (LG.cpp:2546 struct vs
  lgotwist.cpp:96 array → rename one), `VPerm`/`VPermList` (Polynf.cpp:1373
  and poly.cpp:85 → shared header), `CWLatticeBasis` (Coord.cpp:9 and
  nef.cpp:24 → shared header), `subl_int` (Subpoly.h:66 vs Subpoly.cpp:8 →
  remove the .cpp duplicate).
- [ ] **10.3** Rename lowercase struct type names to PascalCase:
  `skelet`→`Skelet`, `prili`→`PriLi`, `smon`→`SMon`, `weights`→`Weights`,
  `wei2/3/4`→`Wei2/3/4`, `ek3fli`→`Ek3Fli`, `symlist`(LG)→`SymList`,
  `triang`→`Triang`, `ratmat`/`ratvec`→`RatMat`/`RatVec`.

### Phase 11 — Naming convention migration (LLVM style)

Large, mechanical, file-by-file. Each file is one commit. Keep a
`sed`-friendly rename map; use `clang-tidy`'s `readability-identifier-naming`
check where possible.

- [ ] **11.1** Functions: `Pascal_Snake` → `camelCase`
  (`Read_CWS_PP` → `readCwsPp`, `Find_Equations` → `findEquations`,
  `Init_FInfoList` → `initFInfoList`). Apply consistently across headers and
  sources. The few existing `camelCase` outliers (`auxString2Int`,
  `checkDimension`, `makeN`, `putN`) are already correct.
- [ ] **11.2** Struct/class/typedef names: lowercase → `PascalCase`
  (see 10.3). `UPPER_CASE` struct names (`CWS`, `INCI`, `DYN_PPL`,
  `MORI_Flags`, `NF_List`, `NEF_Flags`) stay as-is if they are acronyms, or
  move to `PascalCase` if they are words (`MORI_Flags` → `MoriFlags`,
  `NF_List` → `NfList`).
- [ ] **11.3** Member fields: `snake_case` / single-letter → `camelCase`
  (`n_nonIP` → `nNonIp`, `sl_nNF` → `slNNf`). Keep single-letter math
  variables (`d`, `n`, `v`) where they mirror the paper notation and are
  documented.
- [ ] **11.4** Parameters: drop the `_`-prefix convention (`_P` → `p`,
  `_V` → `v`, `_E` → `e`, `_CW` → `cw`). No leading underscore (reserved by
  convention for implementation).
- [ ] **11.5** Local variables: UPPER_CASE locals (`CW`, `W`, `BH`, `VH`,
  `FI`, `IN`) → `camelCase` (`cw`, `w`, `bh`, `vh`, `fi`, `in`). Keep
  single-letter loop counters (`i`, `j`, `k`).
- [ ] **11.6** `constexpr` constants: normalize to `UPPER_CASE` for
  compile-time dimension limits (`NFX_Limit` → `NFX_LIMIT`) or `camelCase`
  for config flags. Macros remain `UPPER_CASE`.

### Phase 12 — I/O and string modernization

- [ ] **12.1** Wrap the `fputs(..., stderr); exit(1);` pattern (hundreds of
  occurrences) into a single `[[noreturn]] void palp::die(std::format_string<
  ...>, ...&&)` helper mirroring the existing `Die` in `cws.cpp:75`.
  Eventually convert to `throw std::runtime_error`.
- [ ] **12.2** Replace the 2 remaining `sprintf` calls
  (`Polynf.cpp:4791`, `poly.cpp:375`) with `snprintf` immediately, then with
  `std::format` / `std::format_to` once C++20 `<format>` is available on all
  target compilers.
- [ ] **12.3** Replace `strtok` in `SingularInput.cpp:70,76` with a small
  `std::string_view` tokenizer.
- [ ] **12.4** Replace `atoi` (~23 calls) with `std::stoi` or
  `std::from_chars`.
- [ ] **12.5** Replace `char *` filename members in `NF_List`
  (`Subpoly.h:118`: `iname`, `oname`, `dbname`) with `std::string`.
- [ ] **12.6** Note: the ~1,882 `printf`/`fprintf` calls that produce CLI
  output are kept as-is to guarantee byte-for-byte output compatibility.
  Internal/diagnostic prints can move to `std::format` opportunistically.

### Phase 13 — C-style casts and `goto`

- [ ] **13.1** Replace 212 C-style casts with `static_cast`/`const_cast`/
  `reinterpret_cast`, or eliminate by changing types. The bulk are `(int)`
  casts in `printf` format args that vanish with `<format>`.
- [ ] **13.2** Remove 8 `goto` statements (`Polynf.cpp:1171`,
  `Coord.cpp:324,565`, `Subadd.cpp:1140,1147,1153,1163,1206`) by refactoring
  to structured returns / small lambdas / `do/while(false)`.

### Phase 14 — `exit()` removal from library code

- [ ] **14.1** Introduce `palp::die()` (Phase 12.1) as the single error
  helper. Library functions call `die()` which throws
  `std::runtime_error`; driver `main()` functions catch and `exit(1)`.
- [ ] **14.2** Convert the ~1,000+ `exit(1)` calls in library functions to
  `die()` / `throw`. File-by-file. This is the largest single refactor and
  may be deferred until after RAII/naming if it risks destabilizing the
  test suite.

### Phase 15 — Dead code, comments, and documentation cleanup

- [ ] **15.1** Remove large commented-out code blocks (see ISSUES.md #68):
  `MoriCone.cpp` (~80 lines of diagnostic blocks), `Polynf.cpp:2692-2703`,
  `Subadd.cpp:2399-2402`, `Subdb.cpp:801-804`, `lgotwist.cpp:1313-1322`,
  `cws.cpp:2453-2458`.
- [ ] **15.2** Translate German comments in `Polynf.cpp:5562–6042` to English
  (ISSUES.md #40).
- [ ] **15.3** Remove dead configuration flags `LLong_EEV`/`TEST_EEV`
  (`Vertex.cpp:323-326`) and the `TRIANG_CHECKSUM` block
  (`MoriCone.cpp:844-853`, ISSUES.md #63).
- [ ] **15.4** Address remaining open low-priority ISSUES.md items: #46
  (`binco` cap — architectural), #40 (German comments), #43 (`static EqList`
  in Subdb.cpp), #44 (`static int InputOK` in Coord.cpp).

### Phase 16 — Consolidation and final hardening

- [ ] **16.1** Move shared local typedefs (`CWLatticeBasis`, `Matrix`,
  `VPerm`) into a single `palp/types.h` header.
- [ ] **16.2** Consolidate `lgotwist.cpp` file-scope globals (~20 mutable
  variables) into the `LgoTwistState` struct or an application context.
  Make `mask[]` and `prime[]` `constexpr`.
- [ ] **16.3** Add `const`-correctness: `const T&` for input-only
  parameters, `const` member functions where applicable.
- [ ] **16.4** Enable `clang-tidy` with `readability-identifier-naming`,
  `modernize-*`, `cppcoreguidelines-*` checks and fix remaining warnings.
- [ ] **16.5** Final full Release + ASan + UBSan run: 196/196 for
  `POLY_Dmax` = 4, 5, 6, 11.

---

## Issue / fix tracking

All bugs found during the audit are recorded in `ISSUES.md`. New issues
discovered during Phase 6–16 work should be added there and fixed in the
appropriate phase. The numbering is stable: closed issues keep their numbers.

## Relevant files

- `CMakeLists.txt`: build configuration; line 17 pins `CMAKE_CXX_STANDARD`
  (currently 17, to become 20 in Phase 6.1).
- `include/palp/Global.h`: central header; lines 4–8 legacy C headers;
  lines 32–33 dead `Long`/`LLong` macros; lines 138–232 struct definitions.
- `include/palp/Subpoly.h`: line 3 legacy `<limits.h>`; large struct
  definitions (`DataBase`, `NF_List`, `FInfoList`).
- `ISSUES.md`: bug and issue tracker.
- `AUDIT.md`: raw audit findings (generated by the Phase 6 codebase audit).