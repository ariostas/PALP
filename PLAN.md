# PALP Modernization Plan: C → C++17

## Goal

Migrate PALP from C to C++17 while preserving byte-for-byte identical CLI
output and behavior. Bugs are documented in `ISSUES.md` and fixed in a
separate phase after migration.

## Principles

- **Behavior preservation first**: every step must keep the full test suite
  passing.
- **One file per step**: each step converts a single `.c` file to `.cpp` (or
  a single cross-cutting concern). Verify with `cmake --build build && ctest`
  after each step.
- **Commit after every step**: once verification passes, the changes for that
  step must be committed before moving on. Use a descriptive commit message
  (e.g., "Migrate Rat.c to Rat.cpp"). This keeps history reviewable and makes
  it easy to bisect if a later step breaks a test.
- **POLY_Dmax stays compile-time**: keep `constexpr` equivalent of current
  macros; do not switch to dynamic allocation in the migration phase.
- **No new features**: don't refactor algorithms; only modernize the language
  constructs.
- **Exact output compatibility**: every test in `tests/` must produce identical
  output before and after each step.

## Test harness (run before starting)

1. Build baseline:
   ```bash
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build
   ```
2. Run tests:
   ```bash
   cd build && ctest
   ```
3. Build and run sanitizer builds:
   ```bash
   cmake -S . -B build/asan -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_C_FLAGS="-O1 -g -fno-omit-frame-pointer -fsanitize=address" \
     -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address"
   cmake --build build/asan
   cd build/asan && ctest

   cmake -S . -B build/ubsan -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_C_FLAGS="-O1 -g -fno-omit-frame-pointer -fsanitize=undefined -fno-sanitize-recover=undefined" \
     -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=undefined"
   cmake --build build/ubsan
   cd build/ubsan && ctest
   ```
4. Record golden outputs from all test scripts for regression checking (the
   test scripts already embed expected outputs, so this is implicit).

---

## Phases

### Phase 0 — Infrastructure

- [x] #### Step 0.1 — Update CMakeLists.txt for C++ support

- Change `LANGUAGES C` to `LANGUAGES C CXX`.
- Add `set(CMAKE_CXX_STANDARD 17)` and `set(CMAKE_CXX_STANDARD_REQUIRED ON)`.
- Add a helper function or macro that can compile a PALP source as either `.c`
  or `.cpp`, so both can coexist during migration.
- Add `target_compile_options(... PRIVATE -Wall -Wextra)` for C++ targets.
- Keep the existing C targets working unchanged.
- **Verify**: `cmake --build build && cd build && ctest` — all pass.

- [x] #### Step 0.2 — Reorganize headers into `include/palp/`

- Create `include/palp/` directory.
- Copy (not move yet) `Global.h`, `Rat.h`, `LG.h`, `Nef.h`, `Mori.h`,
  `Subpoly.h` into `include/palp/` with `#pragma once` added.
- Add `include/palp/` to `target_include_directories` for all targets.
- Update all `#include "Global.h"` etc. to `#include "palp/Global.h"` in all
  `.c` files.
- Remove the old top-level header files (or leave as forwarding shims during
  migration).
- **Verify**: build + test.

- [x] #### Step 0.3 — Create C++ compatibility shim

- Create `include/palp/palp_types.h` defining the integer type aliases in a
  C++-friendly way:
  ```cpp
  #include <cstdint>
  using Long = long;
  using LLong = long long;
  // ... etc
  ```
- This lets `.cpp` files include PALP headers without C-isms breaking.
- **Verify**: a trivial `.cpp` file that includes `palp/Global.h` compiles.

---

### Phase 1 — Foundation libraries (leaf modules)

These modules have no dependencies on other PALP modules (only on standard C
library and each other in a simple chain). Convert in dependency order.

- [x] #### Step 1.1 — `Rat.c` → `Rat.cpp`

- [x] Rename `Rat.c` → `Rat.cpp` (mixed C/C++ linkage was temporarily ensured
  via `extern "C"`; removed after full C++ migration).
- [x] Remove `register` keyword (all instances).
- [x] Replace `#ifdef TEST` test blocks: none present in `Rat.c`/`Rat.cpp`.
- [x] Convert `Rat`/`LRat` structs: add C++ constructors, keep as POD-like
  structs.
- **Verify**: build + test (especially `tests/2.*` which exercise `poly.x`
  which uses `Rat`).

- [x] #### Step 1.2 — `Vertex.c` → `Vertex.cpp`

- [x] Rename `Vertex.c` → `Vertex.cpp`.
- [x] Replace `malloc`/`free` for temp arrays (`CEq`, `CEq_I`, `F_I`) with
  `std::vector` or `std::make_unique`.
- [x] Keep `INCI` macros and typedefs as-is (performance-critical bit operations).
- [x] Keep `exit(0)` calls for now (bug fix is Phase 5).
- **Verify**: build + test (all `tests/2.*`, `tests/3.*`).

- [ ] #### Step 1.3 — `Coord.c` → `Coord.cpp`

- [x] Rename `Coord.c` → `Coord.cpp`.
- [x] Replace some `char c[999]` local buffers with `std::array<char, 999>` or
  `std::vector`.
- [ ] Replace remaining `fscanf`-based parsing: add return-value checks but keep
  identical output behavior (check return, on failure produce same error message).
- [x] `ReadCwsPp`, `Read_PP`, `Read_CWS` now check `fscanf` return values and exit
  with a descriptive error message on failure.
- [ ] Keep `static int InputOK` as-is for now (documented in ISSUES.md #43).
- **Verify**: build + test (all `tests/2.*`, `tests/4.*`).

- [x] #### Step 1.4 — `Polynf.c` → `Polynf.cpp` (largest: ~3223 lines)

This is the highest-risk step. Take extra care.

- [x] Rename `Polynf.c` → `Polynf.cpp`.
- [x] Replace `#if (VERT_Nmax < 129)` stack-vs-heap conditionals with RAII
  allocations (completed in earlier work; none remain).
- [x] `volatile` was not present in `Polynf.cpp`.
- [x] `register` keyword was not present in `Polynf.cpp`.
- [x] `static int` counters (`ConifoldSing` statistics, `Print_Fiber_PolyData`'s
  `f`, `GL_Lattice_Basis`'s `x`, `Fano5d`'s `FPc`) are left as `static` and
  documented (ISSUES.md #16): making them explicit parameters would require
  changing public/internal call signatures and output behavior.
- [x] Replace `#ifdef TEST`/`#ifdef TEST_OUT`/`#ifdef TEST_ImpPhase`/etc. debug
  blocks with `constexpr bool` flags at the top of `Polynf.cpp`.
- [x] Remove `puts("PM")` debug print (ISSUES.md #33).
- **Verify**: build + test for ALL dimensions (4, 5, 6, 11) — run full
  `ctest` with `DIM` variations.

- [x] #### Step 1.5 — `LG.c` → `LG.cpp` (~1113 lines)

- [x] Rename `LG.c` → `LG.cpp`.
- [x] Remove `register` keyword (all instances).
- [x] Replace `char c[999]` in `Read_WZ_PP` with `std::array<char, 999>`.
- [x] Replace `AllocPoCoLi`/`Free_PoCoLi` manual pointer arithmetic with
  `std::vector<int> e` and `std::vector<Pint> c` inside `PoCoLi`. Rewrite
  `PoincarePoly` local `B` to use `AllocPoCoLi` instead of stack arrays.
- [x] Replace `#ifdef TEST`/`#ifdef TEST_PP`/`#ifdef TEST_PD` debug blocks with
  `constexpr bool` flags at the top of `LG.cpp`; remove mid-file `#define TEST`/
  `#undef TEST` toggles. Keep `#if (WZinput)` macro as-is (always 1).
- [x] `static int MaxPoNum` (inside `TEST_WeightMakePoints`) and `static int M`
  (inside `Add_Mono_2_Poly` debug block) remain in compile-time-disabled debug
  code; documented as harmless (ISSUES.md #15).

- [ ] #### Step 1.6 — `Nefpart.c` → `Nefpart.cpp` (~837 lines)

- [x] Rename `Nefpart.c` → `Nefpart.cpp`.
- [ ] Remove `scanf("%c", &c)` interactive debug pause (ISSUES.md #38) or gate
  behind `#ifndef NDEBUG`.
- [ ] Replace bubble sort (`Bubble_PTL`) with `std::sort` — verify identical
  ordering (the sort key is the partition string; `std::sort` with a custom
  comparator must produce the same order).
- [ ] Replace `calloc`/`malloc` with `std::vector`.
- **Verify**: build + run all `tests/6.*` scripts.

---

### Phase 2 — Driver programs (depend on Phase 1 libraries)

- [ ] #### Step 2.1 — `poly.c` → `poly.cpp`

- [x] Rename `poly.c` → `poly.cpp`.
- [ ] Replace `malloc` for large compile-time-sized structs (`CWS`, `EqList`,
  `PolyPointList`, `PairMat`, `FaceInfo`) with `std::make_unique` or stack
  allocation (they are fixed-size).
- [ ] Replace `FILE *inFILE, *outFILE` global definitions: keep as global for now
  (the libraries extern-reference them), but add a comment that this is
  temporary (ISSUES.md #40).
- [ ] Replace `VPermList *VP = (VPermList*) malloc(...)` with `std::make_unique`.
- Keep `printf`/`fprintf` as-is for exact output format preservation.
- **Verify**: build + run ALL `tests/2.*` and `tests/3.*` scripts.

- [ ] #### Step 2.2 — `cws.c` → `cws.cpp` (~1935 lines)

- [x] Rename `cws.c` → `cws.cpp`.
- [ ] Replace `WDIM=800000` and `TWDIM=16384` magic constants with `constexpr int`.
- [ ] Replace VLA `int IN[AMBI_Dmax*(AMBI_Dmax+1)]` with `std::array` or
  `std::vector`.
- [ ] Remove dead `FileRW()` function (ISSUES.md #6).
- [ ] Replace `sprintf(command, ...)` with `snprintf` (behavior-preserving, safer).
- [ ] Replace `atoi("3")` / `atoi("4")` pointless conversions with direct
  integer literals.
- **Verify**: build + run all `tests/4.*` scripts.

- [ ] #### Step 2.3 — `class.c` → `class.cpp`

- [x] Rename `class.c` → `class.cpp`.
- [ ] Replace `char Blank=0` pointer-to-empty-string hack with `std::string` and
  `.c_str()` where needed. Careful: many functions take `char*` — keep
  `.c_str()` pointers valid for the duration of use.
- [ ] Replace `scanf("%s", &hc)` interactive help with a bounded read
  (`std::cin >> hc` or `scanf` with length check).
- [ ] Keep `FILE *inFILE, *outFILE` global definition here (class.x defines them).
- **Verify**: build + run any class-related tests (check if class tests exist
  in `tests/`; if not, create a basic smoke test).

- [ ] #### Step 2.4 — `nef.c` → `nef.cpp`

- [x] Rename `nef.c` → `nef.cpp`.
- [ ] Replace local typedefs (`AmbiLatticeBasis`, `CWLatticeBasis`, `Pstat`) with
  proper struct definitions; check if they duplicate types from `Coord.c`.
- [ ] Replace VLA `Long PM[EQUA_Nmax][VERT_Nmax]` with `std::vector`.
- [ ] Fix `int long nl` → `long int nl` (or just `long nl`) (ISSUES.md #22).
- **Verify**: build + run all `tests/6.*` scripts.

- [ ] #### Step 2.5 — `mori.c` → `mori.cpp`

- [x] Rename `mori.c` → `mori.cpp`.
- [ ] Replace all `malloc` calls with `std::make_unique` or stack allocation.
  Fix the memory leaks (ISSUES.md #8) by ensuring RAII cleanup.
- [ ] Keep `FILE *inFILE, *outFILE` global definition here (mori.x defines them).
- **Verify**: build + run all `tests/7.*` scripts.

---

### Phase 3 — Secondary libraries (depend on Phase 1 + 2)

- [ ] #### Step 3.1 — `E_Poly.c` → `E_Poly.cpp` (~1544 lines)

- [x] Rename `E_Poly.c` → `E_Poly.cpp`.
- [ ] Remove local `#define min`/`#define max` (ISSUES.md #28); use `std::min`/
  `std::max` from `<algorithm>`.
- [ ] Fix `int h[POLY_Dmax][POLY_Dmax] = {{0},{0}}` → `int h[POLY_Dmax][POLY_Dmax]
  = {}` (ISSUES.md #4).
- [ ] Replace `realloc` in `DYNadd_for_completion` with `std::vector` growth
  (ISSUES.md #3).
- [ ] Replace `Die()` function with a `[[noreturn]]` C++ function.
- **Verify**: build + run all `tests/6.*` scripts (nef uses E_Poly).

- [ ] #### Step 3.2 — `Nefpart.c` → `Nefpart.cpp` (~837 lines)

- [x] Rename `Nefpart.c` → `Nefpart.cpp`.
- [ ] Remove `scanf("%c", &c)` interactive debug pause (ISSUES.md #38) or gate
  behind `#ifndef NDEBUG`.
- [ ] Replace bubble sort (`Bubble_PTL`) with `std::sort` — verify identical
  ordering (the sort key is the partition string; `std::sort` with a custom
  comparator must produce the same order).
- [ ] Replace `calloc`/`malloc` with `std::vector`.
- **Verify**: build + run all `tests/6.*` scripts.

- [ ] #### Step 3.3 — `MoriCone.c` → `MoriCone.cpp` (~1792 lines)

- [x] Rename `MoriCone.c` → `MoriCone.cpp`.
- [ ] Fix `assert(++m < binco)` (ISSUES.md #1): move `++m` out of assert:
  ```cpp
  ++m; assert(m < binco);
  ```
- [ ] Replace `Inci64` macros (`makeN`, `putN`, `getN`, etc.) with `constexpr`
  inline functions.
- [ ] Remove dead `#ifdef OLD_code` and `#ifdef FIRST_TRY__TOO_COMPLICATED...`
  blocks (ISSUES.md #35).
- [ ] Remove stray `printf` debug (ISSUES.md #34).
- [ ] Replace `exit(0)` / `exit(1)` with `palp::die()`.
- **Verify**: build + run all `tests/7.*` scripts.

- [ ] #### Step 3.4 — `SingularInput.c` → `SingularInput.cpp` (~591 lines)

- [x] Rename `SingularInput.c` → `SingularInput.cpp`.
- [ ] Replace `char filename[20]` with `std::string`.
- [ ] Replace `char string[maxline]` VLA with `std::vector<char>`.
- [ ] Add `malloc` NULL check for `he` (ISSUES.md #10).
- [ ] Keep `system()` call but add note (ISSUES.md #11); in migration phase, just
  ensure the command string is properly constructed.
- [ ] Replace `mkstemp` usage: keep as-is (it's the right function) or use C++17
  `std::filesystem::temp_directory_path()` for the directory.
- **Verify**: build + run all `tests/7.*` scripts.

- [ ] #### Step 3.5 — `Subpoly.c` → `Subpoly.cpp` (~1614 lines)

- [x] Rename `Subpoly.c` → `Subpoly.cpp`.
- [ ] Replace `subl_int` typedef with explicit `int64_t`.
- [ ] Replace `exit(0)` with `palp::die()`.
- [ ] Replace `drop_point[POLY_Dmax]` "silence compiler" zero-init with `= {}`
  (ISSUES.md, code smell).
- [ ] Replace `malloc`/`calloc` with `std::vector` or `std::make_unique`.
- **Verify**: build + run class-related tests.

- [ ] #### Step 3.6 — `Subadd.c` → `Subadd.cpp` (~1408 lines)

- [x] Rename `Subadd.c` → `Subadd.cpp`.
- [ ] Replace `fscanf(F, "%c%c%c%c", &A, &B, &C, &D)` with `fread` for binary I/O
  (more correct and faster).
- [ ] Replace `unsigned char auxUC[POLY_Dmax*VERT_Nmax]` VLA with `std::vector`.
- [ ] Replace `NF_List *AuxNFLptr = NULL` global (ISSUES.md #17) with explicit
  parameter passing where feasible; otherwise document.
- [ ] Replace `fgetUI`/`fputUI` with `fread`/`fwrite`-based versions.
- **Verify**: build + run class-related tests.

- [ ] #### Step 3.7 — `Subdb.c` → `Subdb.cpp` (~1915 lines)

- [x] Rename `Subdb.c` → `Subdb.cpp`.
- [ ] Replace `goto END_SL` / `goto END_VN` with structured flow control (lambdas
  or early-return helper functions).
- [ ] Replace `static unsigned char uc[NUC_Nmax]` and `static int ms3`
  (ISSUES.md #14) with a context struct passed as parameter.
- [ ] Replace `sprintf(com, ...)` with `snprintf` (ISSUES.md #26).
- [ ] Replace `static EqList E` (ISSUES.md #42) with local variable.
- **Verify**: build + run class-related tests.

---

### Phase 4 — Standalone & cleanup

- [x] #### Step 4.1 — `lgotwist.c` → `lgotwist.cpp` (standalone, not in CMake)

- [x] Rename `lgotwist.c` → `lgotwist.cpp` (no CMake target).
- [ ] Remove duplicated rational arithmetic (ISSUES.md #36); use `Rat.cpp`'s
  functions or keep standalone but at least use the same types.
- [ ] Replace `#define abs/min/max/mod` (ISSUES.md #31) with standard library
  functions.
- [ ] Replace global variables (ISSUES.md, lgotwist.c global state) with a
  `LgoTwistContext` struct.
- [ ] Remove `register` keyword.
- [ ] Replace `#ifdef __MSDOS__` platform checks with C++17 equivalents.
- **Verify**: build standalone (not in test suite; manual verification).

- [x] #### Step 4.2 — Remove `.c` files and legacy build files

- [x] Once all `.cpp` files pass all tests, remove the original `.c` files from
  `CMakeLists.txt`.
- [x] Remove the deprecated `GNUmakefile` and root `Makefile` (CMake is the
  primary and only supported build system).
- [x] Move all C++ sources under `src/` and `Rat.h` under `include/palp/`.
- [x] Remove `extern "C"` shims from all headers and source files now that the
  project is C++-only.
- [x] Remove `extern "C"` shims from all headers and source files now that the
  project is C++-only.
- [x] Remove `extern "C"` shims from all headers and source files now that the
  project is C++-only.
- [x] Add proper include guards (now `#pragma once`) to all headers
  (`Global.h`, `LG.h`, `Nef.h`, `Mori.h`, `Subpoly.h`, `Rat.h`,
  `palp_types.h`).
- [x] Consolidate `min`/`max` definitions: introduce `palp::min`/`palp::max`
  templates in `Global.h`, remove macros from `LG.h`, `Subpoly.h`, and
  `E_Poly.cpp`. Keep standalone `lgotwist.cpp` macros for now.
- [ ] Clean up any remaining forwarding shims from Step 0.2.
- **Verify**: full clean build + test.

- [ ] #### Step 4.3 — Final cleanup

- [x] Remove `palp_types.h` shim if all types are now properly C++.
  *Done earlier; no separate `palp_types.h` exists.*
- [x] Consolidate `min`/`max` definitions: delete from `LG.h`, `Subpoly.h`,
  `E_Poly.cpp` (all instances of `#define min/max`).
  *`palp::min`/`palp::max` templates are in `Global.h`; no `#define min/max`
  remain in those files. Standalone `lgotwist.cpp` keeps its own for now.*
- [ ] Remove the `-DNDEBUG` workaround in CMake (line 15 of `CMakeLists.txt`)
  if asserts have been replaced with proper error handling in Phase 5.
  If not yet, keep it.
- [ ] Convert remaining header-level `#define` constants that are safe to
  `constexpr`/`using` without breaking `#if` array-size logic (e.g. type aliases
  and scalar limits not used in preprocessor conditionals).
- [ ] Run full test suite across all dimensions:
  ```bash
  for DIM in 4 5 6 11; do
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPOLY_Dmax=$DIM
    cmake --build build
    cd build && DIM=$DIM ctest && cd ..
  done
  ```
- [ ] Run ASAN and UBSAN builds one final time.
- **Verify**: everything green.

### Step 4.4 — `#define` → `constexpr`/`using` sweep (in progress)

A cross-cutting pass to replace file-local `#define` constants with C++17
`constexpr`/`using` where it does not change preprocessor-controlled array
sizes or conditional-compilation behavior. Symbols that are still required by
`#if`/`#ifdef` guards are left as macros until their guarded blocks are
converted to `if constexpr` or removed.

Files already converted:
- `src/cws.cpp`: `Only_IP_CWS`, `TRANS_INFO_FOR_IP_WEIGHTS`, `NFmax`,
  `SIMPLEX_POINT_Nmax`, `OSL`, `WDIM`, `lcm`, `TWDIM`, `mod`, `ALLOWHALF`,
  `CHAT`, plus dead-code removal.
- `src/SingularInput.cpp`: `DijkEQ`, `T_DIV`, `DIVclassBase`,
  `TEST_PRINT_SINGULAR_IO`, `NORM_SIMP_NUM`, plus `const` correctness for
  `DivClassBasis` string arguments.
- `src/Vertex.cpp`: `SHOW_NEW_CEq`, `LLong_EEV`, `TEST_EEV`,
  `VERT_WITH_MAX_DISTANCE`, `LONG_EQ_FIRST`, `TEST_GLZ_EQ`, `MAX_BAD_EQ`
  (converted to `constexpr bool` together with its `#if` block).
- `src/Subdb.cpp`: `SUBTRACT_H_FROM_SL`, `Hod_Dif_max`, `Hod_Min_max`, all
  local `TEST`/`TEST_OUT` toggles.
- `src/Polynf.cpp`: `SORT_CWS`, `FIB_PERM`, `SSR_PRINT`,
  `ALL_FANOS_BUT_INEFFICIENT`, `FANO_CONIFOLD`, `SL_Long`, `KPF`,
  `RelativeSimplexVolume`, `No_OLD_FACE_LIST`, `SQnum_Max`, `TESTfano`,
  `FanoProjNPmax`, `FPcirNmax`, `PrintFanoProjCand`, `INCIbits`, `Fputs`
  to an inline function, `SMOOTH`, `NON_REF`, `BARY_PRINT`,
  `ZEROSUM_PRINT`, `KP_PRINT`, `SHOW_NFX_LIMIT`. Active
  preprocessor-dependent flags whose values are used for array sizes
  (`NFX_Limit`, `X_Limit`, `VPM_Limit`) and `KP_VALUE`/`KP_EXIT` remain
  as macros.
- `src/LG.cpp`: `SHOW_b01_TWIST`, `Tout`, `DET_WARN_ONLY`,
  `ABBREV_POLY_PRINT`, `NO_COORD_IMPROVEMENT` (kept as macro to keep dead
  blocks disabled), `TEST_LG`, `TEST_PP`, `TEST_PD`, `StandardOutput`.
- `src/E_Poly.cpp` / `src/nef.cpp`: `WRITE_CWS` replaced with local
  `constexpr bool write_cws`; macro removed from `include/palp/Nef.h`.
- `src/poly.cpp` / `src/nef.cpp`: `OSL`.
- `src/MoriCone.cpp`: `ANfan`, `ANtri`, `Inci64_AND`, `Inci64_EQ`,
  `Inci64_0`, `Inci64_1`, `Inci64_EQ_0`, `Inci64_OR`, `Inci64_PN`,
  `Inci64_D2`, `Inci64_M2`, plus fix of `assert(++m<binco)` side-effect bug.
- `src/Subadd.cpp`: compression constants `UCM`, `USM`, `Nint_XLong`, `NX`,
  `UNIT_OFF`, `LL_BASE`, and `MirTest` to an inline function. Dead
  `MOVE_SAVE_FILE` block removed.
- `include/palp/Global.h`: `GL_Long`, `MAXLD`, `INT_Nbits`,
  `LONG_LONG_Nbits`, `I_NUI`, all single-integer `INCI_*` bit helpers,
  `EQUA_Nmax` fallback, `CEQ_Nmax`, and the large-`INCI` branch of
  `INCI_M2`.
- `include/palp/Nef.h`: `Nef_Max`, `NP_Max`, `MAXSTRING`, `Pos_Max`,
  `FIB_POINT_Nmax` converted to `constexpr`. `W_Nmax` remains a macro
  because both `Nef.h` and `LG.h` define it for array sizes.
- `include/palp/LG.h`: `Pint`.
- `include/palp/Subpoly.h`: `NUC_Nmax`, `MAX_REC_DEPTH`, `Along`, `UPint`,
  `FORCE_SAVE_TIME`, `GOOD_SAVE_TIME`, `WRITE_DIM`, `MIN_NEW`,
  `MIN_W_SAVE_TIME`, `WATCHREF`, `SAVE_INC`, `CperR_MAX`, `BASE_MAX`,
  `BLOCK_LENGTH`, `subl_int`, `NB_MAX`, `SL_Nmax`, `File_Ext_NCmax`,
  `SAVE_FILE_EXT`; preprocessor consistency check replaced with
  `static_assert`.

Still to convert in `.cpp` files:
- `src/Coord.cpp`: `NO_COORD_IMPROVEMENT` is defined, so all
  `#ifndef NO_COORD_IMPROVEMENT` blocks are dead; removing them requires care
  because `Wperm_to_GLZ`/`CWS_to_PermCWS` are declared in prototypes.
- `src/LG.cpp`: `COEFF_Nmax` depends on local variables and is used in array
  sizes; keep as macro for now.
- `src/MoriCone.cpp`: the local geometric helper macros `BZangle`,
  `SameRayBZ`, `BZR`, `BZRx`, `BZRE` capture many local variables from the
  enclosing functions and are used in many contexts; intentionally left as
  expression macros to avoid large risky refactor.
- `src/Subadd.cpp`: the `#ifdef` feature flags `TEST_UCnf`,
  `ADD_LIST_LENGTH`, `INCREMENTAL_TIME`, `INCREMENTAL_WRITE`,
  `ACCEL_PEntComp`, `USE_UNIT_ENCODE` are all currently enabled and deeply
  interleaved with function bodies. Converting them to `if constexpr`
  requires careful untangling and is deferred.
- `src/lgotwist.cpp`: standalone; lower priority.

Still to convert in headers:
- `include/palp/Global.h`: `POLY_Dmax`, `POINT_Nmax`, `VERT_Nmax`,
  `FACE_Nmax`, `SYM_Nmax`, `EQUA_Nmax`, `AMBI_Dmax`, `FIB_Nmax`,
  `CD2F_Nmax` remain compile-time sizing macros; `MULTIPLYING` also
  remains because it controls a global `#if` branch. Converting these to
  `constexpr` would require dynamic allocation of all dependent arrays and
  is out of scope for this sweep.
- `include/palp/LG.h`: `WZinput` (used in `#if`), `W_Nmax`.
- `include/palp/Nef.h`: `W_Nmax`.
- `include/palp/Subpoly.h`: `USE_TMP_DIR` (used in `#if`). `FTELL`/`FSEEK`
  are now inline wrapper functions.
- `include/palp/Rat.h`: `ARG_FUN` converted to a `using` alias.
- `include/palp/Global.h`: the large-`VERT_Nmax` `INCI_*` branch is not
  currently compiled (`VERT_Nmax <= 64`) and can be converted later.

### Step 4.5 — Use `std::array` for fixed-size local buffers (pending)

Some migration changes introduced `std::vector` where the size is actually
known at compile time. Once the `#define` constants they depend on become
`constexpr` (or where they already are), these local buffers should be
replaced with `std::array` to avoid heap allocation and express the fixed-size
intent:

- `src/Coord.cpp`: `std::vector<int> IN(AMBI_Dmax*(AMBI_Dmax+1))` in
  `ReadCwsPp`, `Read_PP`, and `Read_CWS`.
- `src/Vertex.cpp`: `std::vector<INCI> CEq_I(CEQ_Nmax)` and
  `std::vector<INCI> F_I(EQUA_Nmax)` in `Finish_IP_Check`,
  `Find_Equations`, and `IP_Check`.

`PoCoLi` in `include/palp/LG.h` should remain `std::vector` because its size
is determined at runtime by the Poincaré polynomial being computed.

---

### Phase 5 — Bug fixes (separate from migration)

Each bug fix is a separate step with a regression test added to the test
suite. See `ISSUES.md` for the full list. Apply in order of severity.

#### Step 5.1 — Fix `assert(++m < binco)` side-effect bug

- [ ] File: `MoriCone.cpp` (was `MoriCone.c:408`).
- [ ] Move increment out of assert.
- [ ] Add regression test: a polytope that triggers the binco path.

#### Step 5.2 — Fix `exit(0)` → `exit(1)` for error paths

- [ ] All files: replace `exit(0)` on error conditions with `exit(1)`.
- [ ] This is a large mechanical change; do it in one step with careful testing.
- [ ] Tests check stdout, not exit codes, so output should be identical.

#### Step 5.3 — Fix memory leaks

- [ ] `mori.cpp`: free all allocations (ISSUES.md #8).
- [ ] `cws.cpp` `RgcWeights`: free `X` (ISSUES.md #9).
- [ ] `Polynf.cpp` `Eval_Poly_NF`/`Fano5d`: free on error paths (ISSUES.md #12, #13).
- [ ] Run ASAN build to verify zero leaks.

#### Step 5.4 — Fix `realloc` losing old pointer

- [ ] `E_Poly.cpp` `DYNadd_for_completion` (ISSUES.md #3).
- [ ] Use `std::vector` which handles this correctly.

#### Step 5.5 — Fix partial array initialization

- [ ] `E_Poly.cpp` (ISSUES.md #4): `= {{0},{0}}` → `= {}`.

#### Step 5.6 — Add missing `fscanf` return-value checks

- [ ] All files (ISSUES.md #27): wrap `fscanf` calls, check return, produce
  identical error messages on failure.

#### Step 5.7 — Fix `system()` command injection

- [ ] `SingularInput.cpp` (ISSUES.md #11): sanitize `getenv("TMPDIR")` or use
  `std::filesystem::temp_directory_path()`.

#### Step 5.8 — Fix buffer overflow risks

- [ ] `LG.cpp` `char c[999]` (ISSUES.md #23): use `std::string` with bounds-checked
  reading.
- [ ] `lgotwist.cpp` `s->p[s->N++]` (ISSUES.md #24): add bounds check.
- [ ] `cws.cpp` `char command[100]` (ISSUES.md #25): use `snprintf`.
- [ ] `Subdb.cpp` `char com[35]` (ISSUES.md #26): use `snprintf`.

#### Step 5.9 — Fix integer overflow risks

- [ ] `Polynf.cpp` volume computation (ISSUES.md #18): add overflow checks or use
  wider types.
- [ ] `lgotwist.cpp` `lcm` macro (ISSUES.md #19): divide before multiply.

#### Step 5.10 — Fix `assert` used as control flow

- [ ] Replace critical asserts (those guarding data correctness, not just internal
  consistency) with explicit `if` checks + error handling.
- [ ] `Polynf.cpp:894` `assert(g > 0)` (ISSUES.md #5).
- [ ] `Vertex.cpp:154` Euler characteristic check (ISSUES.md #47).

#### Step 5.11 — Fix global `inFILE`/`outFILE` state

- [ ] Replace with a `PalpContext` struct holding `inFILE`/`outFILE` and pass
  explicitly to library functions (ISSUES.md #40, #41, #43).
- [ ] This is a large refactor; do it last and carefully.

#### Step 5.12 — Remove dead code

- [ ] Remove `FileRW()` (ISSUES.md #6).
- [ ] Remove `#ifdef OLD_code` blocks (ISSUES.md #35).
- [ ] Remove debug `puts("PM")`, `printf` (ISSUES.md #33, #34).
- [ ] Remove `#define TEST`/`#undef TEST` toggling (ISSUES.md #37).

---

## Tracking completed work

Mark each step with `- [x]` as it is completed and committed. Keep this list
up to date so the current state of the migration is always visible at a
glance. The next unmarked item is the current step.

## Verification Protocol (after every step)

1. **Build**: `cmake --build build` must succeed with `-Wall -Wextra` and no
   warnings.
2. **Tests**: `cd build && ctest` must pass 100%.
3. **ASAN**: `cmake --build build/asan && ctest --test-dir build/asan` must
   pass with no address sanitizer errors.
4. **UBSAN**: `cmake --build build/ubsan && ctest --test-dir build/ubsan` must
   pass with no undefined behavior sanitizer errors.
5. **Spot-check**: for steps touching core libraries, diff output of a
   representative input across dims 4, 5, 6, 11 against pre-migration golden
   output.
6. **Git**: after each successful step, commit with a descriptive message
   (e.g., "Migrate Rat.c → Rat.cpp"). The commit is part of the step; no
   subsequent step may begin until the previous one is committed and the
   repository is in a clean state.

---

## File dependency graph (conversion order)

```
Rat.c        (leaf)           → Step 1.1
Vertex.c     (← Rat, Global)  → Step 1.2
Coord.c      (← Rat, Global)  → Step 1.3
Polynf.c     (← Rat, Global)  → Step 1.4
LG.c         (← Rat, Global)  → Step 1.5
poly.c       (← LG, Global)   → Step 2.1
cws.c        (← LG, Global)  → Step 2.2
class.c      (← Subpoly)      → Step 2.3
nef.c        (← Nef, LG)      → Step 2.4
mori.c       (← Mori, LG)    → Step 2.5
E_Poly.c     (← Nef, Rat)     → Step 3.1
Nefpart.c    (← Nef)          → Step 3.2
MoriCone.c   (← Rat, Mori)    → Step 3.3
SingularInput.c (← Mori)     → Step 3.4
Subpoly.c    (← Rat, Subpoly.h) → Step 3.5
Subadd.c     (← Subpoly)      → Step 3.6
Subdb.c      (← Subpoly, Rat) → Step 3.7
lgotwist.c   (standalone)     → Step 4.1
```

---

## Risk assessment

| Step | Risk | Reason |
|------|------|--------|
| 1.4 (Polynf) | High | Largest file (3223 lines), complex math, many conditionals |
| 2.2 (cws) | High | 1935 lines, many code paths, hardcoded tables |
| 3.7 (Subdb) | Medium | File I/O, goto flow control, database logic |
| 1.5 (LG) | Medium | Polynomial arithmetic, manual memory management |
| 3.3 (MoriCone) | Medium | Complex geometry, assert side-effect bug |
| All others | Low-Medium | Standard C-to-C++ conversion |