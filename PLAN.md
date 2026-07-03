#PALP Modernization Plan : C → C++ 17

## Goal

Migrate PALP from C to C++17 while preserving byte-for-byte identical CLI
output and behavior. Bugs are documented in `ISSUES.md` and fixed in a
separate phase after migration.

## Principles

- **Behavior preservation first**: every step must keep the full test suite passing.
- **One file per step**: convert a single source or cross-cutting concern at a time.
- **Run `clang-format -i -style=LLVM` before every commit**.
- **Commit after every verified step** with a descriptive message.
- **POLY_Dmax stays compile-time**;
no dynamic allocation in the migration phase.-
    **No new features ** : modernize language constructs only,
    not algorithms.- **Exact output compatibility ** : every `tests/*.sh` output must match the pre-migration baseline.

## Test harness

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build

# Sanitizer builds (Phase 5)
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

- [x] Step 0.1 — `CMakeLists.txt`: CXX-only, C++17, `-Wall -Wextra`, `palp_sources()` helper.
- [x] Step 0.2 — Headers moved to `include/palp/` with `#pragma once`; include paths updated.
- [x] Step 0.3 — `palp_types.h` shim added and integrated via `Global.h`.

### Phase 1 — Foundation libraries

- [x] Step 1.1 — `Rat.c` → `Rat.cpp`: constructors, removed `register`.
- [x] Step 1.2 — `Vertex.c` → `Vertex.cpp`: RAII temp arrays, kept `INCI` macros.
- [x] Step 1.3 — `Coord.c` → `Coord.cpp`: `std::array` buffers, `fscanf` checks.
- [x] Step 1.4 — `Polynf.c` → `Polynf.cpp`: largest file; RAII, `constexpr` debug flags, removed dead code.
- [x] Step 1.5 — `LG.c` → `LG.cpp`: `std::array`/`std::vector`, `PoCoLi` rewritten, debug flags modernized; later split `MakeMobius` packed buffer into separate `data_storage` and `mt_storage` to fix heap overflow / misaligned store.

### Phase 2 — Driver programs

- [x] Step 2.1 — `poly.c` → `poly.cpp`: `std::unique_ptr` for core structs.
- [x] Step 2.2 — `cws.c` → `cws.cpp`: `constexpr`, `std::array`, removed dead `FileRW()`, `snprintf`.
- [x] Step 2.3 — `class.c` → `class.cpp`: bounded `scanf`, `std::unique_ptr<_P>`.
- [x] Step 2.4 — `nef.c` → `nef.cpp`: `using`/struct aliases, `std::vector`/`std::unique_ptr`, `Die()` modernization.
- [x] Step 2.5 — `mori.c` → `mori.cpp`: `std::unique_ptr` for core structs.

### Phase 3 — Secondary libraries (mostly done)

- [x] Step 3.2 — `Nefpart.c` → `Nefpart.cpp`: `std::sort`, RAII, removed debug pause.
- [x] Step 3.4 — `SingularInput.c` → `SingularInput.cpp`: `std::string`/`std::vector`.
- [x] Step 3.5 — `Subdb.c` → `Subdb.cpp`: filename buffers modernized, binary I/O preserved.
- [x] Step 3.6 — `Subpoly.c` → `Subpoly.cpp`: feature flags → `constexpr`, all local allocations RAII.
- [x] Step 3.7 — `Subdb.c` → `Subdb.cpp` (continued): removed `goto`/static state, `snprintf`, converted remaining heap allocations, fixed `Open_DB`/`Close_DB` ownership.

### Phase 4 — Standalone & cleanup

- [x] Step 4.2 — Removed `.c` files from CMake, deleted `GNUmakefile`/root `Makefile`, moved sources to `src/`, removed `extern "C"` shims, consolidated `palp::min`/`palp::max`.
- [x] Step 4.4 — `#define` → `constexpr`/`using` sweep across `cws`, `SingularInput`, `Vertex`, `Subdb`, `Polynf`, `LG`, `E_Poly`/`nef`, `poly`, `MoriCone`, `Subadd`, `Global.h`, `Nef.h`, `LG.h`, `Subpoly.h`, `Rat.h`.
- [x] Step 4.5 — Fixed-size local buffers → `std::array` in `Coord.cpp` and `Vertex.cpp`.

---

## Remaining migration work

### Phase 3 — Secondary libraries (leftovers)

#### Step 3.1 — `E_Poly.cpp` (partial)

- [x] Renamed, removed local `min`/`max` macros, fixed partial array init, modernized `Die()`, converted all local `malloc`/`calloc` to RAII.
- [x] Convert `DYN_PPL.L` (`Vector *` with `realloc`) to `std::vector` growth.

#### Step 3.3 — `MoriCone.cpp` (partial)

- [x] Renamed, fixed `assert(++m < binco)` side-effect bug, converted safe heap allocations.
- [x] Replace `Inci64_*` macros with `constexpr` inline functions (performance-critical bit ops).
- [x] Remove stray `printf` debug output inside `#if TRACE_TRIANGULATION` blocks.
- [ ] Replace `exit(0)`/`exit(1)` with `palp::die()` (deferred to Phase 5).

#### Step 3.6 — `Subadd.cpp`

- [x] Renamed, compression constants modernized, dead `MOVE_SAVE_FILE` block removed.
- [x] Replace `fscanf(F, "%c%c%c%c", &A, &B, &C, &D)` with `fread` for binary I/O.
- [x] Replace `unsigned char auxUC[POLY_Dmax*VERT_Nmax]` VLA with `std::vector`.
- [x] Replace `NF_List *AuxNFLptr = NULL` global with explicit parameter passing.
- [x] Replace `fgetUI`/`fputUI` with `fread`/`fwrite`-based versions.

### Phase 4 — Standalone & cleanup

#### Step 4.1 — `lgotwist.cpp` (standalone, not in CMake)

- [x] Renamed.
- [x] Switched to C++ standard headers (`<cstdio>`, `<cstdlib>`, `<cstring>`).
- [x] Remove `register` keyword.
- [x] Replace `#define abs/min/max/mod` with inline helper functions in an anonymous namespace.
- [x] Replace `#ifdef __MSDOS__` platform sizing branch with `constexpr int` constants.
- [x] Replace small I/O/control globals with a `LgoTwistContext` struct (`infi`, `outfi`, `stdi`, `bugcount`, `invertible`).
- [ ] Remove duplicated rational arithmetic (ISSUES.md #36); reuse `Rat.cpp` logic or align types.
- **Verify**: build standalone manually.

#### Step 4.3 — Final cleanup

- [ ] Remove `-DNDEBUG` workaround in `CMakeLists.txt` (line 16) once asserts are replaced with proper error handling in Phase 5.
- [ ] Convert any remaining header-level `#define` constants that are safe to `constexpr`/`using` without breaking `#if` array-size logic.
- [x] Run full multi-dimension test sweep (POLY_Dmax = 4, 5, 6, 11): 196/196 passed for each dimension.
- [x] Run ASAN build and CTest: 196/196 passed (after fixing `MakeMobius` heap overflow).
- [x] Run UBSAN build and CTest: 196/196 passed.

#### Step 4.4 — `#define` → `constexpr`/`using` sweep (leftovers)

Symbols left as macros because they are required by `#if`/`#ifdef` array sizing or conditional compilation:

- `include/palp/Global.h`: `POLY_Dmax`, `POINT_Nmax`, `VERT_Nmax`, `FACE_Nmax`, `SYM_Nmax`, `EQUA_Nmax`, `AMBI_Dmax`, `FIB_Nmax`, `CD2F_Nmax`, `MULTIPLYING`.
- `include/palp/LG.h`: `WZinput`, `W_Nmax`.
- `include/palp/Nef.h`: `W_Nmax`.
- `include/palp/Subpoly.h`: `USE_TMP_DIR`.
- `src/Coord.cpp`: `NO_COORD_IMPROVEMENT` (always defined; removing dead blocks needs prototype cleanup).
- `src/LG.cpp`: `COEFF_Nmax` (depends on local variables, used for array sizes).
- `src/MoriCone.cpp`: `BZangle`, `SameRayBZ`, `BZR`, `BZRx`, `BZRE` (capture many local variables).
- `src/Subadd.cpp`: `TEST_UCnf`, `ADD_LIST_LENGTH`, `INCREMENTAL_TIME`, `INCREMENTAL_WRITE`, `ACCEL_PEntComp`, `USE_UNIT_ENCODE` (active, deeply interleaved).
- `src/lgotwist.cpp`: tunable constants and helper macros converted; large algorithmic arrays remain file-scope globals for now.
- `include/palp/Global.h`: large-`VERT_Nmax` `INCI_*` branch (currently dead since `VERT_Nmax <= 64`).

---

## Phase 5 — Bug fixes (after migration)

Each fix is a separate step with a regression test where possible. See `ISSUES.md` for details.

- [x] **5.1** Fix `assert(++m < binco)` side effect in `MoriCone.cpp` (no remaining `assert` with `++`/`--` side effects in MoriCone; other files still have some).
- [x] **5.2** Replace `exit(0)` with `exit(1)` on all error paths.
- [x] **5.3** Fix memory leaks in `Fano5d` (`Polynf.cpp`).
- [x] **5.3a** Fix memory leaks in `mori.cpp` main loop (`CW`, `E`, `DE`, `_P`, `_DP`; `PM`/`DPM` are stack arrays).
- [x] **5.4** Fix `realloc` losing old pointer in `E_Poly.cpp` `DYNadd_for_completion`.
- [x] **5.5** Fix partial array initialization in `E_Poly.cpp`.
- [ ] **5.6** Add missing `fscanf` return-value checks everywhere.
- [x] **5.7** Fix `system()` command injection in `SingularInput.cpp`.
- [x] **5.8** Fix buffer-overflow risks (`LG.cpp` `MakeMobius` split packed buffer; `lgotwist.cpp`, `cws.cpp`, `Subdb.cpp` no ASAN failures in current tests).
- [x] **5.9** Fix integer-overflow risks (`Polynf.cpp` volume product checks, `lgotwist.cpp` `Lcm`).
- [x] **5.10** Replace critical `assert`s used as control flow with explicit `if` checks (`Vertex.cpp`).
- [ ] **5.11** Replace global `inFILE`/`outFILE` state with a `PalpContext` struct.
- [x] **5.12** Remove dead code (`FileRW`, `OLD_code`, debug prints, `TEST` toggles).
- [x] **5.13** Fix remaining `assert`-as-control-flow bugs (`Polynf.cpp` division-by-zero, `Vertex.cpp` Euler check, `MoriCone.cpp` Inci64 limit, `LG.cpp` phase GCD checks).
- [x] **5.14** Replace critical `assert`s that guard mathematical correctness with hard errors (`Vertex.cpp`, `MoriCone.cpp`).

---

## File dependency graph

```text
Rat.c           → Step 1.1
Vertex.c        → Step 1.2
Coord.c         → Step 1.3
Polynf.c        → Step 1.4
LG.c            → Step 1.5

poly.c          → Step 2.1  (uses LG, Global)
cws.c           → Step 2.2  (uses LG, Global)
class.c         → Step 2.3  (uses Subpoly)
nef.c           → Step 2.4  (uses Nef, LG)
mori.c          → Step 2.5  (uses Mori, LG)

E_Poly.c        → Step 3.1  (uses Nef, Rat)
Nefpart.c       → Step 3.2  (uses Nef)
MoriCone.c      → Step 3.3  (uses Rat, Mori)
SingularInput.c → Step 3.4  (uses Mori)
Subdb.c         → Step 3.5  (uses Subpoly, Rat)
Subpoly.c       → Step 3.6  (uses Rat, Subpoly.h)
Subadd.c        → Step 3.6  (uses Subpoly)
Subdb.c         → Step 3.7  (continued)

lgotwist.c      → Step 4.1  (standalone)
```

## Risk assessment

| Step | Risk | Reason |
|------|------|--------|
| 1.4 `Polynf` | High | Largest file (~3223 lines), complex math, many conditionals |
| 2.2 `cws` | High | 1935 lines, many code paths, hardcoded tables |
| 3.7 `Subdb` | Medium | File I/O, database logic, ownership semantics |
| 1.5 `LG` | Medium | Polynomial arithmetic, manual memory management |
| 3.3 `MoriCone` | Medium | Complex geometry, assert side-effect bug |
| All others | Low–Medium | Standard C-to-C++ conversion |
