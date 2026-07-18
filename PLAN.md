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

- [x] `CMakeLists.txt`: CXX-only, C++17, `-Wall -Wextra`, `palp_sources()` helper.
- [x] Headers moved to `include/palp/` with `#pragma once`; include paths updated.
- [x] `palp_types.h` shim added and integrated via `Global.h`.

### Phase 1 — Foundation libraries

- [x] `Rat.c` → `Rat.cpp`.
- [x] `Vertex.c` → `Vertex.cpp`.
- [x] `Coord.c` → `Coord.cpp`.
- [x] `Polynf.c` → `Polynf.cpp`.
- [x] `LG.c` → `LG.cpp` (including `MakeMobius` split-buffer fix).

### Phase 2 — Driver programs

- [x] `poly.c` → `poly.cpp`.
- [x] `cws.c` → `cws.cpp`.
- [x] `class.c` → `class.cpp`.
- [x] `nef.c` → `nef.cpp`.
- [x] `mori.c` → `mori.cpp`.

### Phase 3 — Secondary libraries

- [x] `E_Poly.c` → `E_Poly.cpp`.
- [x] `Nefpart.c` → `Nefpart.cpp`.
- [x] `MoriCone.c` → `MoriCone.cpp`.
- [x] `SingularInput.c` → `SingularInput.cpp`.
- [x] `Subdb.c` → `Subdb.cpp`.
- [x] `Subpoly.c` → `Subpoly.cpp`.
- [x] `Subadd.c` → `Subadd.cpp`.

### Phase 4 — Standalone & cleanup

- [x] Removed `.c` files from CMake, deleted legacy Makefiles, moved sources to `src/`,
  removed `extern "C"` shims, consolidated `palp::min`/`palp::max`.
- [x] `#define` → `constexpr`/`using` sweep across all files where safe.
- [x] Fixed-size local buffers → `std::array`/`std::vector` where applicable.
- [x] `lgotwist.c` → `lgotwist.cpp`: standard headers, `constexpr` constants,
  `LgoTwistContext`, local rational helpers deduplicated (partial; full `Rat.cpp`
  sharing deferred to Phase 5.11).
- [x] Removed `-DNDEBUG` workaround from `CMakeLists.txt`.
- [x] Release, ASAN, and UBSAN builds pass 196/196 tests for POLY_Dmax = 4, 5, 6, 11.

### Phase 5 — Bug fixes

- [x] **5.1** `assert(++m < binco)` side effect moved out of `assert()` in `MoriCone.cpp`.
- [x] **5.2** `exit(0)` on error paths replaced with `exit(1)` across the codebase.
- [x] **5.3/5.3a** Memory leaks in `Fano5d` (`Polynf.cpp`) and `mori.cpp` main loop fixed via RAII.
- [x] **5.4** `realloc` losing old pointer in `E_Poly.cpp` resolved by `std::vector`.
- [x] **5.5** Partial array initialization in `E_Poly.cpp` fixed.
- [x] **5.6** `fscanf` return-value checks added everywhere except `nef.cpp` `Read_WPCICY`
  (checked there would change output).
- [x] **5.7** `system()` command injection in `SingularInput.cpp` replaced with `fork/execvp`.
- [x] **5.8** Buffer-overflow risks fixed (`LG.cpp` `MakeMobius`, `lgotwist.cpp` skeleton bounds,
  `cws.cpp`/`Subdb.cpp` buffers).
- [x] **5.9** Integer-overflow risks fixed (`Polynf.cpp` volume products, `lgotwist.cpp` `Lcm`).
- [x] **5.10/5.13/5.14** Critical `assert()`s used as control flow or mathematical guards
  converted to explicit `fputs(..., stderr); exit(1);` errors in all active files.
- [x] **5.12** Dead code removed (`FileRW`, `OLD_code`, debug prints, `TEST` toggles).

---

## Remaining work

### Phase 5 — Bug fixes (remaining)

- [x] **5.11** Replace global `inFILE`/`outFILE` state with explicit `FILE *in`/`FILE *out` parameters.
  Completed in five granular commits (5.11.1–5.11.5):

  - **5.11.1** Introduced `PalpContext` in `Global.h` and threaded explicit
    `FILE *out` through `Rat.cpp` printing helpers, adding backward-compatible
    overloads for global `outFILE`.
  - **5.11.2** Converted pure-output helpers across `E_Poly.cpp`, `Nefpart.cpp`,
    `Polynf.cpp`, `Subadd.cpp`, `Subdb.cpp`, `Subpoly.cpp`, `Vertex.cpp`, and
    internal helpers in `cws.cpp`/`nef.cpp`/`MoriCone.cpp`/`LG.cpp` to take
    explicit `FILE *out` parameters.
  - **5.11.3** Converted pure input helpers in `Coord.cpp`, `LG.cpp`,
    `MoriCone.cpp`, `nef.cpp`, `Subdb.cpp`, `Subpoly.cpp`, and `Polynf.cpp` to
    take explicit `FILE *in` parameters. Filter-mode (`inFILE = NULL`) was
    replaced with local `in = NULL`/`in = stdin` handling without global state.
  - **5.11.4** Removed the global `outFILE` declaration and updated all driver
    `main()` files (`poly.cpp`, `class.cpp`, `mori.cpp`, `nef.cpp`, `cws.cpp`)
    to own local `FILE *out`.
  - **5.11.5** Removed the global `inFILE` declaration and updated all driver
    `main()` files to own local `FILE *in`, passing it through all remaining
    input paths.

  All 196 CTest tests pass for release and ASan builds. The only remaining
  `inFILE`/`outFILE` references are inside commented-out code blocks.

### Phase 4 — Final cleanup (remaining)

- [ ] Convert any remaining header-level `#define` constants that are safe to
  `constexpr`/`using` without breaking `#if` array-size logic. Currently kept as
  macros: `POLY_Dmax`, `POINT_Nmax`, `VERT_Nmax`, `FACE_Nmax`, `SYM_Nmax`,
  `EQUA_Nmax`, `AMBI_Dmax`, `FIB_Nmax`, `CD2F_Nmax`, `MULTIPLYING` (`Global.h`);
  `WZinput`, `W_Nmax` (`LG.h`/`Nef.h`); `USE_TMP_DIR` (`Subpoly.h`); local
  conditional flags in `Coord.cpp`, `LG.cpp`, `MoriCone.cpp`, `Subadd.cpp`, and
  `lgotwist.cpp`.
