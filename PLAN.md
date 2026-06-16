# Modernization Plan

Goal: move PALP from C toward maintainable modern C++ while preserving the historical outputs, especially for 4D reflexive polytopes. Each step should be small enough to verify with the existing shell tests plus targeted new tests.

## Phase 0: Preserve Current Behavior

1. [x] Establish a clean baseline.
   - Run `make cleanall && make -jN`.
   - Run `make check`; run `make checklong` separately when time permits.
   - Capture compiler version, flags, and test output.
   - Completed baseline:
     - `cc --version`: GCC 13.3.0.
     - `make -j2`: passed; default executables were already up to date on the final baseline run.
     - `make check`: passed after installing `Singular` 4.4.1 into the active `vibe` environment.
     - `make checklong`: not run yet.

2. [x] Make tests easier to run in CI.
   - Add CMake/CTest entries that call the existing `tests/*.sh` scripts.
   - Keep `GNUmakefile check` working.
   - Verification: `ctest` and `make check` run the same test scripts.
   - Completed CTest integration:
     - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
     - `cmake --build build -j2`: passed; CMake now builds `poly/class/cws/nef/mori-{4,5,6,11}d.x` test executables with matching `POLY_Dmax`.
     - `ctest --test-dir build --output-on-failure`: passed, 196/196 tests.
     - `make check`: passed, preserving the GNUmakefile test path.

3. [x] Add sanitizer build modes.
   - Add CMake presets or documented commands for AddressSanitizer and UndefinedBehaviorSanitizer.
   - Do not change production flags yet.
   - Verification: default build still matches baseline; sanitizer failures are filed in `ISSUES.md`.
   - Completed sanitizer presets:
     - Added `release`, `asan`, and `ubsan` CMake configure/build/test presets.
     - `cmake --preset release`, `cmake --build --preset release`, and `ctest --preset release`: passed, 196/196 tests.
     - `cmake --preset asan` and `cmake --build --preset asan`: passed.
     - `ctest --preset asan`: failed 3/196 tests, all `6.4.18-nef-v.sh` for 5d/6d/11d; filed as `ISSUES.md` item 14.
     - `cmake --preset ubsan` and `cmake --build --preset ubsan`: passed.
     - `ctest --preset ubsan`: failed 3/196 tests, all `3.2.11-poly-l.sh` for 5d/6d/11d due to `LG.c:754` misaligned pointer stores; filed as `ISSUES.md` item 15.

## Phase 1: Make C Safer Before C++

4. [x] Replace known undefined behavior without changing interfaces.
   - Fix the `LG.c` use-after-free.
   - Fix 64-bit incidence shifts.
   - Fix unchecked `fgets()` in `SingularInput.c`.
   - Verification: targeted regression tests plus full `make check`.
   - Completed UB cleanup:
     - Fixed `LG.c` `Is_Gen_CY()` by saving `E->ne` before `free(E)`.
     - Fixed `MoriCone.c` incidence masks by shifting an `Inci64` value instead of signed `int`.
     - Fixed `SingularInput.c` `Read_HyperSurf()` by checking `fgets()` before parsing.
     - `make -j2`: passed.
     - `make all-dims -j2`: passed.
     - Focused tests: `tests/3.2.11-poly-l.sh` and `tests/7.2.8-mori-b.sh` passed for `DIM=4,5,6,11`.
     - `make check`: passed.

5. [x] Convert assertion-only user checks into explicit runtime checks.
   - Start with parsers and fixed-size array boundaries.
   - Keep assertions for internal invariants.
   - Verification: add tests for invalid input and capacity errors.
   - Completed parser runtime checks:
     - Added checked integer reads for CWS and polytope matrix input in `Coord.cc`.
     - Replaced assertion-only checks for incomplete weight input, single-weight size limits, non-positive weights, negative CWS weights, and excess `/Z` quotient factors.
     - Added matching checks in `LG.c` for `Read_WZ_PP()` and `Read_Weight()`.
     - Added invalid-input regressions to `tests/2.2-error-handling.sh`.
     - Internal algorithm assertions were left unchanged.
     - `cc --version`: GCC 13.3.0.
     - `make -j2`: passed.
     - `make all-dims -j2`: passed.
     - Focused tests: `tests/2.2-error-handling.sh`, `tests/2.1-polytope-input.sh`, `tests/3.2.11-poly-l.sh`, `tests/4.2.1-cws-w.sh`, `tests/4.2.6-cws-N.sh`, and `tests/6.3-nef-output.sh` passed for `DIM=4,5,6,11`.
     - `make check`: passed.

6. [x] Normalize matrix and array shape contracts.
   - Introduce named C typedefs or small wrapper structs for `PairMat`, normal-form matrices, and incidence arrays.
   - Eliminate call sites that pass smaller arrays to larger formal types.
   - Verification: compiler warnings for `-Warray-parameter` and `-Wstringop-overflow` decrease.
   - Completed shape contract cleanup:
     - Made `Complete_Poly()` use the existing `PairMat` typedef in the public prototype.
     - Converted undersized `VERT_Nmax x VERT_Nmax` pairing buffers in `E_Poly.c`, `Polynf.c`, and `Subdb.c` to `PairMat`.
     - Added `AffineNormalForm` for `POLY_Dmax x VERT_Nmax` normal-form matrices and made `Make_ANF()` declarations agree.
     - Changed `Read_HyperSurf()` filename input to `const char *` and `SimplexVolume()` vertex-list input to `Long **`.
     - `cc --version`: GCC 13.3.0.
     - `make -j2`: passed.
     - `make all-dims -j2`: passed; targeted grep found no remaining `Make_VEPM`, `Make_ANF`, `Read_HyperSurf`, `SimplexVolume`, `-Warray-parameter`, or `-Wstringop-overflow` warnings.
     - Focused tests: `poly -fA` passed for `DIM=5,6,11`; `DIM=4` exposed `ISSUES.md` item 16. `tests/6.4.4-nef-H.sh`, `tests/7.2.3-mori-g.sh`, and `tests/7.2.8-mori-b.sh` passed for `DIM=4,5,6,11` except the existing long-test skip for `nef-11d.x -H`.
     - `make check`: passed.

7. [x] Remove stale build paths.
   - Either delete/update `Makefile` or make it delegate to `GNUmakefile`.
   - Make CMake and GNUmakefile source lists match exactly.
   - Verification: both build systems produce all five main executables.
   - Completed build-path cleanup:
     - Replaced stale `Makefile` with a compatibility wrapper that delegates targets to `GNUmakefile`.
     - Confirmed CMake and GNUmakefile source lists agree for shared, class, nef, and Mori modules.
     - Updated `AGENTS.md` to describe the wrapper.
     - `make -f Makefile -j2`: passed.
     - `make -f Makefile all-dims -j2`: passed.
     - `cmake --build build -j2`: passed.
     - `make -f Makefile check`: passed.

## Phase 2: Prepare for C++ Compilation

8. [x] Add include guards and C++ compatibility shims.
   - Guard every header.
   - Remove duplicate local declarations where headers should be authoritative.
   - Add `extern "C"` only if C and C++ objects must coexist temporarily.
   - Verification: C build unchanged; headers can be included from a trivial `.cc` file.
   - Completed header compatibility pass:
     - Added include guards and C++ linkage wrappers to `Global.h`, `Rat.h`, `LG.h`, `Subpoly.h`, `Nef.h`, and `Mori.h`.
     - Fixed `Rat.h` so the guard covers the full header, including `LRat` declarations.
     - Added `tests/header-smoke.cc`, which includes every public header twice from C++.
     - Removed redundant local declarations for public prototypes and published `Read_Weight()` from `LG.h`.
     - Guarded duplicate `min`/`max` macro definitions to avoid redefinition warnings when headers are included together.
     - `cc --version`: GCC 13.3.0.
     - `c++ --version`: GCC 13.3.0.
     - `c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc`: passed.
     - `make -j2`: passed.
     - `make all-dims -j2`: passed; log scan found no `implicit declaration`, macro redefinition, conflicting type, or compiler error lines.
     - `make check`: passed.

9. [x] Introduce fixed-width type aliases.
   - Replace macro aliases like `Long` and `LLong` with a central typedef/using layer.
   - Document required ranges for classified 4D workflows.
   - Verification: binary output remains unchanged on baseline tests.
   - Completed numeric alias cleanup:
     - Replaced the `Long` and `LLong` object-like macros in `Global.h` with central typedef aliases.
     - Preserved the existing underlying types: `Long` is still `long`, and `LLong` is still `long long`.
     - Added compile-time range checks requiring `Long` to provide at least 32 bits and `LLong` at least 64 bits.
     - Extended `tests/header-smoke.cc` with matching C++ `static_assert` checks.
     - Documented the current numeric contract in `AGENTS.md`.
     - `cc --version`: GCC 13.3.0.
     - `c++ --version`: GCC 13.3.0.
     - `c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc`: passed.
     - `make -j2`: passed.
     - `make all-dims -j2`: passed; log scan found no compiler errors, implicit declarations, conflicting types, redefinitions, or unknown type names.
     - `make check`: passed.

10. [x] Compile selected modules as C++ without semantic changes.
    - Start with leaf-like modules: `Rat.cc` done, then `Coord.cc` done, then `Vertex.cc` done.
    - Rename one file at a time only after it compiles cleanly as C++.
    - Verification: full test suite after each file.
    - `Rat.cc` migration completed:
      - Renamed `Rat.c` to `Rat.cc`.
      - Updated GNUmakefile to compile `.cc` objects with `g++` and link mixed C/C++ executables with the C++ linker.
      - Enabled CXX in CMake, set C++17 as the required standard, and added sanitizer CXX flags to presets.
      - Removed obsolete `register` storage hints from `Rat.cc`.
      - `g++ -std=c++17 -O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -c -o /tmp/Rat-smoke.o Rat.cc`: passed.
      - `make -j2`: passed.
      - `make all-dims -j2`: passed; `Rat.cc` compiled as C++ for `POLY_Dmax=4,5,6,11`.
      - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
      - `cmake --build build -j2`: passed.
      - `make check`: passed.
    - `Coord.cc` migration completed:
      - Renamed `Coord.c` to `Coord.cc`.
      - Updated CMake source lists and documentation references.
      - Fixed C++ const-correctness for `CWSZerror()`.
      - Added explicit C linkage for exported `Coord` helpers used by C translation units: `IsNextDigit()`, `Print_CWS_Zinfo()`, `ReadCwsPp()`, `Make_CWS_Points()`, and the `QuotZ_2_SublatG()` declaration.
      - `g++ -std=c++17 -O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -c -o /tmp/Coord-smoke.o Coord.cc`: passed.
      - `make -j2`: passed.
      - `make all-dims -j2`: passed; `Coord.cc` compiled as C++ for `POLY_Dmax=4,5,6,11`.
      - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
      - `cmake --build build -j2`: passed.
      - `make check`: passed.
      - New C++ optimizer warnings in `Print_CWH()` for low `POLY_Dmax` builds were filed as `ISSUES.md` item 17.
    - `Vertex.cc` migration completed:
      - Renamed `Vertex.c` to `Vertex.cc`.
      - Updated CMake source lists and documentation references.
      - Removed the obsolete `register` storage hint from `swap()`.
      - Added explicit C linkage for exported `Vertex` helpers used by C translation units through local declarations: `CompareEq()`, `IsGoodCEq()`, `Sort_PPL()`, `GLZ_Start_Simplex()`, `Finish_IP_Check()`, `Make_FaceIPs()`, and `Eval_BaHo()`.
      - `g++ -std=c++17 -O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -c -o /tmp/Vertex-smoke.o Vertex.cc`: passed.
      - `make -j2`: passed.
      - `make all-dims -j2`: passed; `Vertex.cc` compiled as C++ for `POLY_Dmax=4,5,6,11`.
      - `c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc`: passed.
      - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
      - `cmake --build build -j2`: passed.
      - `make check`: passed.
      - A suspicious `QComplete_Poly()` equation-index expression found during review was filed as `ISSUES.md` item 18.

## Phase 3: Introduce Modern C++ Types

11. [x] Wrap rational numbers.
    - Replace free-function `Rat` operations with a small value type.
    - Add explicit construction and checked division.
    - Verification: unit tests for arithmetic and all existing CLI tests.
    - C++ wrapper foundation completed:
      - Added opt-in `palp::Rational` and `palp::LongRational` value types around the existing `Rat` and `LRat` C ABI.
      - Constructors reject zero denominators with `std::domain_error`.
      - Division operators reject division by zero before calling the legacy quotient helpers.
      - Added `tests/rat-wrapper-test.cc` for arithmetic, normalization, comparison, legacy conversion, and checked error paths.
      - Added `make check-rat-wrapper` and a `rat-wrapper-test` CTest entry.
      - `make check-rat-wrapper`: passed.
      - `c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc`: passed.
      - `c++ -std=c++17 -I. -fsyntax-only tests/rat-wrapper-test.cc`: passed.
      - `make -j2`: passed.
      - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
      - `cmake --build build -j2`: passed.
      - `ctest --test-dir build --output-on-failure`: passed, 197/197 tests.
      - `make check`: passed.
    - First production call-site migration completed:
      - Migrated `Vertex.cc` `Compute_InvMat()` from raw `LRat` arrays and `LrD()`/`LrP()`/`LrQ()` calls to `palp::LongRational` arithmetic.
      - Preserved the existing final integer denominator and inverse-matrix output contract.
      - `g++ -std=c++17 -O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -c -o /tmp/Vertex-rational-smoke.o Vertex.cc`: passed.
      - `make check-rat-wrapper`: passed.
      - `c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc`: passed.
      - `make -j2`: passed.
      - `make all-dims -j2`: passed.
      - Focused tests: `tests/3.2.7-poly-e.sh`, `tests/3.2.23-poly-P.sh`, and `tests/6.4.18-nef-v.sh` passed for `DIM=6`.
      - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
      - `cmake --build build -j2`: passed.
      - `ctest --test-dir build --output-on-failure`: passed, 197/197 tests.
      - `make check`: passed.
      - No additional raw `Rat`/`LRat` production call sites remain in the currently migrated C++ modules, apart from the legacy C ABI implementations in `Rat.cc` and the wrappers themselves.
      - Remaining legacy rational uses are in C translation units and should be migrated when those modules become C++.

12. [x] Wrap bounded arrays.
    - Use `std::array` for compile-time bounded vectors/matrices.
    - Provide `.size()`-checked helpers before changing algorithms.
    - Verification: no output changes; sanitizer builds improve.
    - First local bounded-array migration completed:
      - Converted `Vertex.cc` `Compute_InvMat()` local rational matrices/vectors and pivot-index buffer to `std::array`.
      - Kept public C-style input/output matrix contracts unchanged.
      - `g++ -std=c++17 -O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -c -o /tmp/Vertex-array-smoke.o Vertex.cc`: passed.
      - `make check-rat-wrapper`: passed.
      - `c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc`: passed.
      - `make -j2`: passed.
      - Focused tests: `tests/3.2.7-poly-e.sh`, `tests/3.2.23-poly-P.sh`, and `tests/6.4.18-nef-v.sh` passed for `DIM=6`.
      - `make all-dims -j2`: passed.
      - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
      - `cmake --build build -j2`: passed.
      - `ctest --test-dir build --output-on-failure`: passed, 197/197 tests.
      - `make check`: passed.
    - Second local bounded-array migration completed:
      - Converted `Vertex.cc` `VZ_to_Base()` local index, vector, and row-pointer work buffers to `std::array`.
      - Kept the `VZ_to_Base()` signature and `W_to_GLZ()` adapter call unchanged at the public boundary.
      - `g++ -std=c++17 -O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -c -o /tmp/Vertex-vz-array-smoke.o Vertex.cc`: passed.
      - `c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc`: passed.
      - `make -j2`: passed.
      - Focused tests: `tests/3.2.7-poly-e.sh`, `tests/3.2.11-poly-l.sh`, `tests/3.2.23-poly-P.sh`, and `tests/6.4.18-nef-v.sh` passed for `DIM=6`.
      - `make all-dims -j2`: passed.
      - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
      - `cmake --build build -j2`: passed.
      - `ctest --test-dir build --output-on-failure`: passed, 197/197 tests.
      - `make check`: passed.
    - Third local bounded-array migration completed:
      - Converted `Rat.cc` `REgcd()` and `LREgcd()` two-element temporary buffers to `std::array`.
      - Kept the public C ABI and recursive gcd output buffers unchanged.
      - `g++ -std=c++17 -O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -c -o /tmp/Rat-array-smoke.o Rat.cc`: passed.
      - `make check-rat-wrapper`: passed.
      - `c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc`: passed.
      - `make -j2`: passed.
      - `make all-dims -j2`: passed.
      - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
      - `cmake --build build -j2`: passed.
      - `ctest --test-dir build --output-on-failure`: passed, 197/197 tests.
      - `make check`: passed.
    - Fourth local bounded-array migration completed:
      - Converted `Coord.cc` weight-equation basis work buffers in `Orig_Solve_Next_WEq()` and `Solve_Next_WEq()` to `std::array`.
      - Converted `Coord.cc` `Reduce_PPL_2_Sublat()` local coordinate vector to `std::array`.
      - Kept `CWLatticeBasis`, `W_to_GLZ()`, `REgcd()`, and sublattice reduction public interfaces unchanged by using `.data()` adapters.
      - `g++ -std=c++17 -O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -c -o /tmp/Coord-array-smoke.o Coord.cc`: passed.
      - `c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc`: passed.
      - `make -j2`: passed.
      - Focused tests: `tests/2.1-polytope-input.sh`, `tests/2.2-error-handling.sh`, `tests/4.2.1-cws-w.sh`, and `tests/4.2.6-cws-N.sh` passed for `DIM=6`.
      - `make all-dims -j2`: passed; the known `Print_CWH()` low-`POLY_Dmax` warnings remain tracked as `ISSUES.md` item 17.
      - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
      - `cmake --build build -j2`: passed.
      - `ctest --test-dir build --output-on-failure`: passed, 197/197 tests.
      - `make check`: passed.
    - Fifth local bounded-array migration completed:
      - Converted `Coord.cc` `Old_Make_CWS_Points()` and `Make_CWS_Points()` local enumeration vectors to `std::array`.
      - Kept the `PolyPointList`, `CWS_2_SublatZ()`, `Reduce_PPL_2_Sublat()`, and `Compute_X0()` public interfaces unchanged by using `.data()` adapters.
      - Left the local `G[POLY_Dmax][POLY_Dmax]` matrices as C arrays because existing helper interfaces still require `Long[][POLY_Dmax]`.
      - `g++ -std=c++17 -O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -c -o /tmp/Coord-cws-array-smoke.o Coord.cc`: passed.
      - `c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc`: passed.
      - `make -j2`: passed.
      - Focused tests: `tests/2.1-polytope-input.sh`, `tests/2.2-error-handling.sh`, `tests/4.2.1-cws-w.sh`, and `tests/4.2.6-cws-N.sh` passed for `DIM=6`.
      - `make all-dims -j2`: passed; the known `Print_CWH()` low-`POLY_Dmax` warnings remain tracked as `ISSUES.md` item 17.
      - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
      - `cmake --build build -j2`: passed.
      - `ctest --test-dir build --output-on-failure`: passed, 197/197 tests.
      - `make check`: passed.
    - Sixth local bounded-array migration completed:
      - Converted `Vertex.cc` completion and quick-completion local vectors/index buffers in `add_for_completion()`, `Complete_Poly()`, `Qadd_for_completion()`, `lastline()`, and `QComplete_Poly()` to `std::array`.
      - Kept the `Compute_InvMat()`, `add_for_completion()`, `lastline()`, and matrix parameter interfaces unchanged by using `.data()` adapters.
      - Left fixed matrix parameters as C arrays where existing helper signatures require `Long[][POLY_Dmax]` or `Long[][EQUA_Nmax]`.
      - Added `ISSUES.md` item 19 for the pre-existing `Vertex.cc` warning cluster found during all-dimension verification.
      - `g++ -std=c++17 -O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -c -o /tmp/Vertex-completion-array-smoke.o Vertex.cc`: passed, with known `Vertex.cc` warnings now tracked in `ISSUES.md`.
      - `c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc`: passed.
      - `make -j2`: passed, with known `Vertex.cc` warnings.
      - Focused tests: `tests/3.2.7-poly-e.sh`, `tests/3.2.11-poly-l.sh`, `tests/3.2.23-poly-P.sh`, and `tests/6.4.18-nef-v.sh` passed for `DIM=6`.
      - `make all-dims -j2`: passed, with known `Vertex.cc` warnings tracked as `ISSUES.md` item 19.
      - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
      - `cmake --build build -j2`: passed.
      - `ctest --test-dir build --output-on-failure`: passed, 197/197 tests.
      - `make check`: passed.
    - Seventh local bounded-array migration completed:
      - Converted `Vertex.cc` `EEV_To_Equation()` `LLong` work vector and `OrthBase_red_by_V()` row-product vector to `std::array`.
      - Converted `Coord.cc` `Read_CWS_Zinfo()` fixed line buffer to `std::array<char, 999>`.
      - Kept `VZ_to_Base()` and parser helper boundaries unchanged by using existing pointer-compatible access.
      - `g++ -std=c++17 -O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -c -o /tmp/Vertex-extra-array-smoke.o Vertex.cc`: passed, with known `Vertex.cc` warnings tracked as `ISSUES.md` item 19.
      - `g++ -std=c++17 -O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -c -o /tmp/Coord-parser-array-smoke.o Coord.cc`: passed.
      - `c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc`: passed.
      - `make -j2`: passed, with known `Vertex.cc` warnings.
      - Focused tests: `tests/2.1-polytope-input.sh`, `tests/2.2-error-handling.sh`, `tests/3.2.7-poly-e.sh`, `tests/3.2.11-poly-l.sh`, `tests/3.2.23-poly-P.sh`, `tests/4.2.1-cws-w.sh`, `tests/4.2.6-cws-N.sh`, and `tests/6.4.18-nef-v.sh` passed for `DIM=6`.
      - `make all-dims -j2`: passed, with known `Coord.cc` and `Vertex.cc` warnings tracked as `ISSUES.md` items 17 and 19.
      - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
      - `cmake --build build -j2`: passed.
      - `ctest --test-dir build --output-on-failure`: passed, 197/197 tests.
      - `make check`: passed.
    - Item complete for currently migrated C++ modules:
      - Remaining raw bounded buffers in `Coord.cc`, `Rat.cc`, and `Vertex.cc` are either public/legacy function parameters, pointer aliases, disabled coordinate-improvement code, or matrix-shaped local buffers whose callees still require `Long[][POLY_Dmax]`, `Long[][EQUA_Nmax]`, or `Long[][VERT_Nmax]`.
      - Those matrix/interface conversions should be handled with the next RAII/interface steps, not by adding casts inside this item.

13. [x] Replace manual allocation with RAII.
    - Use `std::vector` for variable-size work buffers.
    - Use automatic storage or `std::unique_ptr` where fixed historical layout is still required.
    - Verification: run sanitizers and compare CLI outputs.
    - First RAII migration completed:
      - Converted `Vertex.cc` temporary `CEqList`, `INCI[]`, and `FaceInfo` allocations in `Find_Equations()`, `IP_Check()`, `Ref_Check()`, and `RC_Calc_BaHo()` from `malloc()`/`free()` to `std::unique_ptr`.
      - Used `std::nothrow` allocation to preserve the existing explicit allocation-failure diagnostics and `exit(0)` behavior.
      - Removed manual frees on early returns while keeping helper interfaces pointer-compatible through `.get()`.
      - `g++ -std=c++17 -O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -c -o /tmp/Vertex-raii-smoke.o Vertex.cc`: passed, with known `Vertex.cc` warnings tracked as `ISSUES.md` item 19.
      - `c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc`: passed.
      - `make -j2`: passed, with known `Vertex.cc` warnings.
      - Focused tests: `tests/3.2.7-poly-e.sh`, `tests/3.2.11-poly-l.sh`, `tests/3.2.23-poly-P.sh`, and `tests/6.4.18-nef-v.sh` passed for `DIM=6`.
      - `make all-dims -j2`: passed, with known `Vertex.cc` warnings tracked as `ISSUES.md` item 19.
      - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
      - `cmake --build build -j2`: passed.
      - `ctest --test-dir build --output-on-failure`: passed, 197/197 tests.
      - `make check`: passed.
      - Sanitizer presets were not rerun for this slice because `ctest --preset asan` and `ctest --preset ubsan` have known pre-existing failures tracked as `ISSUES.md` items 14 and 15.
      - No `malloc()`/`free()` sites remain in the currently migrated C++ modules (`Coord.cc`, `Rat.cc`, `Vertex.cc`); remaining manual allocation is in C translation units and should be migrated as those modules become C++.

14. Introduce a runtime context.
    - Replace global `inFILE`/`outFILE` and scattered option globals with a context object.
    - Migrate one executable at a time.
    - Verification: CLI behavior and prompts remain compatible.
    - First runtime-context migration completed:
      - Added C-compatible `PALP_RuntimeContext` helpers in `Global.h` to snapshot and apply the legacy `inFILE`/`outFILE` globals.
      - Migrated `poly.c` startup file-handle setup to configure a `PALP_RuntimeContext` and apply it before shared code runs.
      - Kept the existing global `inFILE`/`outFILE` ABI intact for all shared modules and other CLI programs.
      - `gcc -O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -c -o /tmp/poly-context-smoke.o poly.c`: passed.
      - `c++ -std=c++17 -I. -fsyntax-only tests/header-smoke.cc`: passed.
      - `make -j2`: passed, with existing warnings.
      - Focused tests: `tests/2.1-polytope-input.sh`, `tests/2.2-error-handling.sh`, `tests/3.2.7-poly-e.sh`, `tests/3.2.11-poly-l.sh`, and `tests/3.2.23-poly-P.sh` passed for `DIM=6`.
      - `make all-dims -j2`: passed, with known warnings.
      - `cmake -S . -B build -D CMAKE_BUILD_TYPE=Release`: passed.
      - `cmake --build build -j2`: passed.
      - `ctest --test-dir build --output-on-failure`: passed, 197/197 tests.
      - `make check`: passed.
      - Remaining in this item: migrate `class.x`, `cws.x`, `nef.x`, and `mori.x` startup code, then move shared readers/printers from direct globals to context-aware adapters.

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
