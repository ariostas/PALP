#PALP Potential Issues / Bugs

Found during code inspection for the C → C++17 migration. Each item has a file,
line number, severity, and description. Fixes are applied in Phase 5 (after
migration), one regression-tested fix at a time.

Line numbers refer to the original C source as of the start of the
modernization effort. They will shift as migration proceeds; the function
names and descriptions should remain locatable.

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
  - `#37` Duplicated rational arithmetic in `lgotwist.cpp` partially deduplicated
    (duplicate `gcd` removed; full `Rat.cpp` sharing deferred until Phase 5.11).
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
- **Status**: Open
- **Description**: `static int MaxPoNum` and `static int M` — persistent state
  across calls.

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
- **Status**: Open
- **Description**: `NF_List *AuxNFLptr = NULL` — global mutable pointer used as
  a "dirty trick" for statistics. Hidden coupling between functions.

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
- **File**: `Polynf.cpp`
- **Lines**: 348, 470 (original)
- **Severity**: Low
- **Status**: Open
- **Description**: `volatile` used to prevent compiler optimization. Should be
  replaced with `std::atomic`, memory barriers, or the root cause addressed.
