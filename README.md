# PALP: a Package for Analyzing Lattice Polytopes

This repository is a fork of the [PALP GitLab repository](https://gitlab.com/stringstuwien/PALP),
maintained for a C to C++20 modernization of the codebase. The original package was
written by Maximilian Kreuzer and Harald Skarke. More information about PALP can be
found at the [PALP website](http://hep.itp.tuwien.ac.at/~kreuzer/CY/CYpalp.html)
and in the [PALP online documentation](http://palp.itp.tuwien.ac.at/wiki/index.php/PALP_online_documentation).

## Building

PALP is built with CMake. To configure a Release build:

```bash
cmake -S . -B build -D CMAKE_BUILD_TYPE=Release
```

Useful CMake options:

```bash
-D CMAKE_INSTALL_PREFIX=[your custom path]      # default system install prefix
-D POLY_Dmax=[dimension]                        # default: 6
-D PALP_TEST_DIMENSIONS="4;5;6;11"              # dimensions for which to build executables and run tests
-D PALP_TESTS_LONG=ON                           # include long-running tests (default: OFF)
```

`POLY_Dmax` controls the maximum dimension the executables are compiled for.
Because many internal arrays are sized from this constant, larger values use
more memory. The build produces dimension-suffixed executables such as
`poly-4d.x`, `poly-5d.x`, ..., `poly-11d.x` as well as unsuffixed copies
compiled with `POLY_Dmax`.

**Note for the 4D reflexive classification:** `class.x` must be compiled with
`POLY_Dmax=4` for that task. Use `class-4d.x` (or build with
`-D POLY_Dmax=4 -D PALP_TEST_DIMENSIONS=4`).

Build the package with:

```bash
cmake --build build
```

Install the binaries with:

```bash
cmake --build build --target install
```

## Running the tests

After configuring, run the CTest suite with:

```bash
ctest --test-dir build --output-on-failure
```

By default this runs 196 tests across the configured dimensions
(`PALP_TEST_DIMENSIONS`). To also run long tests, configure with
`-D PALP_TESTS_LONG=ON`.

## Programs and basic usage

PALP provides several command-line programs. Each program accepts `-h` for
help and can be used as a filter (`-f` or `-`) when listed as the first option.

| Program | Purpose |
|-------- | ------- |
| `poly.x` | Compute polytope data: vertices, facets, dual points, Hodge numbers, symmetries, normal forms, etc. |
| `cws.x`  | Generate and analyze weight systems and combined weight systems (CWS). |
| `class.x`| Classify reflexive polytopes and work with binary lists/databases of normal forms. |
| `nef.x`  | Compute Hodge numbers of nef-partitions of reflexive polytopes. |
| `mori.x` | Compute star triangulations, Mori cones, and intersection rings. |
| `lgotwist.x` | Compute twisted sectors for Landau–Ginzburg models. |

All programs read input either as a weight system on one line, e.g.

```text
5 1 1 1 1 1
```

or as a matrix of lattice points, e.g.

```text
3 2
2 0
0 2
0 0
```

The first form gives `dimension #points` followed by the coordinates.
Rows and columns may be transposed.

### Quick example with `poly.x`

Analyze the quintic hypersurface polytope:

```bash
echo '5 1 1 1 1 1' | ./poly-4d.x -g
```

Output:

```text
5 1 1 1 1 1 M:126 5 N:6 5 H:1,101 [-200]
```

This prints the numbers of lattice points of the M- and N-lattice polytopes,
their vertices, and the Hodge numbers of the corresponding Calabi–Yau
hypersurface.

## Reproducing the classification of 4D reflexive polytopes

The famous Kreuzer–Skarke list of 473,800,776 four-dimensional reflexive
polytopes can be reproduced with PALP's `cws.x` and `class.x` programs. The
algorithm is described in
[hep-th/9805190](https://arxiv.org/abs/hep-th/9805190) and
[hep-th/0002240](https://arxiv.org/abs/hep-th/0002240); a pedagogical
walkthrough for dimension 3 is given in section 3.2 of
[math/0204356](https://arxiv.org/abs/math/0204356).

The high-level strategy is:

1. Generate all 4D IP weight systems.
2. Combine them into combined weight systems (CWS).
3. Convert the resulting CWS to a binary list of normal-form polytopes with
   `class.x`, closing the list under mirror symmetry.

A minimal local reproduction runs roughly as follows. First make sure you are
using executables compiled with `POLY_Dmax=4`, e.g. `cws-4d.x` and
`class-4d.x`.

```bash
# 1. Generate 4D IP weight systems (in dimension <=4 these are automatically
#    reflexive; the 'r' suffix in the output merely marks reflexivity).
./cws-4d.x -w4 > weights.txt
```

Note: this step enumerates a large number of weight systems and may take a
substantial amount of time and memory.

```bash
# 2. Combine the weight systems into CWS.
./cws-4d.x -c4 > cws4d.txt
```

For dimension 4, `cws.x -c4` creates all relevant combined weight systems by
default.

```bash
# 3. Convert the CWS list to a binary list of normal-form polytopes.
./class-4d.x -a -pi cws4d.txt -po class4d.bin
```

The `-a` option tells `class.x` to read ASCII input and write binary output in
the standard normal form; `-pi` and `-po` specify the input and output binary
files.

```bash
# 4. Check for missing mirrors.
./class-4d.x -M -pi class4d.bin > missing_mirrors.txt
```

The `-M` option prints any polytopes whose mirrors are not already present in
the binary list. After a *complete* generation this file should be **empty**;
that emptiness is the criterion used to verify that no polytope has been missed.
If it is not empty, the missing mirrors correspond to weight systems that were
not in the original input. Add those weight systems to the input list, rerun
`class.x`, and check again. Repeat until `missing_mirrors.txt` is empty. The
total number of distinct normal forms in the final closed list is the
Kreuzer–Skarke result: **473,800,776**.

For production use the recommended workflow is to use a binary database rather
than a single binary file:

```bash
./class-4d.x -di class4d_db -do class4d_db_complete
```

A database is a collection of files with a common prefix (`class4d_db`) rather
than one monolithic file; this was the format used for the original KS
classification and is what `class.x` optimizes for on very large runs.

### Important caveats

* **Compile with `POLY_Dmax=4`**: the classification code in `class.x` is
  intended to be run with `POLY_Dmax=4`. Use `class-4d.x`, which is compiled
  that way by default.
* **Time and memory**: the full classification can take days or weeks and uses
  many gigabytes of RAM. Start with a smaller subset (e.g. a single Hodge
  number or a constrained weight search) before attempting the full run.
* **Recoverability**: `class.x` writes intermediate results to
  `<output>.aux`. If a run is interrupted, resume with `-r` using the same
  output prefix: `class-4d.x -r -po class4d.bin`.
* **Sublattice polytopes**: by default `class.x` also finds reflexive
  subpolytopes on sublattices. Restrict to the original lattice with `-o` or
  `-o0` if only original-lattice polytopes are required.
* **Output format**: binary output from `class.x` is not human-readable. Use
  `class-4d.x -b -pi class4d.bin > class4d.txt` to dump the list of normal
  forms or weight systems to ASCII.

For further details on `class.x` and `cws.x`, see the
[PALP online documentation](http://palp.itp.tuwien.ac.at/wiki/index.php/PALP_online_documentation).

## License

PALP was originally written by Maximilian Kreuzer and Harald Skarke, and is
distributed under the GPLv3 license.
