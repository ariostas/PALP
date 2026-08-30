# PALP: a Package for Analyzing Lattice Polytopes 

This is a fork of the [PALP GitLab repository](https://gitlab.com/stringstuwien/PALP). The `upstream` branch tracks the `main` branch of the original repository, and other branches of this repository contain changes to allow the implementation of Python bindings. More information about the PALP package can be found at the [PALP website](http://hep.itp.tuwien.ac.at/~kreuzer/CY/CYpalp.html).

## Building and installing

Building PALP requires a C compiler and CMake 3.14 or newer. Configure the build tree with

```bash
cmake -S . -B build
```

This defaults to a `Release` build (`-O3`), matching what the old `GNUmakefile` produced; pass `-D CMAKE_BUILD_TYPE=Debug` if you want an unoptimized build instead.

The following options can be passed at configure time:

| Option | Default | Description |
| --- | --- | --- |
| `POLY_Dmax` | `6` | Maximum polytope dimension the executables are compiled for |
| `CMAKE_INSTALL_PREFIX` | `/usr/local` on Unix | Where the `install` step places the executables |
| `PALP_BUILD_TESTS` | `OFF` | Also build the test suite (see [Testing](#testing)) |

`POLY_Dmax` replaces the old procedure of editing `Global.h` by hand. Note that the programs still work when it is set higher than needed, but they may be considerably slower, so it is worth setting it to the dimension you actually analyse. For example:

```bash
cmake -S . -B build -D POLY_Dmax=4 -D CMAKE_INSTALL_PREFIX=~/.local
```

Then build the package with

```bash
cmake --build build
```

This produces `poly.x`, `class.x`, `cws.x`, `nef.x` and `mori.x` in `build/`. Optionally, they can then be installed into `<prefix>/bin` with

```bash
cmake --install build
```

PALP can also be consumed by another CMake project via `add_subdirectory`, in which case only the `objects`, `class_objects`, `nef_objects` and `mori_objects` libraries are defined; the executables, the install rules and the tests are skipped.

## Testing

The shell test suite in `tests/` is wired up to CTest. Each executable has to be compiled with the `POLY_Dmax` it is tested at, so the suite builds its own set of programs (`poly-4d.x`, `poly-5d.x`, ...) and is therefore opt-in:

```bash
cmake -S . -B build -D PALP_BUILD_TESTS=ON
cmake --build build --target check
```

The `check` target corresponds to `make check`, and `checklong` to `make checklong`; the latter additionally runs the handful of checks that take a very long time. The dimensions to build and test can be changed with `-D PALP_TEST_DIMENSIONS="4;6"` (defaults to `4;5;6;11`).

`ctest` can also be run directly from the build directory, which makes it easy to select a subset of the suite:

```bash
cd build
ctest --label-exclude long            # same as the "check" target
ctest -L dim6                         # only the 6d executables
ctest -R nef                          # only the nef tests
```

The original `GNUmakefile` is still present and works as before.

## Usage

Please consult the [PALP online documentation](http://palp.itp.tuwien.ac.at/wiki/index.php/PALP_online_documentation) for detailed instructions.

## License

PALP was originally written by Maximilian Kreuzer and Harald Skarke, and is distributed under the GPLv3 license.
