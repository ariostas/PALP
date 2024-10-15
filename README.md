# PALP: a Package for Analyzing Lattice Polytopes 

This is a fork of the [PALP GitLab repository](https://gitlab.com/stringstuwien/PALP). The `upstream` branch tracks the `main` branch of the original repository, and other branches of this repository contain changes to allow the implementation of Python bindings. More information about the PALP package can be found at the [PALP website](http://hep.itp.tuwien.ac.at/~kreuzer/CY/CYpalp.html).

## Building and installing

PALP can now be built using CMAKE. First, generate the build configuration with
```bash
cmake -S . -B build -D CMAKE_BUILD_TYPE=Release
```

Optionally, the install location and maximum polytope dimension can be specified with
```bash
-D CMAKE_INSTALL_PREFIX=[your custom path]
-D POLY_Dmax=[your custom dimension, defaults to 6]
```

Then, build the package with
```bash
cmake --build build
```

Optionally, the binaries can then be installed with
```bash
cmake --build build --target install
```

## Usage

Please consult the [PALP online documentation](http://palp.itp.tuwien.ac.at/wiki/index.php/PALP_online_documentation) for detailed instructions.
