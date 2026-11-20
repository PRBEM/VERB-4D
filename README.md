<!--
SPDX-FileCopyrightText: 2015 UCLA
SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences

SPDX-License-Identifier: BSD-3-Clause
-->

# VERB-4D

**Versatile Electron Radiation Belt in 4D (VERB-4D)** is a computational solver that models radiation belt, ring current, and plasmasphere dynamics in Earth's magnetosphere.

We recommend to use the [verb4dpy](https://git.gfz-potsdam.de/spsw/verb/verb4d/verb4dpy) package to create all necessary inputs for the solver; [Install with uv](#install-with-uv) installs both. For a first run see [First simulation with the `verb-4d` command](#first-simulation-with-the-verb-4d-command); tutorials and documentation are in verb4dpy.

For questions, please contact <yshprits@gfz.de> or <bhaas@gfz.de>.

## Get the source
```
git clone git@git.gfz-potsdam.de:spsw/verb/verb4d/verb-4d.git      # via ssh
git clone https://git.gfz-potsdam.de/spsw/verb/verb4d/verb-4d.git  # via https
cd verb-4d
```
The default builds need no submodule. Fetch only the ones your build options need:

| Build and options | Submodule | Fetch with |
| ----------------- | --------- | ---------- |
| [Install with uv](#install-with-uv), any options | none (uv fetches pybind11 itself) | – |
| [CMake](#build-with-cmake) with `-DBLAS_TYPE=LAPACK` (default) or `MKL` | none | – |
| CMake with `-DBLAS_TYPE=OpenBLAS`, OpenBLAS built by `install_openblas.sh` | `OpenBLAS` | `git submodule update --init OpenBLAS` |
| CMake with `-DPYTHON_BINDINGS=ON` | `extern/pybind11`, unless CMake finds an installed pybind11 (e.g. `-Dpybind11_DIR=$(python -m pybind11 --cmakedir)`) | `git submodule update --init extern/pybind11` |
| `USE_PPFV=ON` (uv or CMake) | none: Eigen 3.4 and xtensor 0.24.7 are taken from the system or downloaded by CMake while configuring (needs internet) | – |

`git clone --recursive` fetches all submodules at once.

## Requirements
* A C++17 compiler, e.g. gcc 10.2 or above. If CMake picks the wrong one, set `CXX`, e.g. `export CXX=/u/local/compilers/gcc/10.2.0/bin/g++`.
* LAPACK and BLAS (e.g. `apt install liblapack-dev` or `brew install lapack`), or OpenBLAS (see the table above), or Intel's MKL.
* OpenMP (optional; without it the solver runs serially).
* For [Install with uv](#install-with-uv): [uv](https://docs.astral.sh/uv/). uv fetches CMake if the system has none or one older than 3.15.
* For [Build with CMake](#build-with-cmake): CMake 3.14 or above.
* A Fortran compiler only to build OpenBLAS; without one, OpenBLAS is compiled from C.

Windows: *coming soon*

## Install with uv
Into a new environment, from the clone:
```
uv venv
uv pip install .
```
Without a clone:
```
uv pip install "verb4d-solver @ git+https://git.gfz-potsdam.de/spsw/verb/verb4d/verb-4d.git"
```

This builds the solver with CMake (through [scikit-build-core](https://scikit-build-core.readthedocs.io)) and installs
it together with the Python bindings and [verb4dpy](https://git.gfz-potsdam.de/spsw/verb/verb4d/verb4dpy) (the
version this solver release pins).

### Build options
The package is built with these [compilation flags](#list-of-compilation-flags):

| Flag | Package default |
| ---- | --------------- |
| `CMAKE_BUILD_TYPE` | `Release` |
| `PYTHON_BINDINGS` | `ON` |
| `BLAS_TYPE` | `LAPACK` (the system LAPACK; the only one tested for this install) |
| `FAST_CONVECTION` | `ON` (unlike the plain CMake build, where it is `OFF`) |
| `USE_PPFV`, `DATA_ASSIMILATION`, `LU_CACHING`, `SAVE_PSD_LOST_CONV` | `OFF` |

Change one with `-C` (pip takes the same option):
```
uv pip install . -Ccmake.define.FAST_CONVECTION=OFF -Ccmake.define.USE_PPFV=ON
```
In Python:
```python
import verb4d_solver

verb4d_solver.executable_path()  # the installed VERB4D_Solver
verb4d_solver.build_options  # {"FAST_CONVECTION": False, "USE_PPFV": True, ...}
verb4d_solver.Diffusion_2D(...)  # the Python bindings
```

### First simulation with the `verb-4d` command
The install also brings the `verb-4d` command of verb4dpy. Download the model data once, then run a three-hour
plasmasphere simulation and make a movie of it:
```
source .venv/bin/activate
python -m verb4dpy.data download   # model data, once (~/.local/share/verb4dpy on Linux)
verb-4d run-default plasmasphere 2013-03-17 2013-03-17T03:00
verb-4d plot plasmasphere-movie plasmasphere_20130317T0000_20130317T0300
verb-4d help                       # list all commands; `verb-4d help COMMAND` for one
```
The simulation is created in `./plasmasphere_20130317T0000_20130317T0300`, the movie in its `Analysis/` folder.

More in verb4dpy: the [README](https://git.gfz-potsdam.de/spsw/verb/verb4d/verb4dpy) (config files, the interactive
shell, comparisons), the [tutorial notebooks](https://git.gfz-potsdam.de/spsw/verb/verb4d/verb4dpy/-/tree/main/docs/tutorials)
(installation and first run, configuration, plotting, comparisons, satellite data, the Python API) and its
[documentation](https://spsw.git-pages.gfz-potsdam.de/verb/verb4d/verb4dpy/).

## Build with CMake
Builds only the solver executable (and, with `-DPYTHON_BINDINGS=ON`, the bindings library) in a build folder, e.g.
for development or on a cluster. New source files must be added to `VERB4D_SOURCES` in `CMakeLists.txt`; then run
CMake again.

### List of compilation flags

| Command line argument | Description | Possible values |
| --------------------- | ----------- | --------------- |
| -DCMAKE_BUILD_TYPE | Standard cmake option enabling debugging symbols or optimization flags | Release, Debug |
| -DBLAS_TYPE | Library type of LAPACK and BLAS routines | LAPACK (default), OpenBLAS, MKL |
| -DOPENBLAS_LIBPATH | Path to the compiled OpenBLAS library (default `~/.local/lib/libopenblas.so`) |  |
| -DDATA_ASSIMILATION | Bool whether the data assimilation routines will be compiled | ON, OFF |
| -DDATA_ASSIMILATION_DEBUG | Bool whether the data assimilation debug output will be saved (!huge file sizes!) | ON, OFF |
| -DFAST_CONVECTION | Bool to decide convection step size in Convection_2D | ON, OFF |
| -DPYTHON_BINDINGS | Bool whether the python-bindings-library will be compiled | ON, OFF |
| -DSAVE_PSD_LOST_CONV | Bool whether the lost PSD during Convection_2d will be saved | ON, OFF |
| -DLU_CACHING | Bool whether LU caching will be used for local diffusion. It's faster, but needs more memory. | ON, OFF |
| -DUSE_PPFV | Bool whether the positive-preserving finite volume method should be compiled | ON, OFF |

All flags are `OFF` by default. The [table above](#get-the-source) lists the submodules some of them need.

### Linux / OSX
#### CMake
Configure a release build with the BLAS library of your choice, then compile:
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBLAS_TYPE=OpenBLAS
cmake --build build -j
```
The executable `VERB4D_Solver` is in `./build`.

#### Shell script
Compile the VERB-4D solver with OpenBLAS into `./build` by calling
```
./install_verb.sh
```
It expects the OpenBLAS library in `~/.local/lib` (see below).

### Windows
*coming soon*

### Verify the installation
Once built, check that the executable and its BLAS/LAPACK and OpenMP dependencies were linked correctly, without needing any `parameters.ini` or input files:
```
./build/VERB4D_Solver --check
```
For the package installed with uv or pip:
```
"$(python -c 'import verb4d_solver; print(verb4d_solver.executable_path())')" --check
```
`--version` prints build info (BLAS backend, enabled compile options, OpenMP thread count, commit hash if available) and `--help` prints usage. Running the solver itself additionally requires a `parameters.ini` file and associated input files in the working directory.

### Compile and install the OpenBLAS library (optional)
Needs the `OpenBLAS` submodule (`git submodule update --init OpenBLAS`).

Linux / OSX: build the shared OpenBLAS library and install it to `~/.local/lib` by calling
```
./install_openblas.sh
```
Windows: *coming soon*


## Development History
- **Originally developed**: 2015 at UCLA (University of California, Los Angeles)
- **Principal Investigator**: Yuri Shprits
- **Key developers**: Dmitri Subbotin, Alexander Drozdov, Adam Kellerman
- **License**: BSD-3-Clause

## Contributors (Alphabetically Sorted)

- Josh Adler
- Hayley Allison
- Nikita Aseev
- Alexander Drozdov
- Bernhard Haas
- Julia Himmelsbach
- Alec Jen
- Adam Kellerman
- Kyung-Chan Kim
- Joseph Lin
- Anthony Miyaguchi
- Derek Nguyen
- Maxim Rakhuba
- Anthony Saikin
- Yuri Shprits
- Dmitri Subbotin
- Derek Nguyen
- Dedong Wang
- Michael Wutzig
- Irina Zhelavskaya

## Acknowledgments

We acknowledge Xin Tao and all developers of the [Sayram code](https://github.com/xtaohub/Sayram-2D/) for providing the C++ code for the positive-preserving finite volume scheme.

We also acknowledge indirect contributors from the scientific community for their helpful discussions and contributions to the development of VERB-4D.
