<!--
SPDX-FileCopyrightText: 2015 UCLA
SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences

SPDX-License-Identifier: BSD-3-Clause
-->

# VERB-4D solver
Source code for the VERB-4D solver. We recommend to use the [verb4d_py](https://git.gfz-potsdam.de/spsw/verb/verb4d/verb4dpy) package to create all necessary inputs for the solver. The main file is in `./src/VERB4D_Solver.cpp`.

## Requirements
Download cmake 3.14 or above. Note that any new source files should be added to
the VERB4D_SOURCES line in CMakeLists.txt, and the build files must be
regenerated with the below commands.

The following dependencies are required and will be automatically searched via
cmake. Use homebrew or the package manager of your choice to install them:

* OpenMP (optional)
* Lapack, OpenBLAS (added as a submodule) or Intel's MKL

A C++ compiler is necessary, a Fortran compiler (only for OpenBLAS) recommended.
If no Fortran compiler is installed, OpenBLAS will compile from C code instead.

#### Linux
The gcc and gfortran compilers can easily be installed if not already available.
If you have a compilation problem, install the latest gcc (minimum 10.2.0).
If CMAKE cannot find the correct g++ version, you can set the environmental variable CXX to the correct binary, e.g.

`export CXX=/u/local/compilers/gcc/10.2.0/bin/g++`

#### Windows
*coming soon*

## Build instructions
Clone the verb4d_solver repository

#### via ssh 
```
git clone git@git.gfz-potsdam.de:spsw/verb/verb4d/verb4d_solver.git
```

#### via https
```
git clone https://git.gfz-potsdam.de/spsw/verb/verb4d/verb4d_solver.git
```
Add the `--recursive` flag to also clone the OpenBLAS and Pybind11 submodule. Enter the verb4d_solver directory with `cd verb4d_solver`. If you cloned verb4d_solver non-recursively, you can update the submodules by executing `git submodule update --init`.

### Compile the VERB-4D solver

#### List of compilation flags

| Command line argument | Description | Possible values
| --------------------- | ----------- | -------------- |
| -DCMAKE_BUILD_TYPE | Standard cmake option enabling debugging symbols or optimization flags | Release, Debug |
| -DBLAS_TYPE | Library type of LAPACK and BLAS routines | LAPACK, OpenBLAS, MKL |
| -DOPENBLAS_LIBPATH | Path towards compiled openblas.lib |  |
| -DDATA_ASSIMILATION | Bool whether the data assimilation routines will be compiled | True, False |
| -DDATA_ASSIMILATION_DEBUG | Bool whether the data assimilation debug output will be saved (!huge file sizes!) | True, False |
| -DFAST_CONVECTION | Bool to decide convection step size in Convection_2D | True, False |
| -DPYTHON_BINDINGS | Bool whether the python-bindings-library will be compiled | True, False |
| -DSAVE_PSD_LOST_CONV | Bool whether the lost PSD during Convection_2d will be saved | True, False |
| -DLU_CACHING | Bool whether LU caching will be used for local diffusion. It's faster, but needs more memory. | True, False |
| -DUSE_PPFV | Bool wheter the positive-preserving finite volume method should be compiled. | True, False |

#### Linux / OSX
##### CMake
Manually compile the VERB-4D solver by creating a build folder
```
mkdir build && cd build
```
and call CMake for a release built and the BLAS library of your choice, either Lapack (default), OpenBLAS or MKL.
```
cmake .. -DCMAKE_BUILD_TYPE=Release -DBLAS_TYPE=OpenBLAS
make -j
```
For the cmake call, when choosing OpenBLAS, an additional parameter `-DOPENBLAS_LIBPATH=<path to libopenblas.so>` can be specified (default `~/.local/lib/libopenblas.so`)
##### Shell script
Automatically compile the VERB-4D solver (using OpenBLAS) by calling
```
./install_verb.sh
```
The executable `VERB4D_Solver` can be found in `./build`.

#### Windows
##### CMake (PowerShell)
*coming soon*

### Compile and install the OpenBLAS library (optional)
#### Linux / OSX
Build the shared OpenBLAS library and install it to `~/.local/lib` by calling
```
./install_openblas.sh
```
#### Windows
*coming soon*

## Contributors and Acknowledgements

We thank all the contributors of the VERB-4D Solver.

We acknowledge Xin Tao and all developers of the [Sayram code](https://github.com/xtaohub/Sayram-2D/) for providing the C++ code for the positive-preserving finite volume scheme.