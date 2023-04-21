# VERB-4D solver
Source code for the VERB-4D solver. Need to be used together with Matlab or Python setup scripts that
create all necessary input. The main file is in `./src/VERB4D_Solver.cpp`

## Requirements
Download cmake 3.7 or above. Note that any new source files should be added to
the VERB4D_SOURCES line in CMakeLists.txt, and the build files must be
regenerated with the below commands.

The following dependencies are required and will be automatically searched via
cmake. Use homebrew or the package manager of your choice to install them:

* OpenMP
* Lapack, OpenBLAS (added as a submodule) or Intel's MKL
* MatLab (optional)

A C++ compiler is necessary, a Fortran compiler (only for OpenBLAS) recomended.
If no Fortran compiler is installed, OpenBLAS will compile from C code instead.

#### Linux
The gcc and gfortran compilers can easily be installed if not already available.
If you have a compilation problem, install the latest gcc (minimum 10.2.0).
If CMAKE cannot find the correct g++ version, you can setup environmental variable CXX:

export CXX=/u/local/compilers/gcc/10.2.0/bin/c++

#### Windows
##### Visual Studio
The MSVC C++ compiler can be installed using [Visual Studio](https://visualstudio.microsoft.com/downloads/).
When choosing OpenBLAS, the only option in Visual Studio is to compile OpenBLAS from C code,
which might be less optimized than building from the Fortran source code.

##### MinGW and Ninja
The recommended way of building OpenBLAS on Windows is by using the [Ninja](https://ninja-build.org/)
build system instead of Visual Studio. Additionally [MSYS2](https://www.msys2.org/) is necessary
to install MinGW's gcc and gfortran compiler. To install both compilers, open 
an MSYS2 MinGW terminal and run
```
pacman -S mingw-w64-x86_64-gcc
pacman -S mingw-w64-x86_64-gcc-fortran
pacman -Syu
```
Make sure that the location of `ninja.exe` and the MinGW binary folder, e.g. `C:\msys64\mingw64\bin\`, are on your PATH.
#### OSX
You will need to install gcc or g++ in order to support compilation
with OpenMP. The version of gcc on OSX 10.10 is actually linked to clang. You
will need to install a version of gcc from homebrew or ports. To test that it
works, run the command `gcc-5 --version` and verify that it is indeed gcc and
not clang.

## Build instructions
Clone the verb4d_solver repository

#### via ssh 
```
git clone git@rbm9.gfz-potsdam.de:verb/verb4d_solver.git
```

#### via https
```
git clone https://rbm9.gfz-potsdam.de/verb/verb4d_solver.git
```
Enter the verb4d_solver directory, checkout the openblas branch and update the
OpenBLAS submodule (this option is currently possible for verb4d_solver openblas branch).
```
cd verb4d_solver
git checkout openblas
git submodule update --init
```

### Compile and install the OpenBLAS library (optional)
#### Linux / OSX
Build the shared OpenBLAS library and install it to `~/.local/lib` by calling
```
./install_openblas.sh
```
### Windows
For the Visual Studio build, open a PowerShell terminal and execute
```
install_openblas_vs.ps1
```
For the Ninja build, make sure that both gcc and gfortran are installed (e.g. `gfortran --version`).
Call the Ninja install script
```
install_openblas_ninja.ps1
```
Both build types compile OpenBLAS as a dynamic libray and copy the library files to `<drive>:\Users\<user>\AppData\Local\OpenBLAS\`.
This location is automatically appended to the user's PATH environment variable, so VERB-4D can use it at runtime.

### Compile the VERB-4D solver
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
For the cmake call, when choosing OpenBLAS, an additional paramter `-DOPENBLAS_LIBPATH=<path to libopenblas.so>` can be specified (default `~/.local/lib/libopenblas.so`)
##### Shell script
Automatically compile the VERB-4D solver (using OpenBLAS) by calling
```
./install_verb.sh
```
The executable `VERB4D_Solver` can be found in `./build`.

#### Windows
##### CMake (PowerShell)
Manually compile the VERB-4D solver by creating a build folder
```
mkdir build && cd build
```
and call CMake for a release built and the BLAS library of your choice, either Lapack (default), OpenBLAS or MKL.
```
cmake .. -DBLAS_TYPE=OpenBLAS
cmake --build . -j --config Release
```
For the first cmake call, when choosing OpenBLAS, an additional paramter `-DOPENBLAS_LIBPATH=<path to openblas.lib>` can be specified <br>
 (default `<drive>:\Users\<user>\AppData\Local\OpenBLAS\openblas.lib`)
##### PowerShell script
Automatically compile the VERB-4D solver (using OpenBLAS) by calling
```
install_verb.ps1
```
The executable `VERB4D_Solver.exe` can be found in `.\build\Release`.

### Use the VERB-4D solver

Copy the executable to a simulation folder that contains a `parameters.ini` file and run it to start a simulation.


## Old Notes (Possibly outdated)
ReadMe - COMPILING ON MAC

Note: All tabbed items should be executed on the command line

Xcode does not come with openmp support so a different compiler is necessary.
However you will still need Xcode for some tools in order to build the new
compiler

download Xcode with developer tools - located in app store (free) (developer
tools might need to be downloaded separately)

in terminal write: xcode-select --install

download gcc-5.1.0 as explained at
https://solarianprogrammer.com/2015/05/01/compiling-gcc-5-mac-os-x/ NOTE: some
of the links seem to not work anymore. Just copy and paste the name of the
desired file to download into google and you should be able to find it.

Get the verb4d folder from gitlab

In the terminal go to the src directory which contains VERB4D_Solver.cpp from
verb4d cd ~/…(your path if any).../verb4d/VERB4D_Solver/src

add paths for both gcc and omp to work export PATH=/usr/gcc-5.1.0/bin:$PATH
Assuming you kept the gcc-5.1.0 file in your downloads like in the above
tutorial export PATH=~/Downloads/gcc-5.1.0/libgomp:$PATH

Might need to change <malloc.h> to <stdlib.h> in the file Matrix.h

Create the executable make

You now have an executable name VERB4D_Solver

To test that it works:

Copy the executable “VERB4D_Solver” into the 1.0 1D_Alpha_Example folder cp
VERB4D_Solver ../../Examples/1.0\ 1D_Alpha_Example/VERB4D_Solver

Open Conv_Dif.m from that file in matlab

Scroll down to the run portion and change system('VERB4D_Solver.exe’) to
system('./VERB4D_Solver')

press run (play button) at top of page

COMPILING not on mac

Must have OPENMP2.0 or above - comes with visual studio etc.  Must have c++
compiler such as gcc or g++ etc.  Change the CC variable in Makefile to
whatever compiler you are using

MATLAB CAPABILITY

If you have matlab installed on your machine then you can opt to read/write
files into .mat format which is faster to use and more precise than text files.

Make sure MATLAB_CAPABLE is set to true in Matrix.h

In all Conv_Dif.m files there is a variable named use_matlab. Set this to true
as well.

Also to compile the solver you must link to the matlab/c++ shared libraries:
for mac they can be found at /Applications/MATLAB_R2015a.app/bin/maci64
assuming you downloaded matlab into applications folder. Look at the Makefile
for an idea of how to link the proper Matlab libraries for compiling.

For correct path type matlabroot into matlab, and then add on the necessary
extensions such as extern/include and bin/os

If you do not have matlab capabilities use the MakefileNoMatlab to compile
after setting MATLAB_CAPABLE to false in Matrix.h.

NOTES:

For the standard testing of an example include_boundary should be true and
inversion method should be lapack. The ADI inversion methods are not working
properly as of Sept 2015
