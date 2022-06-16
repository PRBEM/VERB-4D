# VERB4D
Source code for the VERB-4D solver. Need to be used together with Matlab or Python setup scripts that
create all necessary input. The main file is in `./src/VERB4D_Solver.cpp`

## Build instructions
Download cmake 2.8 and above. Note that any new source files should be added to
the VERB4D_SOURCES line in CMakeLists.txt, and the build files must be
regenerated with the below commands.

### OSX / Linux
The following dependencies are required and will be automatically searched via
cmake. Use homebrew or the package manager of your choice to install them:

* MatLab (optional)
* OpenBLAS (added as a submodule)

On OSX, you will need to install gcc or g++ in order to support compilation
with OpenMP.  The version of gcc on OSX 10.10 is actually linked to clang. You
will need to install a version of gcc from homebrew or ports. To test that it
works, run the command `gcc-5 --version` and verify that it is indeed gcc and
not clang.

Clone the repository recursively to also get the OpenBLAS library

#### via ssh 
```
git clone --recurse-submodules git@rbm9.gfz-potsdam.de:verb/verb4d_solver.git
```

#### via https
```
git clone --recurse-submodules https://rbm9.gfz-potsdam.de/verb/verb4d_solver.git
```

#### Compile and install the OpenBlas library
Enter the verb4d_solver directory and build the OpenBLAS library and install it to `~/.local/lib` by calling
```
cd verb4d_solver
./install_openblas.sh
```

#### Compile and install the VERB-4D solver via script
Compile the VERB4D solver and install it to `~/.local/bin` by calling
```
source install_verb.sh
```
You can now call the executable `VERB4D_Solver` from a simulation folder that contains a `parameters.ini` file.

#### Compile manually
Create a build folder inside `verb4d_solver`, run CMake with a Release build type and call make to compile
```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release`
make
```

### Windows
The CMake instructions will also work for windows. You will want to create a
new folder named build inside of the root directory. In the gui, set the
generator name to "Visual Studio 12" for a Visual Studio 2013 *.sln project.

On the command line, you can run:

```
$ mkdir build
$ cd build
$ cmake .. -G "Visual Studio 14"
```

Replace `Visual Studio 14` with the version of your choice. For example,
`Visual Studio 14` corresponds to VS 2015.

Note that the solution includes three projects.

ZERO_CHECK will rerun cmake. You can/should execute this after changing
something on your CMake files.

ALL_BUILD is simply a target which builds all and everything project in the
active solution, one can compare it to "make all".

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
