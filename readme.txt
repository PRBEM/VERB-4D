ReadMe - COMPILING ON MAC

Note: All tabbed items should be executed on the command line

Xcode does not come with openmp support so a different compiler is necessary. However you will still need Xcode for some tools in order to build the new compiler

download Xcode with developer tools - located in app store (free)
(developer tools might need to be downloaded separately)

in terminal write:
	xcode-select --install

download gcc-5.1.0 as explained at https://solarianprogrammer.com/2015/05/01/compiling-gcc-5-mac-os-x/
NOTE: some of the links seem to not work anymore. Just copy and paste the name of the desired file to download into google and you should be able to find it.

Get the verb4d folder from gitlab

In the terminal go to the src directory which contains VERB4D_Solver.cpp from verb4d
	cd ~/…(your path if any).../verb4d/VERB4D_Solver/src

add paths for both gcc and omp to work
	export PATH=/usr/gcc-5.1.0/bin:$PATH
Assuming you kept the gcc-5.1.0 file in your downloads like in the above tutorial
	export PATH=~/Downloads/gcc-5.1.0/libgomp:$PATH

Might need to change <malloc.h> to <stdlib.h> in the file Matrix.h

Create the executable
	make

You now have an executable name VERB4D_Solver

To test that it works:

Copy the executable “VERB4D_Solver” into the 1.0 1D_Alpha_Example folder
	cp VERB4D_Solver ../../Examples/1.0\ 1D_Alpha_Example/VERB4D_Solver

Open Conv_Dif.m from that file in matlab

Scroll down to the run portion and change system('VERB4D_Solver.exe’) to system('./VERB4D_Solver')

press run (play button) at top of page



COMPILING not on mac

Must have OPENMP2.0 or above - comes with visual studio etc.
Must have c++ compiler such as gcc or g++ etc.
Change the CC variable in Makefile to whatever compiler you are using


MATLAB CAPABILITY

If you have matlab installed on your machine then you can opt to read/write files into .mat format which is faster to use and more precise than text files.

Make sure MATLAB_CAPABLE is set to true in Matrix.h

In all Conv_Dif.m files there is a variable named use_matlab. Set this to true as well.

Also to compile the solver you must link to the matlab/c++ shared libraries:
for mac they can be found at /Applications/MATLAB_R2015a.app/bin/maci64 assuming you downloaded matlab into applications folder. Look at the Makefile for an idea of how to link the proper Matlab libraries for compiling.

For correct path type matlabroot into matlab, and then add on the necessary extensions such as extern/include and bin/os

If you do not have matlab capabilities use the MakefileNoMatlab to compile after setting MATLAB_CAPABLE to false in Matrix.h.

NOTES:

For the standard testing of an example include_boundary should be true and inversion method should be lapack. The ADI inversion methods are not working properly as of Sept 2015


