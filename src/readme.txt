ReadMe - COMPILING ON MAC



Note: All tabbed items should be executed on the command line

download Xcode with developer tools - located in app store (free)
(developer tools might need to be downloaded separately)

in terminal write:
	xcode-select --install

download gcc-5.1.0 as explained at https://solarianprogrammer.com/2015/05/01/compiling-gcc-5-mac-os-x/

Get the verb4d folder from gitlab

In the terminal go to the src directory which contains VERB4D_Solver.cpp from verb4d
	cd ~/…(your path if any).../verb4d/VERB4D_Solver/src

add paths for both gcc and omp to work
	export PATH=/usr/gcc-5.1.0/bin:$PATH
Assuming you kept the gcc-5.1.0 file in your downloads like in the above tutorial
	export PATH=~/Downloads/gcc-5.1.0/libgomp:$PATH

Might need to change <malloc.h> to <stdlib.h> in the file MatrixSolver.cpp

Create the executable
	make

You now have an executable name VERB4D_Solver

To test that it works:

Copy the executable “VERB4D_Solver” into the 1.0 1D_Alpha_Example folder
	cp VERB4D_Solver ../../Examples/1.0\ 1D_Alpha_Example/VERB4D_Solver

Open Conv_Dif.m from that file in matlab

Scroll down to the run portion and change system('VERB4D_Solver.exe’) to system('./VERB4D_Solver')

press command enter to run

go to plot - press command enter to run
