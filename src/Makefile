CC=g++-5.1.0
MATPATH = /Applications/MATLAB_R2015a.app/extern/include
FLAGS = -I$(MATPATH) 
VERSION=-std=c++0x
OPTIONS= -g -O -c -fopenmp
LDFLAGS= -Wl,-rpath,/Applications/MATLAB_R2015a.app/bin/maci64
all: Convection_1D_ULTIMATE_QUICKEST6.cpp Convection_2D.cpp Diffusion_1D.cpp Diffusion_2D.cpp Diffusion_ADI1.cpp Diffusion_ADI2.cpp Diffusion_ADI3.cpp Interpolation.cpp Logger.cpp Matrix.cpp MatrixSolver.cpp MonotCubicInterpolator.cpp Parameters.cpp ReadInitialData.cpp UpdatableMatrix.cpp VERB4D_Solver.cpp		
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) Convection_1D_ULTIMATE_QUICKEST6.cpp
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) Convection_2D.cpp
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) Diffusion_1D.cpp
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) Diffusion_2D.cpp
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) Diffusion_ADI1.cpp
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) Diffusion_ADI2.cpp
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) Diffusion_ADI3.cpp
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) Interpolation.cpp
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) Logger.cpp
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) Matrix.cpp
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) MatrixSolver.cpp -llapack -lblas
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) MonotCubicInterpolator.cpp
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) Parameters.cpp
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) ReadInitialData.cpp
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) UpdatableMatrix.cpp
	  $(CC) $(VERSION) $(OPTIONS) $(FLAGS) VERB4D_Solver.cpp
	  $(CC) $(VERSION) -fopenmp -o VERB4D_Solver VERB4D_Solver.o Convection_1D_ULTIMATE_QUICKEST6.o Convection_2D.o Diffusion_1D.o Diffusion_2D.o Diffusion_ADI1.o Diffusion_ADI2.o Diffusion_ADI3.o Interpolation.o Logger.o Matrix.o MatrixSolver.o MonotCubicInterpolator.o Parameters.o ReadInitialData.o UpdatableMatrix.o -llapack -lblas -L"/Applications/MATLAB_R2015a.app/bin/maci64" -lmx -lmex -lmat -lc++ -leng $(LDFLAGS)
clean:
	$(RM) Convection_1D_ULTIMATE_QUICKEST6.o
	$(RM) Convection_2D.o
	$(RM) Diffusion_1D.o
	$(RM) Diffusion_2D.o
	$(RM) Diffusion_ADI1.o
	$(RM) Diffusion_ADI2.o
	$(RM) Diffusion_ADI3.o
	$(RM) Interpolation.o
	$(RM) Logger.o
	$(RM) Matrix.o
	$(RM) MatrixSolver.o
	$(RM) MonotCubicInterpolator.o
	$(RM) Parameters.o
	$(RM) ReadInitialData.o
	$(RM) UpdatableMatrix.o
	$(RM) VERB4D_Solver.o
	$(RM) VERB4D_Solver

