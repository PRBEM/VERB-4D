
#ifndef CONVECTION_1D_ULTIMATE_QUICKEST6_H_
#define CONVECTION_1D_ULTIMATE_QUICKEST6_H_

#include "Matrix.h"
#include "MatrixSolver.h"


bool Convection_1D_ULTIMATE_QUICKEST6 ( Matrix1D < double > &PSD,
		Matrix1D < double > x,
		int x_size,
		double x_LBC, double x_UBC,
		string x_LBC_type, string x_UBC_type,
		Matrix1D < double > Dx,
		Matrix1D < double > Sources, Matrix1D < double > Losses,
		double dt);


#endif
