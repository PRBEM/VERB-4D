/*
 *
 *  Created on: Aug 22, 2011
 *      Author: subbotin
 */

#ifndef CONVECTION_2D_ULTIMATE_QUICKEST_H_
#define CONVECTION_2D_ULTIMATE_QUICKEST_H_

#include "Matrix.h"
#include "MatrixSolver.h"


#include "Convection_1D_ULTIMATE_QUICKEST6.h"

bool Convection_2D( Matrix2D < double > &PSD_PR,
			Matrix2D < double > P, Matrix2D < double > R,
			int P_size, int R_size,
			Matrix1D < double > P_LBC, Matrix1D < double > P_UBC,
			Matrix1D < double > R_LBC, Matrix1D < double > R_UBC,
			string P_LBC_type, string P_UBC_type,
			string R_LBC_type, string R_UBC_type,
			Matrix2D < double > DP, Matrix2D < double > DR,
			Matrix2D < double > Sources, Matrix2D < double > Losses,
			double dt, double min_PSD, double min_V);

#endif
