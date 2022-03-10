/**
 * \file Convection_2D.h
 *
 * Using the Convection_1D_Ultimate_QUICKEST6() for calculating P and R.
 * All equations and formulas for these calculations can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf .
 * Mostly coming from 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme.
 * In the source code mathematical equations are numbered corresponding to the numbering of equations found in the paper
 *
 *
 * \brief Calculates the convection in 2D given a 2D matrix of Phase Space Densities, P, R, boundary conditions, diffusion, sources and losses
 *
 *  Created on: Aug 22, 2011
 *      Author: subbotin
 */

#ifndef CONVECTION_2D_ULTIMATE_QUICKEST_H_
#define CONVECTION_2D_ULTIMATE_QUICKEST_H_

#include "Matrix.h"
#include "MatrixSolver.h"


#include "Convection_1D_ULTIMATE_QUICKEST6.h"


/**
* Function that creates a 2D Convection matrix and returns a bool upon completion
*/
bool Convection_2D( 
	Matrix2D < double > &PSD_PR,
	Matrix2D < double > P, Matrix2D < double > R,
	int P_size, int R_size,
	Matrix1D < double > P_LBC, Matrix1D < double > P_UBC,
	Matrix1D < double > R_LBC, Matrix1D < double > R_UBC,
	std::string P_LBC_type, std::string P_UBC_type,
	std::string R_LBC_type, std::string R_UBC_type,
	Matrix2D < double > DP, Matrix2D < double > DR,
	Matrix2D < double > Sources, Matrix2D < double > Losses,
	double dt, double min_PSD, double min_V
);

#endif
