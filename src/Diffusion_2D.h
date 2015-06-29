/*
 * Diffusion.h
 *
 *  Created on: May 28, 2011
 *      Author: dimath
 */

#ifndef DIFFUSION_2D_H_
#define DIFFUSION_2D_H_

#include "Matrix.h"
#include "MatrixSolver.h"


bool Diffusion_2D(
		  Matrix2D<double> &psd,
		  Matrix2D<double> Q1, Matrix2D<double> Q2,
		  int Q1_size, int Q2_size,
		  Matrix1D<double> Q1_LBC, Matrix1D<double> Q1_UBC,
		  Matrix1D<double> Q2_LBC, Matrix1D<double> Q2_UBC,
		  string Q1_LBC_type, string Q1_UBC_type,
		  string Q2_LBC_type, string Q2_UBC_type,
		  Matrix2D<double> DQ1Q1, Matrix2D<double> DQ2Q2, Matrix2D<double> DQ1Q2, Matrix2D<double> DQ2Q1,
		  Matrix2D<double> G, Matrix2D<double> Sources, Matrix2D<double> Losses, double dt);


#endif /* DIFFUSION_2D_H_ */
