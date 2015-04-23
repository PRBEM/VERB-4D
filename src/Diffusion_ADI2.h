/*
 * Diffusion.h
 *
 *  Created on: May 28, 2011
 *      Author: dimath
 */

#ifndef DIFFUSION_ADI2_H_
#define DIFFUSION_ADI2_H_

#include "Matrix.h"
#include "MatrixSolver.h"


bool Diffusion_2D_ADI2(
				  Matrix2D<double> &psd,
				  Matrix2D<double> x, Matrix2D<double> y,
				  int x_size, int y_size,
				  Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
				  Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
				  string x_LBC_type, string x_UBC_type,
				  string y_LBC_type, string y_UBC_type,
				  Matrix2D<double> Dxx, Matrix2D<double> Dyy, Matrix2D<double> Dxy, Matrix2D<double> Dyx,
				  Matrix2D<double> G, Matrix2D<double> Sources, Matrix2D<double> Losses, double dt);

#endif /* DIFFUSION_ADI2_H_ */
