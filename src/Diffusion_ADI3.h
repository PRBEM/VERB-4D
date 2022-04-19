/**
 * \file Diffusion_ADI3.h
 *
 * \brief Calculating the diffusion in 2D using method from Jihye Shin and Sungsoo Kim (2008)
 *
 *  Created on: May 28, 2011
 *      Author: dimath
 */

#ifndef DIFFUSION_ADI3_H_
#define DIFFUSION_ADI3_H_

#include "Matrix.h"
#include "MatrixSolver.h"


bool Diffusion_2D_ADI3(
	Matrix2D<double> &psd,
	Matrix2D<double> x, Matrix2D<double> y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	std::string x_LBC_type, std::string x_UBC_type,
	std::string y_LBC_type, std::string y_UBC_type,
	Matrix2D<double> Dxx, Matrix2D<double> Dyy, Matrix2D<double> Dxy, Matrix2D<double> Dyx,
	Matrix2D<double> G, Matrix2D<double> Sources, Matrix2D<double> Losses, double dt
);

#endif /* DIFFUSION_ADI3_H_ */
