/**
 * \file Diffusion_2D.h
 *
 * \brief Getting the diffusion in 2D for the given a 2D matrix of Phase Space Densities, x, y, boundary conditions, jacobian, sources, losses and change in time
 *
 *  Created on: May 28, 2011
 *      Author: dimath
 */

#ifndef DIFFUSION_2D_H_
#define DIFFUSION_2D_H_

#include "Matrix.h"
#include "MatrixSolver.h"

/**
* Function that creates a 2D Diffusion Matrix and returns a bool upon completion
*/
bool Diffusion_2D(
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


#endif /* DIFFUSION_2D_H_ */
