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
#include "BoundaryConditionType.hpp"

/**
* Function that creates a 2D Diffusion Matrix and returns a bool upon completion
*/
bool Diffusion_2D(
	Matrix2D<double> &psd,
	const Matrix2D<double>& x, const Matrix2D<double>& y,
	int x_size, int y_size,
	const Matrix1D<double>& x_LBC, const Matrix1D<double>& x_UBC,
	const Matrix1D<double>& y_LBC, const Matrix1D<double>& y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	const Matrix2D<double>& Dxx, const Matrix2D<double>& Dyy, const Matrix2D<double>& Dxy, const Matrix2D<double>& Dyx,
	const Matrix2D<double>& G, const Matrix2D<double>& Sources, const Matrix2D<double>& Losses, double dt
);


#endif /* DIFFUSION_2D_H_ */
