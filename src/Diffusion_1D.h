/**
 * \file Diffusion_1D.h
 *
 * \brief Getting the diffusion in 1D for the given a 1D matrix of Phase Space Densities, x, boundary conditions, jacobian, sources, losses and change in time
 *
 *  Created on: May 28, 2011
 *      Author: dimath
 */

#ifndef DIFFUSION_1D_H_
#define DIFFUSION_1D_H_

#include "Matrix.h"
#include "MatrixSolver.h"

/**
* Function that creates a 1D Diffusion Matrix and returns a bool upon completion
*/
bool Diffusion_1D(
				  Matrix1D<double> &psd,
				  Matrix1D<double> x, int x_size,
				  double x_LBC, double x_UBC,
				  string x_LBC_type, string x_UBC_type,
				  Matrix1D<double> Dxx, Matrix1D<double> G,
				  Matrix1D<double> Sources, Matrix1D<double> Losses,
				  double dt);

#endif /* DIFFUSION_1D_H_ */
