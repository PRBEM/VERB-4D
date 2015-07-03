/**
* \file Convection_1D_ULTIMATE_QUICKEST6.h
*
* \brief Calculates the convection in 1D given a 1D matrix of Phase Space Densities, boundary conditions, diffusion, sources and losses
*/

#ifndef CONVECTION_1D_ULTIMATE_QUICKEST6_H_
#define CONVECTION_1D_ULTIMATE_QUICKEST6_H_





#include "Matrix.h"
#include "MatrixSolver.h"
#include "Logger.h"

/**
* Function that creates a 1D Convection Quickset and returns a bool upon completion
* \param PSD - 1D matrix of phase space densities
* \param x - matrix for determining step size dx = x[1] - x[0]
* \param x_size - size of matrix x
* \param x_LBC, x_UBC, _type - The upper and lower boundary conditions including the type of boundary condition
* \param Ux - diffusion matrix
* \param Sources - Sources matrix
* \param Losses - Losses (loss cone)
*/
bool Convection_1D_ULTIMATE_QUICKEST6 ( Matrix1D < double > &PSD,
		Matrix1D < double > x,
		int x_size,
		double x_LBC, double x_UBC,
		string x_LBC_type, string x_UBC_type,
		Matrix1D < double > Ux,
		Matrix1D < double > Sources, Matrix1D < double > Losses,
		double dt);


#endif
