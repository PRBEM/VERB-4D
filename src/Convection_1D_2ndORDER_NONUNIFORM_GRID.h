/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
* \file Convection_1D_2ndORDER_NONUNIFORM_GRID.h
*
* All equations and formulas for these calculations can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf .
 * Mostly coming from 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme.
 * In the source code mathematical equations are numbered corresponding to the numbering of equations found in the paper
 * 
* \brief Calculates the convection in 1D given a 1D matrix of Phase Space Densities, boundary conditions, diffusion, sources and losses
*/

#pragma once

#include "Matrix.h"
#include "MatrixSolver.h"
#include "Logger.h"
#include "BoundaryConditionType.hpp" 

/**
* Function that implements 1D Convection Quickest and returns a bool upon completion as referenced in the following publication
*
* All equations and formulas for this function can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf
* Mostly coming from 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme.
* In the source code mathematical equations are numbered corresponding to the numbering of equations found in the paper
* 
* \f$ c = u*dt/dx \f$ where c = Courant number
*
* \f$ \phi_{i}^{n+1} - \phi_{i}^{n} = -c(\phi_{r}*\phi_{l}) \f$ where \f$ \phi_{r} \f$ and \f$ \phi_{l} \f$ refer to the right and left faces (between (i and i-1), and (i and i+1)) for PSD
*
* \param PSD - 1D matrix of phase space densities
* \param PSDtimesLoss - 1D matrix of phase space densities times the Coulomb losses
* \param x - matrix for determining step size dx = x[1] - x[0]
* \param x_size - size of matrix x
* \param x_LBC, x_UBC, _type - The upper and lower boundary conditions including the type of boundary condition
* \param Ux - diffusion matrix
* \param Sources - Sources matrix
* \param Losses - Losses (loss cone)
*/
bool Convection_1D_2ndORDER_NONUNIFORM_GRID ( 
		Matrix1D < double >& PSD,
		const Matrix1D < double >& x,
		int x_size,
		double x_LBC, double x_UBC,
		BoundaryConditionType x_LBC_type, 
		BoundaryConditionType x_UBC_type,
		const Matrix1D < double >& Ux,
		double dt_total);