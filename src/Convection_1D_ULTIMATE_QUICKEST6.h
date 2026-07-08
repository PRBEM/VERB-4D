/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * \file Convection_1D_ULTIMATE_QUICKEST6.h
 *
 * \brief One-dimensional convection solver using ULTIMATE QUICKEST6 scheme
 *
 * This file implements the ULTIMATE QUICKEST6 scheme for solving the 1D advection equation.
 * The method combines high-order accuracy with monotonicity preservation through flux limiting.
 * 
 * **References:**
 * - Equations and formulas: http://www.hadian.ir/teaching/CompHydr/3.pdf
 * - Leonard, 1991; Leonard and Niknafs, 1991;
 * - Leonard BP (1988) Universal Limiter for transient interpolation modeling of the advective transport equations: the ULTIMATE conservative difference scheme, NASA technical Memorandum 100916 ICOMP-88-11
 * 
 * It works somehow, edit with a great care! - edited by Hayley to correct and issue with the PSD update. Without this correction, the total PSD will gain and drop periodically with the drift
 * Mostly coming from 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme.
 * In the source code mathematical equations are numbered corresponding to the numbering of equations found in the paper
 *
 *
 * Created on: May 28, 2011
 * Author: dimath
 */

#ifndef CONVECTION_1D_ULTIMATE_QUICKEST6_H_
#define CONVECTION_1D_ULTIMATE_QUICKEST6_H_

#include "Matrix.h"
#include "MatrixSolver.h"
#include "Logger.h"
#include "BoundaryConditionType.hpp"

/**
 * @brief Solves 1D advection equation using ULTIMATE QUICKEST6 scheme
 * 
 * @details This function implements the ULTIMATE QUICKEST6 scheme for the 1D advection equation:
 * 
 * \f[
 * \frac{\partial f}{\partial t} + \frac{\partial}{\partial x}(u \cdot f) = 0
 * \f]
 * 
 * **ULTIMATE QUICKEST6 Method:**
 * 
 * \f$ c = u*dt/dx \f$ where c = Courant number
 * 
 * The discrete update formula is:
 * 
 * \f$ \phi_{i}^{n+1} - \phi_{i}^{n} = -c(\phi_{r}*\phi_{l}) \f$ where \f$ \phi_{r} \f$ and \f$ \phi_{l} \f$ refer to the right and left faces (between (i and i-1), and (i and i+1)) for PSD
 * 
 * 
 * @param PSD [in,out] Phase space density vector (size x_size)
 * @param x [in] Spatial coordinate grid (size x_size, assumed uniform)
 * @param x_size [in] Number of grid points (must be ≥ 11 for full stencil)
 * @param x_LBC [in] Lower boundary condition value
 * @param x_UBC [in] Upper boundary condition value
 * @param x_LBC_type [in] Lower boundary condition type (Periodic/ConstantValue/ConstantDerivative)
 * @param x_UBC_type [in] Upper boundary condition type (Periodic/ConstantValue/ConstantDerivative)
 * @param Ux [in] Velocity field at grid points (size x_size)
 * @param dt_total [in] Total time step (automatically subdivided if needed)
 * 
 * @return true on successful completion, exits on CFL violation
 * 
 * @note Requires uniform grid spacing for current implementation
 * @note Minimum grid size is 11 points due to ghost cells and stencil
 */
bool Convection_1D_ULTIMATE_QUICKEST6 (
	Matrix1D<double> &PSD, 
	const Matrix1D<double>& x,
	int x_size,
	double x_LBC, double x_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	const Matrix1D<double>& Ux,
	double dt_total
);

#endif
