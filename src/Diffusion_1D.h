/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * \file Diffusion_1D.h
 *
 * \brief One-dimensional diffusion solver for phase space density evolution
 *
 * This file implements a finite difference solver for the one-dimensional diffusion equation
 * with sources and losses, commonly used in plasma physics and space weather modeling.
 *
 * Created on: May 28, 2011
 * Author: dimath
 */

#ifndef DIFFUSION_1D_H_
#define DIFFUSION_1D_H_

#include "Matrix.h"
#include "MatrixSolver.h"
#include "BoundaryConditionType.hpp"

/**
 * @brief Solves the one-dimensional diffusion equation with sources and losses
 * 
 * @details This function implements an implicit finite difference scheme to solve the 1D diffusion equation:
 * 
 * \f[
 * \frac{\partial f}{\partial t} = \frac{1}{G} \frac{\partial}{\partial x} \left( G D_{xx} \frac{\partial f}{\partial x} \right) + S - L \cdot f
 * \f]
 * 
 * where:
 * - \f$ f \f$ is the phase space density (psd)
 * - \f$ G \f$ is the Jacobian (metric factor)
 * - \f$ D_{xx} \f$ is the diffusion coefficient
 * - \f$ S \f$ represents source terms
 * - \f$ L \f$ represents loss terms
 * 
 * The equation is discretized with implicit treatment of the diffusion operator for numerical stability. The resulting tridiagonal system is solved using the Thomas algorithm (tridag).
 * 
 * **Mathematical Method:**
 * 1. Apply boundary conditions (Dirichlet or Neumann)
 * 2. Include source and loss terms
 * 3. Discretize the diffusion operator using second-order finite differences `SecondDerivativeApproximation_1D()`
 * 4. Solve the resulting tridiagonal linear system with `tridag()`
 * 
 * @param psd [in,out] Phase space density vector to be updated (size x_size)
 * @param x [in] Spatial coordinate grid points (size x_size)
 * @param x_size [in] Number of grid points in the x-direction
 * @param x_LBC [in] Lower boundary condition value for x
 * @param x_UBC [in] Upper boundary condition value for x
 * @param x_LBC_type [in] Type of lower boundary condition (ConstantValue or ConstantDerivative) for x
 * @param x_UBC_type [in] Type of upper boundary condition (ConstantValue or ConstantDerivative) for x
 * @param Dxx [in] Diffusion coefficient at each grid point (size x_size) for x
 * @param G [in] Jacobian (metric factor) at each grid point (size x_size)
 * @param Sources [in] Source terms at each grid point (size x_size)
 * @param Losses [in] Loss coefficients at each grid point (size x_size)
 * @param dt [in] Time step size
 * 
 * @return true on successful completion, false on error
 * 
 * @note Requires x_size >= 3 for proper finite difference stencil
 * @warning The ConstantDerivative boundary condition currently only works correctly for zero derivative
 */
bool Diffusion_1D(
	Matrix1D<double>& psd,
	const Matrix1D<double>& x, int x_size,
	double x_LBC, double x_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	const Matrix1D<double>& Dxx, const Matrix1D<double>& G,
	const Matrix1D<double>& Sources, const Matrix1D<double>& Losses,
	double dt
);

std::tuple<CalculationMatrix, CalculationMatrix, CalculationMatrix> getCalculationMatrix1DDiffusion(
	const Matrix1D<double>& x, int x_size,
	double x_LBC, double x_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	const Matrix1D<double>& Dxx, const Matrix1D<double>& G,
	const Matrix1D<double>& Sources, const Matrix1D<double>& Losses,
	double dt

);

#endif /* DIFFUSION_1D_H_ */
