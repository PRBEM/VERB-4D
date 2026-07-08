/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * \file Diffusion_2D.h
 *
 * \brief Two-dimensional diffusion solver using direct matrix methods
 *
 * This file implements a 2D diffusion solver for the Fokker-Planck equation
 * with full tensor diffusion coefficients, including cross-derivative terms. The solver
 * uses direct matrix methods (LAPACK) for high accuracy but with higher computational cost.
 *
 * Created on: May 28, 2011
 * Author: dimath
 */

#ifndef DIFFUSION_2D_H_
#define DIFFUSION_2D_H_

#include "Matrix.h"
#include "MatrixSolver.h"
#include "BoundaryConditionType.hpp"

/**
 * @brief Solves the two-dimensional diffusion equation with full tensor coefficients
 *
 * @details This function implements an implicit finite difference scheme for the 2D Fokker-Planck equation:
 *
 * \f[
 * \frac{\partial f}{\partial t} = \frac{1}{G} \left[
 * \frac{\partial}{\partial x} \left( G D_{xx} \frac{\partial f}{\partial x} \right) +
 * \frac{\partial}{\partial y} \left( G D_{yy} \frac{\partial f}{\partial y} \right) +
 * \frac{\partial}{\partial x} \left( G D_{xy} \frac{\partial f}{\partial y} \right) +
 * \frac{\partial}{\partial y} \left( G D_{yx} \frac{\partial f}{\partial x} \right)
 * \right] + S - L \cdot f
 * \f]
 *
 * where:
 * - \f$ f(x,y,t) \f$ is the phase space density
 * - \f$ G(x,y) \f$ is the Jacobian (metric factor)
 * - \f$ D_{xx}, D_{yy} \f$ are diagonal diffusion coefficients
 * - \f$ D_{xy}, D_{yx} \f$ are cross-diffusion coefficients
 * - \f$ S(x,y) \f$ represents source terms
 * - \f$ L(x,y) \f$ represents loss coefficients
 *
 * **Mathematical Method:**
 * 1. Apply boundary conditions using AddBoundaries_2D()
 * 2. Include source and loss terms
 * 3. Discretize the diffusion operator using second-order finite differences `SecondDerivativeApproximation_2D()`
 * 4. Solve matrix with Lapack()
 *
 * @param psd [in,out] Phase space density matrix (x_size × y_size)
 * @param x [in] X-coordinate grid (x_size × y_size)
 * @param y [in] Y-coordinate grid (x_size × y_size)
 * @param x_size [in] Number of grid points in x-direction
 * @param y_size [in] Number of grid points in y-direction
 * @param x_LBC [in] Lower boundary condition values for x-direction (size y_size)
 * @param x_UBC [in] Upper boundary condition values for x-direction (size y_size)
 * @param y_LBC [in] Lower boundary condition values for y-direction (size x_size)
 * @param y_UBC [in] Upper boundary condition values for y-direction (size x_size)
 * @param x_LBC_type [in] Type of lower x boundary condition
 * @param x_UBC_type [in] Type of upper x boundary condition
 * @param y_LBC_type [in] Type of lower y boundary condition
 * @param y_UBC_type [in] Type of upper y boundary condition
 * @param Dxx [in] Diffusion coefficient matrix for \f$ \frac{\partial^2f}{\partial x^2} \f$ terms
 * @param Dyy [in] Diffusion coefficient matrix for \f$ \frac{\partial^2f}{\partial y^2} \f$ terms
 * @param Dxy [in] Cross-diffusion coefficient matrix for \f$ \frac{\partial}{\partial x}(D\frac{\partial f}{\partial y}) \f$ terms
 * @param Dyx [in] Cross-diffusion coefficient matrix for \f$ \frac{\partial}{\partial y}(D\frac{\partial f}{\partial x}) \f$ terms
 * @param G [in] Jacobian/metric factor matrix (x_size × y_size)
 * @param Sources [in] Source term matrix (x_size × y_size)
 * @param Losses [in] Loss coefficient matrix (x_size × y_size)
 * @param dt [in] Time step size
 * @param sub_dt [in] Desired substep size (at runtime `dt` will be evaluated to \f$ \frac{dt}{\lfloor\frac{dt}{sub\_dt}\rfloor} \f$)
 * @param lu_cache [in] LU decomposition cache (only for `LU_CACHING` version)
 * @param index_cache [in] Index cache for LU decomposition (only for `LU_CACHING` version)
 * @param recompute_lu [in] Whether to recompute LU decomposition (only for `LU_CACHING` version)
 *
 * @return true on successful completion, false on error
 *
 * @warning Memory usage scales as \f$\mathcal{O}(n^2)\f$ - monitor for large grids
 */
#ifdef LU_CACHING
	bool Diffusion_2D(
		Matrix2D<double> &psd,
		const Matrix2D<double>& x, const Matrix2D<double>& y,
		int x_size, int y_size,
		const Matrix1D<double>& x_LBC, const Matrix1D<double>& x_UBC,
		const Matrix1D<double>& y_LBC, const Matrix1D<double>& y_UBC,
		BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
		BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
		const Matrix2D<double>& Dxx, const Matrix2D<double>& Dyy, const Matrix2D<double>& Dxy, const Matrix2D<double>& Dyx,
		const Matrix2D<double>& G, const Matrix2D<double>& Sources, const Matrix2D<double>& Losses, double dt, double sub_dt,
		double* lu_cache, long* index_cache, bool recompute_lu
	);
#else
	bool Diffusion_2D(
		Matrix2D<double> &psd,
		const Matrix2D<double>& x, const Matrix2D<double>& y,
		int x_size, int y_size,
		const Matrix1D<double>& x_LBC, const Matrix1D<double>& x_UBC,
		const Matrix1D<double>& y_LBC, const Matrix1D<double>& y_UBC,
		BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
		BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
		const Matrix2D<double>& Dxx, const Matrix2D<double>& Dyy, const Matrix2D<double>& Dxy, const Matrix2D<double>& Dyx,
		const Matrix2D<double>& G, const Matrix2D<double>& Sources, const Matrix2D<double>& Losses, double dt, double sub_dt
	);
#endif


#endif /* DIFFUSION_2D_H_ */
