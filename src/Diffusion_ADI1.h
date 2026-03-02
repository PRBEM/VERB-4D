/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * \file Diffusion_ADI1.h
 *
 * \brief Two-dimensional diffusion solver using Alternating Direction Implicit (ADI) method
 *
 * This file implements the ADI method for solving 2D diffusion equations. The ADI approach
 * splits the 2D problem into a sequence of 1D problems, making it computationally efficient
 * for large grids while maintaining good stability properties.
 *
 * Created on: May 28, 2011
 * Author: dimath
 */

#ifndef DIFFUSION_ADI1_H_
#define DIFFUSION_ADI1_H_

#include "Matrix.h"
#include "BoundaryConditionType.hpp"

/**
 * @brief Solves 2D diffusion equation using Alternating Direction Implicit (ADI) method
 * 
 * @details This function implements the ADI method for the 2D Fokker-Planck equation:
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
 * **ADI Method Overview:**
 * For the x direction and then for the y direction
 * 1. Add boundary conditions `AddBoundary()`
 *
 * 2. Add sources and losses
 *
 * 3. get the seond derivative approximation with diffusion coeficient `SecondDerivativeApproximation_2D()`
 *
 * 4. make `RHS = B*f + C` and then solve with `tridag()`

 * @param psd [in,out] Phase space density matrix to be updated
 * @param x [in] X-coordinate grid (x_size × y_size)
 * @param y [in] Y-coordinate grid (x_size × y_size)
 * @param x_size [in] Number of grid points in x-direction
 * @param y_size [in] Number of grid points in y-direction
 * @param x_LBC [in] Lower x boundary condition values (size y_size)
 * @param x_UBC [in] Upper x boundary condition values (size y_size)
 * @param y_LBC [in] Lower y boundary condition values (size x_size)
 * @param y_UBC [in] Upper y boundary condition values (size x_size)
 * @param x_LBC_type [in] Type of lower x boundary condition
 * @param x_UBC_type [in] Type of upper x boundary condition
 * @param y_LBC_type [in] Type of lower y boundary condition
 * @param y_UBC_type [in] Type of upper y boundary condition
 * @param Dxx [in] X-direction diffusion coefficient matrix
 * @param Dyy [in] Y-direction diffusion coefficient matrix
 * @param Dxy [in] Cross-diffusion coefficient matrix (x→y coupling)
 * @param Dyx [in] Cross-diffusion coefficient matrix (y→x coupling)
 * @param G [in] Jacobian/metric factor matrix
 * @param Sources [in] Source term matrix
 * @param Losses [in] Loss coefficient matrix
 * @param dt [in] Time step size
 * 
 * @return true on successful completion
 * 
 * @warning This implementation has known issues that can produce negative PSD values, using Lapack instead
 */
bool Diffusion_2D_ADI1(
	Matrix2D<double> &psd,
	Matrix2D<double> x, Matrix2D<double> y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> Dxx, Matrix2D<double> Dyy, Matrix2D<double> Dxy, Matrix2D<double> Dyx,
	Matrix2D<double> G, Matrix2D<double> Sources, Matrix2D<double> Losses, double dt
);


/**
 * @brief Creates model matrices for ADI x-direction sweep (currently unused)
 * 
 * @details This function constructs the coefficient matrices for the x-direction
 * step of the ADI method. It treats x-derivatives implicitly and y-derivatives explicitly.
 * 
 * @param matr_A [out] Coefficient matrix for implicit terms
 * @param matr_B [out] Coefficient matrix for explicit terms  
 * @param matr_C [out] Constant term matrix
 * @param x [in] X-coordinate grid
 * @param y [in] Y-coordinate grid
 * @param x_size [in] Grid size in x-direction
 * @param y_size [in] Grid size in y-direction
 * @param x_LBC [in] Lower x boundary conditions
 * @param x_UBC [in] Upper x boundary conditions
 * @param y_LBC [in] Lower y boundary conditions
 * @param y_UBC [in] Upper y boundary conditions
 * @param x_LBC_type [in] Lower x boundary condition type
 * @param x_UBC_type [in] Upper x boundary condition type
 * @param y_LBC_type [in] Lower y boundary condition type
 * @param y_UBC_type [in] Upper y boundary condition type
 * @param Dxx [in] X-direction diffusion coefficients
 * @param Dxy [in] Cross-diffusion coefficients
 * @param Dyx [in] Cross-diffusion coefficients
 * @param G [in] Jacobian factors
 * @param dt [in] Time step size
 * 
 * @return true on successful completion
 * 
 * @note This function is currently not being used in the main ADI implementation
 */
bool MakeModelMatrix_2D_ADI1_x(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	Matrix2D<double> &x, Matrix2D<double> &y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> &Dxx, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
	Matrix2D<double> &G, double dt
);

/**
 * @brief Creates model matrices for ADI y-direction sweep (currently unused)
 * 
 * @details This function constructs the coefficient matrices for the y-direction
 * step of the ADI method. It treats y-derivatives implicitly and x-derivatives explicitly.
 * 
 * @param matr_A [out] Coefficient matrix for implicit terms
 * @param matr_B [out] Coefficient matrix for explicit terms
 * @param matr_C [out] Constant term matrix
 * @param x [in] X-coordinate grid
 * @param y [in] Y-coordinate grid
 * @param x_size [in] Grid size in x-direction
 * @param y_size [in] Grid size in y-direction
 * @param x_LBC [in] Lower x boundary conditions
 * @param x_UBC [in] Upper x boundary conditions
 * @param y_LBC [in] Lower y boundary conditions
 * @param y_UBC [in] Upper y boundary conditions
 * @param x_LBC_type [in] Lower x boundary condition type
 * @param x_UBC_type [in] Upper x boundary condition type
 * @param y_LBC_type [in] Lower y boundary condition type
 * @param y_UBC_type [in] Upper y boundary condition type
 * @param Dyy [in] Y-direction diffusion coefficients
 * @param Dxy [in] Cross-diffusion coefficients
 * @param Dyx [in] Cross-diffusion coefficients
 * @param G [in] Jacobian factors
 * @param dt [in] Time step size
 * 
 * @return true on successful completion
 * 
 * @note This function is currently not being used in the main ADI implementation
 */
bool MakeModelMatrix_2D_ADI1_y(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	Matrix2D<double> &x, Matrix2D<double> &y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
	Matrix2D<double> &G, double dt
);


#endif /* DIFFUSION_ADI1_H_ */
