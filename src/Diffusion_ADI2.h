/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * \file Diffusion_ADI2.h
 *
 * \brief Calculating the diffusion in 2D using multidimensional diffusion process - Xiao, F., Z. Su, H. Zheng, and S. Wang (2009)
 *
 *  Created on: May 28, 2011
 *      Author: dimath
 */

#ifndef DIFFUSION_ADI2_H_
#define DIFFUSION_ADI2_H_

#include "Matrix.h"
#include "BoundaryConditionType.hpp"

/** Create model matrices and solve the system,
 * The method is from
 * Xiao, F., Z. Su, H. Zheng, and S. Wang (2009), Modeling of outer radiation belt electrons by multidimensional diffusion process, J. Geophys. Res., 114, A03201, doi:10.1029/2008JA013580.
 *
 * http://onlinelibrary.wiley.com/store/10.1029/2008JA013580/asset/jgra19599.pdf?v=1&t=ic3m61g2&s=334d50627dee2d592de0ae583e755afce42948a9
 *
 * \todo Fix this function - currently results in negative PSD values - using Lapack instead
 *
 * Method:
 *
 * For the x direction and then for the y direction
 *
 * 1. Add boundary conditions AddBoundary()
 *
 * 2. Add sources and losses
 *
 * 3. get the seond derivative approximation with diffusion coeficient SecondDerivativeApproximation_2D()
 *
 * 4. Implicit part (first derivatives)	\f$ \frac{-1}{G} \frac {\partial{(G D_{xy})}}{\partial x} \frac {\partial f}{\partial y} \f$
 *
 * 5. solve matrix with tridag()
 *
 * @param psd [in,out] - phase space density
 * @param x [in] - one dimensional slice
 * @param y [in] - one dimensional slice
 * @param x_size [in] - size of x slice
 * @param y_size [in] - size of y slice
 * @param x_LBC [in] - lower boundary condition for param x
 * @param x_UBC [in] - upper boundary condition for param x
 * @param y_LBC [in] - lower boundary condition for param y
 * @param y_UBC [in] - upper boundary condition for param y
 * @param x_LBC_type [in] - type of lower boundary condition for param x
 * @param x_UBC_type [in] - type of upper boundary condition for param x
 * @param y_LBC_type [in] - type of lower boundary condition for param y
 * @param y_UBC_type [in] - type of upper boundary condition for param y
 * @param Dxx [in] - 2D Diffusion matrix
 * @param Dyy [in] - 2D Diffusion matrix
 * @param Dxy [in] - 2D Diffusion matrix
 * @param Dyx [in] - 2D Diffusion matrix
 * @param G [in] - 2D used for Jacobian to normalize matrix
 * @param Sources [in] - matrix used for Sources
 * @param Losses [in] - Matrix used for Losses (loss cone)
 * @param dt [in] - change in time of single time step
 */

bool Diffusion_2D_ADI2(
	Matrix2D<double> &psd,
	Matrix2D<double> x, Matrix2D<double> y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> Dxx, Matrix2D<double> Dyy, Matrix2D<double> Dxy, Matrix2D<double> Dyx,
	Matrix2D<double> G, Matrix2D<double> Sources, Matrix2D<double> Losses, double dt);

bool MakeModelMatrix_2D_ADI2_x(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	Matrix2D<double> &x, Matrix2D<double> &y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> &Dxx, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
	Matrix2D<double> &G, double dt);

bool MakeModelMatrix_2D_ADI2_y(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	Matrix2D<double> &x, Matrix2D<double> &y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
	Matrix2D<double> &G, double dt);

#endif /* DIFFUSION_ADI2_H_ */
