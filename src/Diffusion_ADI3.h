/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * \file Diffusion_ADI3.h
 *
 * \brief Calculating the diffusion in 2D using method from Jihye Shin and Sungsoo Kim (2008)
 *
 *  Created on: May 28, 2011
 *      Author: dimath
 */

#ifndef DIFFUSION_ADI3_H_
#define DIFFUSION_ADI3_H_

#include "Matrix.h"
#include "BoundaryConditionType.hpp"

/** Create model matrices and solve the system,
 * The method is from
 * Jihye Shin and Sungsoo Kim (2008), ALTERNATING DIRECTION IMPLICIT METHOD FOR	TWO-DIMENSIONAL FOKKER-PLANCK EQUATION OF DENSE SPHERICAL STELLAR SYSTEMS
 *
 * http://arxiv.org/pdf/0805.0054v1.pdf
 *
 * \todo Fix this function - currently results in negative PSD values - using Lapack instead
 *
 * Method:
 *
 * For the x direction and then for the y direction
 *
 * 1. Add boundary conditions AddBoundary()
 *
 * 2. Add sources and Losses (only for the y direction)
 *
 * 3. get the seond derivative approximation with diffusion coeficient SecondDerivativeApproximation_2D()
 *
 * 4. Multiply matrix B by f (PSD_1D_x) and add matrix C - this will be the RHS for tridag()
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

bool Diffusion_2D_ADI3(
	Matrix2D<double> &psd,
	Matrix2D<double> x, Matrix2D<double> y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> Dxx, Matrix2D<double> Dyy, Matrix2D<double> Dxy, Matrix2D<double> Dyx,
	Matrix2D<double> G, Matrix2D<double> Sources, Matrix2D<double> Losses, double dt);

bool MakeModelMatrix_2D_ADI3_x(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	Matrix2D<double> &x, Matrix2D<double> &y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> &Dxx, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
	Matrix2D<double> &G, double dt);

bool MakeModelMatrix_2D_ADI3_y(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	Matrix2D<double> &x, Matrix2D<double> &y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> &Dyy,
	Matrix2D<double> &G, double dt);
#endif /* DIFFUSION_ADI3_H_ */
