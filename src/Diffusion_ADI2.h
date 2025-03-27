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

bool Diffusion_2D_ADI2(
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


bool MakeModelMatrix_2D_ADI2_x(
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

bool MakeModelMatrix_2D_ADI2_y(
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

#endif /* DIFFUSION_ADI2_H_ */
