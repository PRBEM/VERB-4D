/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * \file Convection_2D.h
 *
 * Using the Convection_1D_Ultimate_QUICKEST6() for calculating P and R.
 * All equations and formulas for these calculations can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf .
 * Mostly coming from 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme.
 * In the source code mathematical equations are numbered corresponding to the numbering of equations found in the paper
 *
 *
 * \brief Calculates the convection in 2D given a 2D matrix of Phase Space Densities, P, R, boundary conditions, diffusion, sources and losses
 *
 *  Created on: Aug 22, 2011
 *      Author: subbotin
 */

#pragma once

#include "Convection_1D_ULTIMATE_QUICKEST6.h"
#include "Matrix.h"
#include "MatrixSolver.h"

#include "Convection_1D_ULTIMATE_QUICKEST6.h"
#include "BoundaryConditionType.hpp"

/**
* Function that creates a 2D Convection matrix and returns a bool upon completion
*/

#ifdef SAVE_PSD_LOST_CONV
bool Convection_2D(
    Matrix2D<double>& PSD_PR,
    Matrix2D<double>& PSD_lost_PR,
    const Matrix2D<double>& P, const Matrix2D<double>& R,
    int P_size, int R_size,
    const Matrix1D<double>& P_LBC, const Matrix1D<double>& P_UBC,
    const Matrix1D<double>& R_LBC, const Matrix1D<double>& R_UBC,
    BoundaryConditionType P_LBC_type, BoundaryConditionType P_UBC_type,
    BoundaryConditionType R_LBC_type, BoundaryConditionType R_UBC_type,
    const Matrix2D<double>& VP, const Matrix2D<double>& VR,
    const Matrix2D<double>& Sources, const Matrix2D<double>& Losses,
    const Matrix2D<double>& G_conv, double dt_total);
#else
bool Convection_2D( 
	Matrix2D<double>& PSD_PR,
	const Matrix2D<double>& P, const Matrix2D<double>& R,
	int P_size, int R_size,
	const Matrix1D<double>& P_LBC, const Matrix1D<double>& P_UBC,
	const Matrix1D<double>& R_LBC, const Matrix1D<double>& R_UBC,
	BoundaryConditionType P_LBC_type, BoundaryConditionType P_UBC_type,
	BoundaryConditionType R_LBC_type, BoundaryConditionType R_UBC_type,
	const Matrix2D<double>& VP, const Matrix2D<double>& VR,
	const Matrix2D<double>& Sources, const Matrix2D<double>& Losses,
	const Matrix2D<double>& G_conv, double dt_total);
#endif

