/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * \file Convection_2D.h
 *
 * \brief Two-dimensional convection solver using dimensional splitting
 *
 * This file implements a 2D convection solver using dimensional splitting with the
 * ULTIMATE QUICKEST6 scheme. The 2D advection problem is solved by applying 1D
 * convection operators sequentially in each coordinate direction.
 *
 * **Mathematical Background:**
 * The 2D advection equation:
 * \f[
 * \frac{\partial f}{\partial t} + \frac{\partial}{\partial P}(V_P \cdot f) + \frac{\partial}{\partial R}(V_R \cdot f) = S - L \cdot f
 * \f]
 * 
 * is solved using operator splitting into 1D problems.
 *
 * **References:**
 * - Leonard, B.P. "The ULTIMATE conservative difference scheme"
 * - Equations: http://www.hadian.ir/teaching/CompHydr/3.pdf
 *
 * Created on: Aug 22, 2011
 * Author: subbotin
 */

#pragma once

#include "Convection_1D_ULTIMATE_QUICKEST6.h"
#include "Matrix.h"
#include "MatrixSolver.h"
#include "BoundaryConditionType.hpp"

/**
 * @brief Solves 2D advection equation using dimensional splitting with ULTIMATE QUICKEST6
 * 
 * @details This function solves the 2D advection equation:
 * 
 * \f[
 * \frac{\partial f}{\partial t} + \frac{\partial}{\partial P}(V_P \cdot f) + \frac{\partial}{\partial R}(V_R \cdot f) = S - L \cdot f
 * \f]
 * Calculation of 2D convection
 *
 * Uses Convection_1D_Ultimate_QUICKEST6.h for calculating P and R if either are of size > 3.
 *
 * All equations and formulas for these calculations can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf .
 * Mostly coming from 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme.
 * In the source code mathematical equations are numbered corresponding to the numbering of equations found in the paper.
 * Refer to the equations/function in Convection_1D_Ultimate_QUICKEST6.h
 * 
 * @param PSD_PR [in,out] Phase space density matrix (P_size × R_size)
 * @param PSD_lost_PR [out] Lost PSD tracking matrix (only for SAVE_PSD_LOST_CONV version)
 * @param P [in] P-coordinate grid matrix (P_size × R_size)
 * @param R [in] R-coordinate grid matrix (P_size × R_size)
 * @param P_size [in] Number of grid points in P-direction
 * @param R_size [in] Number of grid points in R-direction
 * @param P_LBC [in] Lower P boundary condition values (size R_size)
 * @param P_UBC [in] Upper P boundary condition values (size R_size)
 * @param R_LBC [in] Lower R boundary condition values (size P_size)
 * @param R_UBC [in] Upper R boundary condition values (size P_size)
 * @param P_LBC_type [in] Lower P boundary condition type
 * @param P_UBC_type [in] Upper P boundary condition type
 * @param R_LBC_type [in] Lower R boundary condition type
 * @param R_UBC_type [in] Upper R boundary condition type
 * @param VP [in] Velocity field in P-direction (P_size × R_size)
 * @param VR [in] Velocity field in R-direction (P_size × R_size)
 * @param Sources [in] Source term matrix (P_size × R_size)
 * @param Losses [in] Loss coefficient matrix (P_size × R_size)
 * @param G_conv [in] Jacobian/metric factor for convection (P_size × R_size)
 * @param dt_total [in] Total time step size
 * 
 * @return true on successful completion
 * 
 * @note Grid must satisfy minimum size requirements for QUICKEST6 stencil
 * @note CFL condition must be satisfied in both directions
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

