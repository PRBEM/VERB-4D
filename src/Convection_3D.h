/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * \file Convection_3D.h
 *
 * \brief Three-dimensional convection solver using dimensional splitting
 *
 * This file implements a 3D convection solver for phase space transport using dimensional
 * splitting with high-order finite difference schemes. The solver handles transport in
 * P-R-V phase space commonly used in plasma physics and space weather modeling.
 *
 * **Mathematical Background:**
 * The 3D advection equation in phase space:
 * \f[
 * \frac{\partial f}{\partial t} + \frac{\partial}{\partial P}(V_P \cdot f) + 
 * \frac{\partial}{\partial R}(V_R \cdot f) + \frac{\partial}{\partial V}(V_V \cdot f) = S - L \cdot f
 * \f]
 * 
 * **References:**
 * - Leonard, B.P. "The ULTIMATE conservative difference scheme"
 * - Equations: http://www.hadian.ir/teaching/CompHydr/3.pdf
 *
 * Created on: May 28, 2011
 * Author: dimath
 */

#pragma once

#include "Matrix.h"
#include "MatrixSolver.h"
#include "BoundaryConditionType.hpp"
#include "Convection_1D_ULTIMATE_QUICKEST6.h"
#include "Convection_1D_2ndORDER_NONUNIFORM_GRID.h"

/**
 * @brief Solves 3D advection equation using dimensional splitting
 * 
 * @details This function solves the 3D phase space advection equation:
 * 
 * \f[
 * \frac{\partial f}{\partial t} + \frac{\partial}{\partial P}(V_P \cdot f) + 
 * \frac{\partial}{\partial R}(V_R \cdot f) + \frac{\partial}{\partial V}(V_V \cdot f) = S - L \cdot f
 * \f]
 * 
 * where:
 * - \f$ f(P,R,V,t) \f$ is the phase space density
 * - \f$ P, R \f$ are spatial coordinates (e.g., radial distance, magnetic local time)
 * - \f$ V \f$ is velocity coordinate (e.g., particle energy, pitch angle)
 * - \f$ V_P, V_R, V_V \f$ are velocity fields in each direction
 * - \f$ S, L \f$ are source and loss terms
 * 
 * Calculation of 3D convection
 *
 * Uses Convection_1D_Ultimate_QUICKEST6.h for calculating P and R if either are of size > 3.
 * Uses Convection_1D_2ndORDER_NONUNIFORM_GRID.h for calculating V if size > 3.
 *
 * All equations and formulas for these calculations can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf .
 * Mostly coming from 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme.
 * In the source code mathematical equations are numbered corresponding to the numbering of equations found in the paper.
 * Refer to the equations/function in `Convection_1D_Ultimate_QUICKEST6.h`
 *
 * 
 * **Dimensional Splitting Algorithm:**
 * 
 * The 3D problem is decomposed into three sequential 1D problems:
 * 
 * 1. **P-direction sweep:** Apply ULTIMATE QUICKEST6 or 2nd-order scheme
 * 2. **R-direction sweep:** Apply ULTIMATE QUICKEST6 or 2nd-order scheme
 * 3. **V-direction sweep:** Apply appropriate 1D convection scheme with source/loss terms
 * 
 * 
 * @param PSD_PRV [in,out] 3D phase space density (P_size × R_size × V_size)
 * @param P [in] P-coordinate grid (P_size × R_size × V_size)
 * @param R [in] R-coordinate grid (P_size × R_size × V_size)
 * @param V [in] V-coordinate grid (P_size × R_size × V_size)
 * @param P_size [in] Number of grid points in P-direction
 * @param R_size [in] Number of grid points in R-direction
 * @param V_size [in] Number of grid points in V-direction
 * @param P_LBC [in] Lower P boundary conditions (R_size × V_size)
 * @param P_UBC [in] Upper P boundary conditions (R_size × V_size)
 * @param R_LBC [in] Lower R boundary conditions (P_size × V_size)
 * @param R_UBC [in] Upper R boundary conditions (P_size × V_size)
 * @param V_LBC [in] Lower V boundary conditions (P_size × R_size)
 * @param V_UBC [in] Upper V boundary conditions (P_size × R_size)
 * @param P_LBC_type [in] Lower P boundary condition type
 * @param P_UBC_type [in] Upper P boundary condition type
 * @param R_LBC_type [in] Lower R boundary condition type
 * @param R_UBC_type [in] Upper R boundary condition type
 * @param V_LBC_type [in] Lower V boundary condition type
 * @param V_UBC_type [in] Upper V boundary condition type
 * @param VP [in] Velocity field in P-direction (P_size × R_size × V_size)
 * @param VR [in] Velocity field in R-direction (P_size × R_size × V_size)
 * @param VV [in] Velocity field in V-direction (P_size × R_size × V_size)
 * @param Sources [in] Source terms (P_size × R_size × V_size)
 * @param Losses [in] Loss coefficients (P_size × R_size × V_size)
 * @param G_conv [in] Jacobian/metric factors (P_size × R_size × V_size)
 * @param dt [in] Time step size
 * 
 * @return true on successful completion
 * 
 * @note CFL condition must be satisfied in all three directions
 */

bool Convection_3D( 
		Matrix3D < double >& PSD_PRV,
		const Matrix3D < double >& P, 
		const Matrix3D < double >& R, 
		const Matrix3D < double >& V,
		int P_size, int R_size, int V_size,
		const Matrix2D < double >& P_LBC, const Matrix2D < double >& P_UBC,
		const Matrix2D < double >& R_LBC, const Matrix2D < double >& R_UBC,
		const Matrix2D < double >& V_LBC, const Matrix2D < double >& V_UBC,
		BoundaryConditionType P_LBC_type, BoundaryConditionType P_UBC_type,
		BoundaryConditionType R_LBC_type, BoundaryConditionType R_UBC_type,
		BoundaryConditionType V_LBC_type, BoundaryConditionType V_UBC_type,
		const Matrix3D < double >& VP, const Matrix3D < double >& VR, const Matrix3D < double >& VV,
		const Matrix3D < double >& Sources, const Matrix3D < double >& Losses, 
		const Matrix3D<double>& G_conv, double dt);

