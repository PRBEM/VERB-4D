/**
 * \file Convection_3D.h
 *
 * Using the Convection_1D_Ultimate_QUICKEST6_COULOMB() for calculating P and R.
 * All equations and formulas for these calculations can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf .
 * Mostly coming from 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme.
 * In the source code mathematical equations are numbered corresponding to the numbering of equations found in the paper
 *
 * */

#pragma once

#include "Matrix.h"
#include "MatrixSolver.h"
#include "BoundaryConditionType.hpp"

#include "Convection_1D_ULTIMATE_QUICKEST6.h"
#include "Convection_1D_2ndORDER_NONUNIFORM_GRID.h"

/*
* Function that creates a 3D Convection matrix and returns a bool upon completion
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

