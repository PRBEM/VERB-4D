// SPDX-FileCopyrightText: 2026 GFZ Helmholtz Centre for Geosciences
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Matrix.h"
#include "MatrixSolver.h"
#include "BoundaryConditionType.hpp"
#include <xtensor/xtensor.hpp>
#include <Eigen/Dense>
#include "ppfv/Mesh.h"
#include "ppfv/Cases/VERB_ppfv2d.h"
#include "ppfv/Solver.h"

Eigen::VectorXd convertToEigen(const Matrix1D<double> &mat);
xt::xtensor<double, 2> convertToXt(const Matrix2D<double> &mat);

bool Diffusion_2D_PPFV(
		Matrix2D<double> &psd,
		const Matrix1D<double>& x, const Matrix1D<double>& y,
        int &x_size, int &y_size,
		const Matrix1D<double>& x_LBC, const Matrix1D<double>& x_UBC,
		const Matrix1D<double>& y_LBC, const Matrix1D<double>& y_UBC,
		BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
		BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
		const Matrix2D<double>& Dxx, const Matrix2D<double>& Dyy, const Matrix2D<double>& Dxy,
		const Matrix2D<double>& G, const Matrix2D<double>& Sources, const Matrix2D<double>& Losses, double total_time ,double dt
	);