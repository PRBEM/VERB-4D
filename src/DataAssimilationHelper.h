// SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file DataAssimilationHelper.h
**/

#pragma once
#include "DataAssimilation.h"
#include "Diffusion_1D.h"
#include "Matrix.h"
#include "Parameters.h"
#include "BoundaryConditionType.hpp"
#include <string>
#include <vector>

extern const double FILLVAL;
using ParametersIni=::Parameters;
namespace data_assimilation::internal {
    void getModelMatrixDiffusion1D(
        Matrix2D<double> &F,
        const Matrix1D<double>& R,
        const Matrix1D<double>& DLL,
        const Matrix1D<double>& G,
        const double x_LBC,
        const double x_UBC,
        const BoundaryConditionType x_LBC_type,
        const BoundaryConditionType x_UBC_type,
        const double timeStep);
    Matrix2D<double> getModelMatrixConvection2D(
        const Matrix2D<double> &VP,
        const Matrix2D<double> &VR,
        const Matrix2D<double> &Loss,
        const Matrix2D<bool> &saturation_map,
        double timeStep,
        double dP,
        double dR,
        bool use_log,
        double correlation_time);
    Matrix2D<double> getModelMatrixConvection2DNoStabilityCheck(
        const Matrix2D<double> &VP,
        const Matrix2D<double> &VR,
        const Matrix2D<double> &Loss,
        const Matrix2D<bool> &saturation_map,
        double timeStep,
        double dP,
        double dR,
        bool use_log,
        double correlation_factor);
    std::pair<int, double> splitTimeStepCourantCondition(
        double maximumCourantNumber, double timeStep, 
        const Matrix2D<double>& VP, const Matrix2D<double>& VR, 
        double dP, double dR
    );
    Matrix2D<double> getModelMatrixConvection2D(
        const Matrix2D<double>& VP,
        const Matrix2D<double>& VR,
        const Matrix2D<double>& Loss,
        double timeStep,
        double dP,
        double dR
    );

    double conditionalMean(const Matrix1D<double>& a, const Matrix1D<bool>& flag);

    data_assimilation::ObservationSpace convertToObservationSpace(const Matrix2D<double>& data);
    data_assimilation::ObservationSpace convertToObservationSpace(const Matrix1D<double>& data);

}
