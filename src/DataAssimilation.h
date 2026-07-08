// SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file DataAssimilation.h
**/

#pragma once

#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <filesystem>

#include "Matrix.h"
#include "UpdatableMatrix.h"
#include "Parameters.h"

using ParametersIni = ::Parameters;
namespace data_assimilation {
struct Parameters {
    bool runDataAssimilation;
    bool useLog;
    double timeStep;
    double modelError;
    double observationError;
    double correlationTime = std::numeric_limits<double>::infinity();
    DataAssimilationDataSource dataSource = DataAssimilationDataSource::DataServer;
    bool growing_Q_at_boundary = false;
};
struct Observations {
    Matrix1D<double> P;
    Matrix1D<double> R;
    Matrix1D<double> PSD;
};
struct ObservationSpace {
    Matrix1D<double> data;
    Matrix2D<double> H;
};

class DataSource {
public:

    virtual Matrix4D<double> getObservations(
        const double timeStart, const double timeEnd,
        const Matrix2D<double>& P, const Matrix2D<double>& R,
        const Matrix2D<double>& V, const Matrix2D<double>& K) = 0;

    virtual Matrix4D<double> getObservations(
        const double timeStart, const double timeEnd,
        const Matrix1D<double>& R,
        const Matrix2D<double>& V,
        const Matrix2D<double>& K) = 0;
    };

class LocalFilesDataSource : public DataSource {
public:
    LocalFilesDataSource(const std::string& satellite_lst_file, 
            int size_q1, int size_q2, int size_q3, int size_q4);

    LocalFilesDataSource(const std::string& satellite_lst_file, 
            int size_q1, int size_q2, int size_q3);

    virtual Matrix4D<double> getObservations(
        const double timeStart, const double timeEnd,
        const Matrix2D<double>& P, const Matrix2D<double>& R,
        const Matrix2D<double>& V, const Matrix2D<double>& K) override;

    virtual Matrix4D<double> getObservations(
        const double timeStart, const double timeEnd,
        const Matrix1D<double>& R,
        const Matrix2D<double>& V,
        const Matrix2D<double>& K) override;

private:
    UpdatableListMatrix<Matrix4D<double>> _data;
};

struct ProcessedMatFileData {
    Matrix1D<double> MLT;
    Matrix1D<double> R;
    Matrix3D<double> V;
    Matrix2D<double> K;
    Matrix3D<double> PSD;
};

Parameters readParameters(const std::string& filename);

struct DebugOuput2D {
    Matrix2D<double> observations;
    Matrix2D<double> forecast_init;
    Matrix2D<double> forecast_result;
    Matrix2D<double> a_priori_err_cov;
    Matrix2D<double> a_posteriori_err_cov;
    Matrix2D<double> model;
};

Matrix2D<double> runKalmanFilter(
    Matrix1D<double>& forecast,
    Matrix2D<double>& Pa,
    const Matrix2D<double>& M, const Matrix2D<double>& Q,
    const Matrix1D<double>& obs, const Matrix2D<double>& H,
    const Matrix2D<double>& R);

}  // namespace data_assimilation
