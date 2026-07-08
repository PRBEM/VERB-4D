// SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file DataAssimilation.cpp
**/

#include "DataAssimilation.h"

#include <omp.h>

#include <iomanip>

#include "DataAssimilationHelper.h"
#include "MatrixOperations.h"

Matrix2D<double> data_assimilation::runKalmanFilter(
    Matrix1D<double>& forecast,
    Matrix2D<double>& Pa,
    const Matrix2D<double>& M, const Matrix2D<double>& Q,
    const Matrix1D<double>& obs, const Matrix2D<double>& H,
    const Matrix2D<double>& R) {

    // --- 1. A Priori (Forecast) Error Covariance (Pf) ---
    // Pf = M * Pa * M^T + Q
    
    // Pf = M * Pa * M^T
    Matrix2D<double> Pf = abtrans(M * Pa, M); // M * Pa * M^T
    
    // Pf += Q (Add Process Noise Covariance)
    Pf += Q;
    

    // --- 2. Kalman Gain (K) ---
    // K = Pf * H^T * (H * Pf * H^T + R)^-1

    // Pf_times_HT = Pf * H^T
    Matrix2D<double> Pf_times_HT = abtrans(Pf, H);

    // Innovation Covariance = H * Pf * H^T + R
    Matrix2D<double> InnovationCovariance = H * Pf_times_HT;
    InnovationCovariance += R;
    
    // K = trans_solve(InnovationCovariance, Pf_times_HT)
    // Pf_times_HT is used as the output variable K, and is overwritten.
    Matrix2D<double>& K = trans_solve(InnovationCovariance, Pf_times_HT); 
    

    // --- 3. Analysis (Updated) Forecast ---
    // forecast_a = forecast_f + K * (obs - H * forecast_f)
    
    // Innovation (y - H*x_f): obs - H * forecast
    Matrix1D<double> innovation = obs - (H * forecast);
    
    // forecast += K * innovation (Update)
    Matrix1D<double> correction = K * innovation;

    Matrix1D<double> forecast_cache = forecast;
    forecast += correction;

    for (size_t i = 0; i < forecast.size_q1; ++i) {
        if (std::isnan(forecast[i])) {
            throw std::runtime_error("NaN detected in forecast vector after Kalman update.");
        }

    }

    // --- 4. A Posteriori (Analysis) Error Covariance (Pa) ---
    // Pa = (I - K * H) * Pf
    
    // Pa = Pf - K * H * Pf (equivalent to (I - K*H) * Pf)
    Pa = Pf - (K * H * Pf);
    
    
    return Pf;
}

data_assimilation::Parameters data_assimilation::readParameters(const std::string& filename) {
    data_assimilation::Parameters result;

    if (std::filesystem::exists(filename)) {

        ParametersIni parameters(filename);

        parameters.getParameter("run_data_assimilation", result.runDataAssimilation, true);

        if (result.runDataAssimilation) {
            parameters.getParameter("useLog", result.useLog, true);
            parameters.getParameter("time_step", result.timeStep, true);
            parameters.getParameter("model_error", result.modelError, true);
            parameters.getParameter("observation_error", result.observationError, true);
            parameters.getParameter("data_source", result.dataSource, false);
            parameters.getParameter("correlation_time", result.correlationTime, false);
            parameters.getParameter("growing_Q_at_boundary", result.growing_Q_at_boundary, false);

            double correlation_factor = exp(-result.timeStep/result.correlationTime);
            Logger::debug << "\tData assimilation time step: " << result.timeStep << std::endl;
            Logger::debug << "\tCorrelation time: " << result.correlationTime << std::endl;
            Logger::debug << "\tCorrelation factor: " << correlation_factor << std::endl;
        }
    }

    return result;
}

data_assimilation::LocalFilesDataSource::LocalFilesDataSource(const std::string& satellite_lst_file,
        int size_q1, int size_q2, int size_q3, int size_q4) {

    _data = UpdatableListMatrix<Matrix4D<double>>(UpdatableListMatrix<Matrix4D<double>>::MERGE_TYPE::MEAN);
    _data.AllocateMemory(size_q1, size_q2, size_q3, size_q4);
    _data.readFromIniFile(satellite_lst_file, _data, _data, _data, _data);

};

data_assimilation::LocalFilesDataSource::LocalFilesDataSource(const std::string& satellite_lst_file,
        int size_q1, int size_q2, int size_q3) {

    _data = UpdatableListMatrix<Matrix4D<double>>(UpdatableListMatrix<Matrix4D<double>>::MERGE_TYPE::MEAN);
    _data.AllocateMemory(1, size_q1, size_q2, size_q3);
    _data.readFromIniFile(satellite_lst_file, _data, _data, _data, _data);

};

Matrix4D<double> data_assimilation::LocalFilesDataSource::getObservations(
        const double timeStart, [[maybe_unused]] const double timeEnd,
        [[maybe_unused]] const Matrix2D<double>& P, [[maybe_unused]] const Matrix2D<double>& R,
        [[maybe_unused]] const Matrix2D<double>& V, [[maybe_unused]] const Matrix2D<double>& K) {

    _data.update(timeStart, _data, _data, _data, _data);

    for (size_t w = 0; w < _data.size_w; w++) {
        for (size_t x = 0; x < _data.size_x; x++) {
            for (size_t y = 0; y < _data.size_y; y++) {
                for (size_t z = 0; z < _data.size_z; z++) {
                    if (_data[w][x][y][z] == 0) {
                        _data[w][x][y][z] = NAN;
                    }
                }
            }
        }
    }

    return _data;
};

Matrix4D<double> data_assimilation::LocalFilesDataSource::getObservations(
        const double timeStart,
        [[maybe_unused]] const double timeEnd,
        [[maybe_unused]] const Matrix1D<double>& R,
        [[maybe_unused]] const Matrix2D<double>& V,
        [[maybe_unused]] const Matrix2D<double>& K) {

    _data.update(timeStart, _data, _data, _data, _data);

    for (size_t x = 0; x < _data.size_x; x++) {
        for (size_t y = 0; y < _data.size_y; y++) {
            for (size_t z = 0; z < _data.size_z; z++) {
                if (_data[0][x][y][z] == 0) {
                    _data[0][x][y][z] = NAN;
                }
            }
        }
    }

    return _data;
};