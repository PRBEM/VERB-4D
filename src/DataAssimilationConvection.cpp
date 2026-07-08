// SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file DataAssimilationConvection.cpp
**/

#include "DataAssimilationConvection.h"

#include <omp.h>

#include <iomanip>

#include "DataAssimilationHelper.h"
#include "MatrixOperations.h"

data_assimilation::DataAssimilationManagerConvection::DataAssimilationManagerConvection(
    const std::string &parametersFile,
    const Matrix2D<double> &P, 
    const Matrix2D<double> &R, const Matrix2D<double> &V, const Matrix2D<double> &K, const std::string& debug_output_folder)
    : _analysisCovarianceConvection(P.size_q1, R.size_q2, V.size_q1, K.size_q2), _P(P), _R(R), _V(V), _K(K) {

    _assimilationParameters = readParameters(parametersFile);
    _runDataAssimilation = _assimilationParameters.runDataAssimilation;
    _debug_output_folder = debug_output_folder;

    if (_runDataAssimilation) {

        if (_assimilationParameters.dataSource == DataAssimilationDataSource::DataServer) {
            throw std::runtime_error("Deprecated option!");
        } else if (_assimilationParameters.dataSource == DataAssimilationDataSource::LocalFiles) {
            _dataSource = std::make_unique<LocalFilesDataSource>(LocalFilesDataSource("satellite_data.lst", P.size_q1, R.size_q2, V.size_q1, K.size_q2));
        } else {
            throw std::invalid_argument("Encountered invalid 'data_source' argument in parameters_da.ini!");
        }

        _model_matrix = std::vector(V.size_q1, std::vector<Matrix2D<double>>(K.size_q2));
    }
}

void data_assimilation::DataAssimilationManagerConvection::assimilate(
    double start_time, double end_time, Matrix4D<double> &PSD,
    const Matrix4D<double> &VP, const Matrix4D<double> &VR,
    const bool has_VX_updated,
    const Matrix4D<double> &Loss,
    const Matrix4D<double> &SaturationDensity,
    double dt) {

    if (_runDataAssimilation) {

        size_t P_size = PSD.size_w;
        size_t R_size = PSD.size_x;

        Matrix4D<double> observations = _dataSource->getObservations(start_time, end_time, _P, _R, _V, _K);

#ifdef DATA_ASSIMILATION_DEBUG
        data_assimilation::DebugOutput4D debug_output_4d(VP.size_w, VP.size_x, VP.size_y, VP.size_z);
#endif

        Logger::debug << "\tDone Interpolating. Applying Kalman filter...\n";
        if (has_VX_updated) {
            Logger::debug << "\tVelocities have changed; recomputing model matrix..." << std::endl;
        }

        int progress_count = 0;
        int progress_total = _V.size_q1 * _K.size_q2;

#pragma omp parallel for schedule(dynamic, 1) collapse(2)
        for (int iV = _V.size_q1 - 1; iV >= 0; --iV) {
            for (size_t iK = 0; iK < _K.size_q2; ++iK) {

                if (Logger::getDebugLevel() == Logger::DebugLevel::DEBUG_LEVEL_DEBUG) {
                    if (omp_get_thread_num() == 0) {
                        std::cout << "\b\b\b\b\b\b\b\b\b" << '\t'
                                << (int)((double)progress_count / progress_total * 100) << "%" << std::flush;
                    }
                }

                Matrix2D<double> PSD_PR = PSD.yzSlice(iV, iK);
                Matrix2D<bool> saturation_map;
                if (SaturationDensity.initialized) {
                    saturation_map = PSD_PR >= SaturationDensity.yzSlice(iV, iK);
                } else {
                    saturation_map = Matrix2D<bool>(PSD_PR.size_q1, PSD_PR.size_q2);
                    saturation_map = false;
                }
                Matrix2D<double> loss_rate = Loss.yzSlice(iV, iK).divide(PSD_PR);

                data_assimilation::DebugOutput2D debug_output_2d = runKalmanFilterConvection2D(PSD_PR,
                                            _analysisCovarianceConvection[iV][iK],
                                            _model_matrix[iV][iK],
                                            _P, _R,
                                            VP.yzSlice(iV, iK), VR.yzSlice(iV, iK),
                                            has_VX_updated,
                                            loss_rate,
                                            saturation_map,
                                            observations.yzSlice(iV, iK), 
                                            dt, _assimilationParameters);

#ifdef DATA_ASSIMILATION_DEBUG
            debug_output_4d.insert_output_2d(debug_output_2d, iV, iK);
#endif

                for (size_t iP = 0; iP < P_size; iP++) {
                    for (size_t iR = 0; iR < R_size; iR++) {
                        PSD[iP][iR][iV][iK] = PSD_PR[iP][iR];
                    }
                }
                
                progress_count++;
            }
        }

#ifdef DATA_ASSIMILATION_DEBUG
        Logger::debug << "Writing DA Debug output..." << std::endl;
        debug_output_4d.write_files(start_time, _debug_output_folder);
#endif

        Logger::debug << '\n';
    }
}

data_assimilation::DebugOutput2D data_assimilation::runKalmanFilterConvection2D(
    Matrix2D<double> &forecast,
    Matrix2D<double> &analysisErrorCovariance,
    Matrix2D<double> &model_matrix,
    const Matrix2D<double> &P,
    const Matrix2D<double> &R,
    const Matrix2D<double> &VP,
    const Matrix2D<double> &VR,
    const bool has_VX_updated,
    const Matrix2D<double> &Loss,
    const Matrix2D<bool> &saturation_map,
    const Matrix2D<double> &observations,
    const double timeStep,
    const Parameters &parameters) {

    auto dP = P[1][0] - P[0][0];
    auto dR = R[0][1] - R[0][0];

    data_assimilation::DebugOutput2D debug_output;

    auto forecast1D = toMatrix1D(forecast);
    if (has_VX_updated || not model_matrix.initialized) {
        model_matrix = internal::getModelMatrixConvection2D(
            VP, VR, Loss, saturation_map, timeStep, dP, dR, parameters.useLog, parameters.correlationTime);
    }

#ifdef DATA_ASSIMILATION_DEBUG
    debug_output.observations = observations;
    debug_output.forecast_init = Matrix2D<double>(forecast);
    debug_output.model = Matrix2D<double>(model_matrix);
#endif

    auto observationSpace = internal::convertToObservationSpace(observations);

    double modelError = parameters.modelError;
    double observationError = parameters.observationError;
    if (parameters.useLog) {
        forecast1D = log(forecast1D);
        observationSpace.data = log(observationSpace.data);

        // modelError = log(1.0 + modelError);
        // observationError = log(1.0 + observationError);
    }
    Matrix1D<double> D = Matrix1D<double>{forecast1D.size_q1};
    D = modelError;
    // auto Q = parameters.useLog ? diag(D) : diag(forecast1D * modelError);

    double correlation_factor = exp(-timeStep/parameters.correlationTime);
    auto tmp = diag(D);
    for (size_t iP = 0; iP < P.size_q1; iP++) {
        tmp[P.size_q2-1 + iP * P.size_q2][P.size_q2-1 + iP * P.size_q2] = (1.-(correlation_factor*correlation_factor)) * modelError;
        tmp[iP * P.size_q2][iP * P.size_q2] = (1.-(correlation_factor*correlation_factor)) * modelError;
    }

    auto model_matrix_1day = model_matrix^(std::ceil(3./24 / timeStep));

    // sum Q elements
    // float total_sum = 0;
    // for (size_t i = 0; i < tmp.size_q1; i++) {
    //     for (size_t j = 0; j < tmp.size_q1; j++) {
    //         total_sum += tmp[i][j];
    //     }
    // }
    // std::cout << "total sum before: " << total_sum << std::endl;

    auto Q = model_matrix_1day*tmp*transpose(model_matrix_1day) + diag(D)/10;

    // total_sum = 0;
    // for (size_t i = 0; i < Q.size_q1; i++) {
    //     for (size_t j = 0; j < Q.size_q1; j++) {
    //         total_sum += Q[i][j];
    //     }
    // }
    // std::cout << "total sum after: " << total_sum << std::endl;

    // std::cout << "co-variances: " << Q.maxabs() << std::endl;
    // for (size_t i = 0; i < Q.size_q1; i++) {
    //     std::cout << (i - P.size_q2+1) % P.size_q2 << std::endl;
    //     if ((i - P.size_q2+1) % P.size_q2 == 0 || i % P.size_q2 == 0) {
    //         continue; // skip boundary points
    //     }
    //     std::cout << Q[i][i] << std::endl;
    //     double scale_factor = modelError / Q[i][i];
    //     for (size_t j = 0; j < Q.size_q1; j++) {
    //         Q[i][j] *= scale_factor;
    //     }
    //     assert(Q[i][i] - modelError < 1e-16);
    // }
    // std::cout << "with variances: " << Q.maxabs() << std::endl;

    if (not parameters.growing_Q_at_boundary) {
        for (size_t iP = 0; iP < P.size_q1; iP++) {
            Q[P.size_q2-1 + iP * P.size_q2][P.size_q2-1 + iP * P.size_q2] = (1.-(correlation_factor*correlation_factor)) * modelError;
            Q[iP * P.size_q2][iP * P.size_q2] = (1.-(correlation_factor*correlation_factor)) * modelError;
        }
    }

    Matrix1D<double> E = Matrix1D<double>(observationSpace.data.size_q1);
    E = observationError;
    auto Robs = parameters.useLog ? diag(E) : diag(observationSpace.data * observationError);

    // first time step
    if (analysisErrorCovariance.maxabs() == 0) {
        for (size_t iP = 0; iP < analysisErrorCovariance.size_q1; iP++) {
            analysisErrorCovariance[iP][iP] = 5*modelError;
        }
    }

    Matrix2D<double> Pf = runKalmanFilter(
        forecast1D,
        analysisErrorCovariance,
        model_matrix,
        Q,
        observationSpace.data,
        observationSpace.H,
        Robs);

    if (parameters.useLog) {
        forecast1D = pow(std::exp(1), forecast1D);
    }

    forecast = toMatrix2D(forecast1D, P.size_q1, R.size_q2);

#ifdef DATA_ASSIMILATION_DEBUG
    debug_output.forecast_result = Matrix2D<double>(forecast);
    debug_output.a_priori_err_cov = Pf;
    debug_output.a_posteriori_err_cov = Matrix2D<double>(analysisErrorCovariance);
#endif
    return debug_output;
}

data_assimilation::Convection2DAnalysisCovariances::Convection2DAnalysisCovariances(
    size_t P_size, size_t R_size, size_t V_size, size_t K_size)
    : _data(V_size, std::vector<Matrix2D<double>>(K_size)) {
    for (size_t i = 0; i < V_size; ++i) {
        for (size_t j = 0; j < K_size; ++j) {
            _data[i][j].AllocateMemory(P_size * R_size, P_size * R_size);
            _data[i][j] = 0.;
        }
    }
}

const std::vector<Matrix2D<double>>& data_assimilation::Convection2DAnalysisCovariances::operator[](size_t iV) const {
    return _data[iV];
}

std::vector<Matrix2D<double>>& data_assimilation::Convection2DAnalysisCovariances::operator[](size_t iV) {
    return _data[iV];
}