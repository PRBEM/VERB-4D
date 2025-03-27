// SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
//
// SPDX-License-Identifier: BSD-3-Clause

#include "DataAssimilation.h"

#include <omp.h>

#include <iomanip>

#include "DataAssimilationHelper.h"
#include "MatrixOperations.h"

// namespace da = data_assimilation;

data_assimilation::DataAssimilationManagerConvection::DataAssimilationManagerConvection(
    const std::string &parametersFile,
    double timeStart, double timeEnd, const Matrix2D<double> &P, 
    const Matrix2D<double> &R, const Matrix2D<double> &V, const Matrix2D<double> &K, const std::string& debug_output_folder)
    : _analysisCovarianceConvection(P.size_q1, R.size_q2, V.size_q1, K.size_q2), _timeStart(timeStart), _timeEnd(timeEnd), _P(P), _R(R), _V(V), _K(K) {

    _assimilationParameters = readParameters(parametersFile);
    _runDataAssimilation = _assimilationParameters.runDataAssimilation;
    _debug_output_folder = debug_output_folder;

    if (_runDataAssimilation) {

        if (_assimilationParameters.dataSource == DataAssimilationDataSource::DataServer) {
            #if (MATLAB_CAPABLE)
            _dataSource = std::make_unique<DataServerDataSource>(DataServerDataSource("satellite_data.lst"));
            #else
            throw std::invalid_argument("'data_source' was set to 'DataSever', but the solver is compiled without Matlab capabilities!");
            #endif
        } else if (_assimilationParameters.dataSource == DataAssimilationDataSource::LocalFiles) {
            _dataSource = std::make_unique<LocalFilesDataSource>(LocalFilesDataSource("satellite_data.lst", P.size_q1, R.size_q2, V.size_q1, K.size_q2));
        } else {
            throw std::invalid_argument("Encountered invalid 'data_source' argument in parameters_da.ini!");
        }
        _timePrev = timeStart;
        _timeNext = timeStart + _assimilationParameters.timeStep;

        _model_matrix = std::vector(V.size_q1, std::vector<Matrix2D<double>>(K.size_q2));
    }
}

void data_assimilation::DataAssimilationManagerConvection::assimilate(
    double time, Matrix4D<double> &PSD,
    const Matrix4D<double> &VP, const Matrix4D<double> &VR,
    const bool has_VX_updated,
    const Matrix4D<double> &Loss,
    const Matrix4D<double> &SaturationDensity,
    double dt) {

    if (_runDataAssimilation && time >= _timeNext) {

        size_t P_size = PSD.size_w;
        size_t R_size = PSD.size_x;

        Matrix4D<double> observations = _dataSource->getObservations(_timePrev, _timeNext, _P, _R, _V, _K);

#ifdef DATA_ASSIMILATION_DEBUG
        data_assimilation::DebugOuput4D debug_output_4d(VP.size_w, VP.size_x, VP.size_y, VP.size_z);
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

                data_assimilation::DebugOuput2D debug_output_2d = runKalmanFilterConvection2D(PSD_PR,
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
        debug_output_4d.write_files(_timePrev, _debug_output_folder);
#endif

        Logger::debug << '\n';
        _timePrev = _timeNext;
        _timeNext += _assimilationParameters.timeStep;

    }
}

// std::vector<std::vector<data_assimilation::Observations>> data_assimilation::getObservations(
//     double timeStart,
//     double timeEnd,
//     const Matrix2D<double>& V,
//     const Matrix2D<double>& K,
//     const std::vector<pmf::Parameters>& parameters) {

//     std::vector<ProcessedMatFileData> pmfDataSplit;
//     for (auto par : parameters) {
//         ProcessedMatFileData one_instrument = internal::readData(timeStart, timeEnd, par);
//         if (one_instrument.MLT.initialized) {
//             // if MLT was initalized, all of the other variables are initialized
//             pmfDataSplit.push_back(one_instrument);
//         }
//     }
//     std::vector<std::vector<data_assimilation::Observations>> result;
//     if (pmfDataSplit.empty()) {
//         Logger::debug << "\tNo data found...\n";
//     } else {
//         Logger::debug << "\tDone reading Data. Interpolating...\n";
//         result = internal::interpolate(pmfDataSplit, V, K);
//     }
//     return result;
// }

data_assimilation::DebugOuput2D data_assimilation::runKalmanFilterConvection2D(
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

    data_assimilation::DebugOuput2D debug_output;

    auto forecast1D = toMatrix1D(forecast);
    if (has_VX_updated || not model_matrix.initialized) {
        model_matrix = internal::getModelMatrixConvection2D(
            VP, VR, Loss, saturation_map, timeStep, dP, dR, parameters.useLog, parameters.correlationTime);
    }

#ifdef DATA_ASSIMILATION_DEBUG
    debug_output.observations = observations;
    debug_output.forecast_init = Matrix2D<double>(forecast);
    debug_output.error_covariance = Matrix2D<double>(analysisErrorCovariance);
    debug_output.model = Matrix2D<double>(model_matrix);
#endif

    auto observationSpace = internal::convertToObservationSpace(observations);

    double modelError = parameters.modelError;
    double observationError = parameters.observationError;
    if (parameters.useLog) {
        forecast1D = log10(forecast1D);
        observationSpace.data = log10(observationSpace.data);

        modelError = log10(1.0 + modelError);
        observationError = log10(1.0 + observationError);
    }
    Matrix1D<double> D = Matrix1D<double>{forecast1D.size_q1};
    D = modelError;
    auto Q = parameters.useLog ? diag(D) : diag(forecast1D * modelError);

    if (not parameters.growing_Q_at_boundary) {
        double correlation_factor = exp(-timeStep/parameters.correlationTime);
        for (size_t iP = 0; iP < P.size_q1; iP++) {
            Q[P.size_q2-1 + iP * P.size_q2][P.size_q2-1 + iP * P.size_q2] = (1.-(correlation_factor*correlation_factor)) * modelError;
            Q[iP * P.size_q2][iP * P.size_q2] = (1.-(correlation_factor*correlation_factor)) * modelError;
        }
    }

    Matrix1D<double> E = Matrix1D<double>(observationSpace.data.size_q1);
    E = observationError;
    auto Robs = parameters.useLog ? diag(E) : diag(observationSpace.data * observationError);

    runKalmanFilter(
        forecast1D,
        analysisErrorCovariance,
        model_matrix,
        Q,
        observationSpace.data,
        observationSpace.H,
        Robs);

    if (parameters.useLog) {
        forecast1D = pow(10., forecast1D);
    }

    forecast = toMatrix2D(forecast1D, P.size_q1, R.size_q2);

#ifdef DATA_ASSIMILATION_DEBUG
    debug_output.forecast_result = Matrix2D<double>(forecast);
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

void data_assimilation::runKalmanFilter(
    Matrix1D<double>& forecast,
    Matrix2D<double>& Pa,
    const Matrix2D<double>& M, const Matrix2D<double>& Q,
    const Matrix1D<double>& obs, const Matrix2D<double>& H,
    const Matrix2D<double>& R) {

    // forecast error covariance matrix [n x n]
    Matrix2D<double> Pf = abtrans(M * Pa, M);
    Pf += Q;

    Matrix2D<double> Pf_times_HT = abtrans(Pf, H);

    // Kalman matrix [n x m], Pf_times_HT will be overwritten and referenced now as K
    Matrix2D<double>& K = trans_solve(H * Pf_times_HT + R, Pf_times_HT);

    forecast += K * (obs - H * forecast);

    Pa = Pf - K * H * Pf;
}

data_assimilation::Parameters data_assimilation::readParameters(const std::string& filename) {
    data_assimilation::Parameters result;

    if (std::filesystem::exists(filename)) {

        ParametersIni parameters(filename);

        Logger::message << std::endl;
        Logger::writeSeparator();
        Logger::message << std::setw(58) << "Data assimilation" << std::endl;
        Logger::writeSeparator();

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

#if (MATLAB_CAPABLE)

data_assimilation::DataServerDataSource::DataServerDataSource(const std::string& satellite_lst_file) {
    _dataParameters = pmf::readParameters(satellite_lst_file);
};

Matrix4D<double> data_assimilation::DataServerDataSource::getObservations(
        const double timeStart, const double timeEnd,
        const Matrix2D<double>& P, const Matrix2D<double>& R,
        const Matrix2D<double>& V, const Matrix2D<double>& K) {

    std::vector<ProcessedMatFileData> pmfDataSplit;

    for (auto par : _dataParameters) {
        pmfDataSplit.push_back(internal::readData(timeStart, timeEnd, par));
    }
    Logger::debug << "\tDone reading Data. Interpolating and binning...\n";
    std::vector<std::vector<data_assimilation::Observations>> observations_interpolated_V_K = internal::interpolate(pmfDataSplit, V, K);

    Matrix4D<double> result = internal::bin(observations_interpolated_V_K, P, R, "log10");

    return result;
};

#endif

data_assimilation::LocalFilesDataSource::LocalFilesDataSource(const std::string& satellite_lst_file,
        int size_q1, int size_q2, int size_q3, int size_q4) {

    _data = UpdatableListMatrix<Matrix4D<double>>(UpdatableListMatrix<Matrix4D<double>>::MERGE_TYPE::MEAN);
    _data.AllocateMemory(size_q1, size_q2, size_q3, size_q4);
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
