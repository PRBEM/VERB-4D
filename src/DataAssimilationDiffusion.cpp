// SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file DataAssimilationDiffusion.cpp
**/

#include "DataAssimilationDiffusion.h"

#include <omp.h>

#include <iomanip>
#include <sstream>

#include "DataAssimilationHelper.h"
#include "MatrixOperations.h"


namespace data_assimilation {

Diffusion1DAnalysisCovariances::Diffusion1DAnalysisCovariances(
    size_t R_size, size_t V_size, size_t K_size)
    : _data(V_size, std::vector<Matrix2D<double>>(K_size)) {
    
    for (size_t i = 0; i < V_size; ++i) {
        for (size_t j = 0; j < K_size; ++j) {
            _data[i][j].AllocateMemory(R_size, R_size);
            _data[i][j] = 0.;
        }
    }
}

const std::vector<Matrix2D<double>>& Diffusion1DAnalysisCovariances::operator[](size_t iV) const {
    return _data[iV];
}

std::vector<Matrix2D<double>>& Diffusion1DAnalysisCovariances::operator[](size_t iV) {
    return _data[iV];
}

DataAssimilationManagerDiffusion::DataAssimilationManagerDiffusion(
    const std::string &parametersFile,
    const Matrix1D<double> &R,
    const Matrix2D<double> &V,
    const Matrix2D<double> &K,
    const std::string& debug_output_folder)
    : _analysisCovariance(R.size_q1, V.size_q1, K.size_q2), _R(R), _V(V), _K(K) {

    _assimilationParameters = readParameters(parametersFile);
    _runDataAssimilation = _assimilationParameters.runDataAssimilation;
    _debug_output_folder = debug_output_folder;
    
    if (_runDataAssimilation) {
        if (_assimilationParameters.dataSource == DataAssimilationDataSource::DataServer) {
            throw std::runtime_error("Deprecated option!");
        } else if (_assimilationParameters.dataSource == DataAssimilationDataSource::LocalFiles) {
            _dataSource = std::make_unique<LocalFilesDataSource>(LocalFilesDataSource("satellite_data.lst", R.size_q1, V.size_q1, K.size_q2));
        } else {
            throw std::invalid_argument("Encountered invalid 'data_source' argument in parameters_da.ini!");
        }

        // Initialize model matrix storage
        _model_matrix.resize(V.size_q1);
        for(size_t i = 0; i < V.size_q1; ++i) {
            _model_matrix[i].resize(K.size_q2);
        }
    }
}

void DataAssimilationManagerDiffusion::assimilate(
    double start_time,
    double end_time,
    Matrix4D<double> &PSD,
    const Matrix4D<double>& DLL,
    const Matrix4D<double>& G,
    Matrix3D<double>& x_LBC, Matrix3D<double>& x_UBC,
    BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
    double dt) {

    if (_runDataAssimilation) {
        
        size_t R_size = PSD.size_x;
        size_t V_size = _V.size_q1;
        size_t K_size = _K.size_q2;
        size_t total_iterations = V_size * K_size;

        Matrix4D<double> observations = _dataSource->getObservations(start_time, end_time, _R, _V, _K);
        

#ifdef DATA_ASSIMILATION_DEBUG
        data_assimilation::DebugOutput3D debug_output_3d(PSD.size_x, PSD.size_y, PSD.size_z);
#endif

        // Use progress count only for logging, shared by all threads
        int progress_count = 0;
        int progress_total = total_iterations;

        // Use a critical section for printing progress to avoid garbling, or use OMP_CRITICAL
        #pragma omp parallel for schedule(dynamic, 1) collapse(2)
        for (int iV = _V.size_q1 - 1; iV >= 0; --iV) {
            for (size_t iK = 0; iK < _K.size_q2; ++iK) {
                

                #pragma omp atomic
                progress_count++;
                
                if (Logger::getDebugLevel() == Logger::DebugLevel::DEBUG_LEVEL_DEBUG) {
                    if (omp_get_thread_num() == 0) {
                        std::cout << "\b\b\b\b\b\b\b\b\b" << '\t'
                                          << std::fixed << std::setprecision(1) 
                                          << (double)progress_count / progress_total * 100.0 << "%" << std::flush;
                    }
                }

                Matrix1D<double> PSD_R = PSD.wyzSlice(0, iV, iK);

                data_assimilation::DebugOutput2DDiffusion debug_output_2d = runKalmanFilterDiffusion1D(
                    PSD_R,
                    _analysisCovariance[iV][iK],
                    _model_matrix[iV][iK],
                    _R,
                    DLL.wyzSlice(0, iV, iK),
                    G.wyzSlice(0, iV, iK),
                    x_LBC[0][iV][iK],
                    x_UBC[0][iV][iK],
                    x_LBC_type,
                    x_UBC_type,
                    observations.wyzSlice(0, iV, iK), 
                    dt,
                    _assimilationParameters);

#ifdef DATA_ASSIMILATION_DEBUG
                debug_output_3d.insert_output_2d(debug_output_2d, iV, iK);
#endif

                // Copy the updated forecast back to the main PSD array
                for (size_t iR = 0; iR < R_size; iR++) {
                    PSD[0][iR][iV][iK] = PSD_R[iR];
                }
            }
        }

#ifdef DATA_ASSIMILATION_DEBUG
        debug_output_3d.write_files(start_time, _debug_output_folder);
#endif

    }
}

DebugOutput2DDiffusion runKalmanFilterDiffusion1D(
    Matrix1D<double> &forecast,
    Matrix2D<double> &analysisErrorCovariance,
    Matrix2D<double> &model_matrix,
    const Matrix1D<double>& R,
    const Matrix1D<double>& DLL,
    const Matrix1D<double>& G,
    double x_LBC, double x_UBC,
    BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
    const Matrix1D<double> &observations,
    const double timeStep,
    const Parameters &parameters) {

    assert(forecast[forecast.size_q1 - 1] == x_UBC);
    assert(forecast[0] == x_LBC);

    data_assimilation::DebugOutput2DDiffusion debug_output;

    // --- 1. Compute Model Matrix (State Transition Matrix F) ---
    internal::getModelMatrixDiffusion1D(
        model_matrix, R, DLL, G, x_LBC, x_UBC, x_LBC_type, x_UBC_type, timeStep);
    
    Matrix2D<double> model_matrix_reduced = Matrix2D<double>(model_matrix.size_q1-2, model_matrix.size_q2-2);
    for (size_t i = 1; i < model_matrix.size_q1 - 1; ++i) {
        for (size_t j = 1; j < model_matrix.size_q2 -1; ++j) {
            model_matrix_reduced[i-1][j-1] = model_matrix[i][j];
        }
    }


    Matrix1D<double> forecast_reduced = Matrix1D<double>(forecast.size_q1 - 2);
    for (size_t i = 1; i < forecast.size_q1 - 1; ++i) {
        forecast_reduced[i-1] = forecast[i];
    }


#ifdef DATA_ASSIMILATION_DEBUG
    debug_output.observations = observations;
    debug_output.forecast_init = Matrix1D<double>(forecast);
    debug_output.model = Matrix2D<double>(model_matrix);
#endif

    // --- 2. Convert to Observation Space ---
    auto observationSpace = internal::convertToObservationSpace(observations);

    // --- 3. Setup Noise Covariance Matrices (Q and Robs) ---
    double modelError = parameters.modelError;
    double observationError = parameters.observationError;

    // Process Noise Covariance (Q)
    Matrix1D<double> D = Matrix1D<double>{forecast.size_q1};
    D = modelError;

    auto Q = diag(D);

    // Measurement Noise Covariance (Robs)
    Matrix1D<double> E = Matrix1D<double>(observationSpace.data.size_q1);
    E = observationError;

    auto Robs = diag(E);

    // --- 4. Initialize Analysis Error Covariance (P_a) on first step ---
    if (analysisErrorCovariance.maxabs() == 0) {
        for (size_t iP = 0; iP < analysisErrorCovariance.size_q1; iP++) {
            // Set diagonal elements to an initial high uncertainty value (e.g., 5*modelError)
            analysisErrorCovariance[iP][iP] = 5 * modelError;
        }
    }

    Matrix2D<double> H_reduced = Matrix2D<double>(observationSpace.H.size_q1, observationSpace.H.size_q2 - 2);
    for (size_t i = 0; i < observationSpace.H.size_q1; ++i) {
        for (size_t j = 1; j < observationSpace.H.size_q2 -1; ++j) {
            H_reduced[i][j-1] = observationSpace.H[i][j];
        }
    }

    // --- 5. Run Kalman Filter ---
    Matrix2D<double> Pf = runKalmanFilter( // Pf is the a-priori (forecast) error covariance
        forecast,
        analysisErrorCovariance, // This matrix will be overwritten with the a-posteriori error covariance P_a
        model_matrix,
        Q,
        observationSpace.data,
        observationSpace.H,
        Robs);

    assert(forecast[forecast.size_q1 - 1] == x_UBC);
    assert(forecast[0] == x_LBC);

    // Apply Dirichlet Boundary Conditions to error covariance matrices
    // Pf[0][0] = observationError;
    // Pf[Pf.size_q1-1][Pf.size_q1-1] = observationError;
    // analysisErrorCovariance[0][0] = observationError;
    // analysisErrorCovariance[analysisErrorCovariance.size_q1-1][analysisErrorCovariance.size_q1-1] = observationError;

#ifdef DATA_ASSIMILATION_DEBUG
    debug_output.forecast_result = Matrix1D<double>(forecast);
    debug_output.a_priori_err_cov = Pf;
    debug_output.a_posteriori_err_cov = Matrix2D<double>(analysisErrorCovariance);
#endif

    return debug_output;
}

} // namespace data_assimilation