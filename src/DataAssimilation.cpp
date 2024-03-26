#include "DataAssimilation.h"

#include <omp.h>

#include <iomanip>

#include "DataAssimilationHelper.h"
#include "MatrixOperations.h"

// namespace da = data_assimilation;

data_assimilation::DataAssimilationManagerConvection::DataAssimilationManagerConvection(
    const std::string& parametersFile,
    double timeStart, double timeEnd, const Matrix2D<double>& V,
    const Matrix2D<double>& K, int P_size, int R_size, const std::string& debug_output_folder)
    : _timeStart(timeStart), _timeEnd(timeEnd), _V(V), _K(K),
    _analysisCovarianceConvection(P_size, R_size, V.size_q1, K.size_q2), _debug_output_folder(debug_output_folder)
{
    _assimilationParameters = readParameters(parametersFile);
    _runDataAssimilation = _assimilationParameters.runDataAssimilation;
    _debug_output_folder = debug_output_folder;

    if (_runDataAssimilation) {
        _dataParameters = pmf::readParameters("satellite_data.list");

        _timePrev = timeStart;
        _timeNext = timeStart + _assimilationParameters.timeStep;
    }
}

void data_assimilation::DataAssimilationManagerConvection::assimilate(
    double time, Matrix4D<double>& PSD,
    const Matrix4D<double>& P, const Matrix4D<double>& R,
    const Matrix4D<double>& VP, const Matrix4D<double>& VR,
    const Matrix4D<double>& Loss,
    double dt) {
    std::cout << std::setprecision(10) << "input: " << time << "   _timeNext: " << _timeNext << std::endl;
    if (_runDataAssimilation && time >= _timeNext) {
        size_t P_size = PSD.size_w;
        size_t R_size = PSD.size_x;

        auto P2D = P.yzSlice(0, 0);
        auto R2D = R.yzSlice(0, 0);

        auto observations = getObservations(_timePrev, _timeNext, _V, _K,
                                            _dataParameters);

        if (observations.empty()) {
            std::cout << "No data found. Skipping data assimilation...\n";
        } else {
            std::cout << "Done Interpolating. Applying Kalman filter...\n";
            int progress_count = 0;
            int progress_total = _V.size_q1 * _K.size_q2;

#ifdef DATA_ASSIMILATION_DEBUG
            data_assimilation::DebugOuput4D debug_output_4d(P.size_w, P.size_x, P.size_y, P.size_z);
#endif

#pragma omp parallel for shared(progress_total, progress_count) schedule(dynamic, 1) collapse(2)
            for (int iV = _V.size_q1 - 1; iV >= 0; --iV) {
                for (int iK = 0; iK < _K.size_q2; ++iK) {
                    if (omp_get_thread_num() == 0) {
                        std::cout << "\b\b\b\b\b\b\b\b\b" << '\t'
                                  << (int)((double)progress_count / progress_total * 100) << "%" << std::flush;
                    }
                    // int numnan = 0;
                    // set all observations to zero for testing
                    // for (int l = 0; l < observations[iV][iK].PSD.size_q1; l++)
                    //     if(isnan(observations[iV][iK].PSD[l])) numnan++;
                    // observations[iV][iK].PSD[l] = 1e-31;
                    // if(numnan > 0) std::cout << numnan << " of " << observations[iV][iK].PSD.size_q1 << '\n';
                    Matrix2D<double> PSD_PR = PSD.yzSlice(iV, iK);
                    auto debug_output = runKalmanFilterConvection2D(PSD_PR,
                                                                    _analysisCovarianceConvection[iV][iK],
                                                                    P2D, R2D,
                                                                    VP.yzSlice(iV, iK), VR.yzSlice(iV, iK),
                                                                    Loss.yzSlice(iV, iK),
                                                                    observations[iV][iK], dt, _assimilationParameters);

#ifdef DATA_ASSIMILATION_DEBUG
                    debug_output_4d.insert_output_2d(debug_output, iV, iK);
#endif

                    for (auto iP = 0; iP < P_size; iP++) {
                        for (auto iR = 0; iR < R_size; iR++) {
                            PSD[iP][iR][iV][iK] = PSD_PR[iP][iR];
                        }
                    }
                    progress_count++;
                }
            }

#ifdef DATA_ASSIMILATION_DEBUG
            std::cout << _debug_output_folder << std::endl;
            debug_output_4d.write_files(_timePrev, _debug_output_folder);
#endif
        }
        std::cout << '\n';
        _timePrev = _timeNext;
        _timeNext += _assimilationParameters.timeStep;
    }
}

std::vector<std::vector<data_assimilation::Observations>> data_assimilation::getObservations(
    double timeStart,
    double timeEnd,
    const Matrix2D<double>& V,
    const Matrix2D<double>& K,
    const std::vector<pmf::Parameters>& parameters) {
    std::vector<ProcessedMatFileData> pmfDataSplit;
    for (auto par : parameters) {
        ProcessedMatFileData one_instrument = internal::readData(timeStart, timeEnd, par);
        if (one_instrument.MLT.initialized) {
            // if MLT was initalized, all of the other variables are initialized
            pmfDataSplit.push_back(one_instrument);
        }
    }
    std::cout << "Done reading Data. Interpolating...\n";
    std::vector<std::vector<data_assimilation::Observations>> result;
    if (!pmfDataSplit.empty()) {
        result = internal::interpolate(pmfDataSplit, V, K);
    }
    return result;
}

data_assimilation::DebugOuput2D data_assimilation::runKalmanFilterConvection2D(
    Matrix2D<double>& forecast,
    Matrix2D<double>& analysisErrorCovariance,
    const Matrix2D<double>& P,
    const Matrix2D<double>& R,
    const Matrix2D<double>& VP,
    const Matrix2D<double>& VR,
    const Matrix2D<double>& Loss,
    const Observations& observations,
    double timeStep,
    const Parameters& parameters) {
    auto dP = P[1][0] - P[0][0];
    auto dR = R[0][1] - R[0][0];

    data_assimilation::DebugOuput2D debug_output;

    auto forecast1D = toMatrix1D(forecast);
    auto modelMatrix = internal::getModelMatrixConvection2D(
        VP, VR, Loss, timeStep, dP, dR);

    auto observationsBinned = internal::bin(observations, P, R, "log10");

#ifdef DATA_ASSIMILATION_DEBUG
    debug_output.binned_observations = observationsBinned;
    debug_output.forecast_init = Matrix2D<double>(forecast);
    debug_output.error_covariance = Matrix2D<double>(analysisErrorCovariance);
    debug_output.model = Matrix2D<double>(modelMatrix);
#endif

    auto observationSpace = internal::convertToObservationSpace(observationsBinned);

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

    Matrix1D<double> E = Matrix1D<double>{observationSpace.data.size_q1};
    E = observationError;
    auto Robs = parameters.useLog ? diag(E) : diag(observationSpace.data * observationError);

    runKalmanFilter(
        forecast1D,
        analysisErrorCovariance,
        modelMatrix,
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

/* class Convection2DAnalysisCovariances {
public:
    Convection2DAnalysisCovariances(size_t V_size, size_t K_size);
    const std::vector<Matrix2D<double>>& operator[](size_t iV) const;
    std::vector<Matrix2D<double>>& operator[](size_t iV);
private:
    std::vector<std::vector<Matrix2D<double>>> _data;
}; */

data_assimilation::Convection2DAnalysisCovariances::Convection2DAnalysisCovariances(
    size_t P_size, size_t R_size, size_t V_size, size_t K_size)
    : _data(V_size, std::vector<Matrix2D<double>>(K_size)) {
    for (auto i = 0; i < V_size; ++i) {
        for (auto j = 0; j < K_size; ++j) {
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
    // auto Pf = M * Pa * transpose(M) + Q;

    // if (H.maxabs() > 0) {

    // forecast error covariance matrix [n x n]
    Matrix2D<double> Pf = abtrans(M * Pa, M);
    Pf += Q;
    // auto HT = transpose(H);    //transforms observations to model space [n x m]

    // auto Pf_times_HT = Pf * HT;
    Matrix2D<double> Pf_times_HT = abtrans(Pf, H);

    // Kalman matrix [n x m], Pf_times_HT will be overwritten and referenced now as K
    Matrix2D<double>& K = trans_solve(H * Pf_times_HT + R, Pf_times_HT);

    forecast += K * (obs - H * forecast);

    // std::cout << std::endl << "size 1: " << H.size_q1 << "    size 2: " << H.size_q2 << "       max: " << H.maxabs() << std::endl;

    Pa = Pf - K * H * Pf;
    //}
}

data_assimilation::Parameters data_assimilation::readParameters(const std::string& filename) {
    data_assimilation::Parameters result;
    ParametersIni parameters(filename);

    internal::getParameterBool(parameters, "run_data_assimilation",
                               result.runDataAssimilation, true);
    internal::getParameterBool(parameters, "useLog", result.useLog, true);
    parameters.getParameter("time_step", result.timeStep, true);
    parameters.getParameter("model_error", result.modelError, true);
    parameters.getParameter("observation_error", result.observationError, true);

    return result;
}
