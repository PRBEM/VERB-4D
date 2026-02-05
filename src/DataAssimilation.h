#pragma once

#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <filesystem>

#include "Matrix.h"
#include "UpdatableMatrix.h"
#include "PMF.h"
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
};

class DataServerDataSource : public DataSource {
public:
    DataServerDataSource(const std::string& satellite_lst_file);

    virtual Matrix4D<double> getObservations(
        const double timeStart, const double timeEnd,
        const Matrix2D<double>& P, const Matrix2D<double>& R,
        const Matrix2D<double>& V, const Matrix2D<double>& K) override;

private:
    std::vector<pmf::Parameters> _dataParameters;
};

class LocalFilesDataSource : public DataSource {
public:
    LocalFilesDataSource(const std::string& satellite_lst_file, 
            size_t size_q1, size_t size_q2, size_t size_q3, size_t size_q4);

    virtual Matrix4D<double> getObservations(
        const double timeStart, const double timeEnd,
        const Matrix2D<double>& P, const Matrix2D<double>& R,
        const Matrix2D<double>& V, const Matrix2D<double>& K) override;

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
class Convection2DAnalysisCovariances {
   public:
    Convection2DAnalysisCovariances(size_t P_size, size_t R_size, size_t V_size, size_t K_size);
    const std::vector<Matrix2D<double>>& operator[](size_t iV) const;
    std::vector<Matrix2D<double>>& operator[](size_t iV);

   private:
    std::vector<std::vector<Matrix2D<double>>> _data;
};
class DataAssimilationManagerConvection {
   public:
    DataAssimilationManagerConvection(const std::string& parametersFile,
                                      double time_start, double time_end,
                                      const Matrix2D<double>& P,
                                      const Matrix2D<double>& R,
                                      const Matrix2D<double>& V,
                                      const Matrix2D<double>& K,
                                      const std::string& debug_output_folder);

    void assimilate(double time, Matrix4D<double>& PSD,
                    const Matrix4D<double>& VP, const Matrix4D<double>& VR,
                    const bool has_VX_updated,
                    const Matrix4D<double> &SaturationDensity,
                    const Matrix4D<double>& Loss, double dt);

    Convection2DAnalysisCovariances _analysisCovarianceConvection;

   private:
    std::string _debug_output_folder;
    double _timeStart;
    double _timeEnd;
    Matrix2D<double> _P;
    Matrix2D<double> _R;
    Matrix2D<double> _V;
    Matrix2D<double> _K;
    std::vector<std::vector<Matrix2D<double>>> _model_matrix;

    bool _runDataAssimilation;
    double _timePrev;
    double _timeNext;
    Parameters _assimilationParameters;
    std::unique_ptr<DataSource> _dataSource;

};

Parameters readParameters(const std::string& filename);

std::vector<std::vector<Observations>> getObservations(
    double timeStart,
    double timeEnd,
    const Matrix2D<double>& V,
    const Matrix2D<double>& K,
    const std::vector<pmf::Parameters>& parameters);

struct DebugOuput2D {
    Matrix2D<double> observations;
    Matrix2D<double> forecast_init;
    Matrix2D<double> forecast_result;
    Matrix2D<double> error_covariance;
    Matrix2D<double> model;
};

data_assimilation::DebugOuput2D runKalmanFilterConvection2D(
    Matrix2D<double>& forecast, Matrix2D<double>& analysisErrorCovariance, Matrix2D<double>& model_matrix,
    const Matrix2D<double>& P, const Matrix2D<double>& R,
    const Matrix2D<double>& VP, const Matrix2D<double>& VR,
    const bool has_VR_updated, const Matrix2D<double>& Loss, const Matrix2D<bool> &saturation_map,
    const Matrix2D<double>& observations, const double timeStep, const Parameters& parameters);

/**
 * @brief Applies the Kalman filter updating forecast and Pa
 * @param[in, out]  forecast    the current forecast [n x 1]
 * @param[in, out]  Pa          analysis error covariance matrix [n x n]
 * @param[in]       M           model matrix operator [n x n]
 * @param[in]       Q           model error corvariance matrix [n x n]
 * @param[in]       obs         observations [m x 1]
 * @param[in]       H           model operator, transform forecast to observation space [m x n]
 * @param[in]       R           observation error covariance matrix [m x m]
 */
void runKalmanFilter(
    Matrix1D<double>& forecast,
    Matrix2D<double>& Pa,
    const Matrix2D<double>& M, const Matrix2D<double>& Q,
    const Matrix1D<double>& obs, const Matrix2D<double>& H,
    const Matrix2D<double>& R);

#ifdef DATA_ASSIMILATION_DEBUG
    struct DebugOuput4D {
        Matrix4D<double> observations;
        Matrix4D<double> forecast_init;
        Matrix4D<double> forecast_result;
        Matrix4D<double> error_covariance;
        Matrix4D<double> model;

        DebugOuput4D(int P_size, int R_size, int V_size, int K_size) {
            observations = Matrix4D<double>(P_size, R_size, V_size, K_size);
            forecast_init = Matrix4D<double>(P_size, R_size, V_size, K_size);
            forecast_result = Matrix4D<double>(P_size, R_size, V_size, K_size);
            error_covariance = Matrix4D<double>(P_size * R_size, P_size * R_size, V_size, K_size);
            model = Matrix4D<double>(P_size * R_size, P_size * R_size, V_size, K_size);
        };

        void insert_output_2d(DebugOuput2D output_2d, int iV, int iK) {
            for (size_t iP = 0; iP < observations.size_w; iP++) {
                for (size_t iR = 0; iR < observations.size_x; iR++) {
                    observations[iP][iR][iV][iK] = output_2d.observations[iP][iR];
                    forecast_init[iP][iR][iV][iK] = output_2d.forecast_init[iP][iR];
                    forecast_result[iP][iR][iV][iK] = output_2d.forecast_result[iP][iR];
                }
            }
            for (size_t iP = 0; iP < error_covariance.size_w; iP++) {
                for (size_t iR = 0; iR < error_covariance.size_x; iR++) {
                    error_covariance[iP][iR][iV][iK] = output_2d.error_covariance[iP][iR];
                    model[iP][iR][iV][iK] = output_2d.model[iP][iR];
                }
            }
        };

        void write_files(double time_step, std::string debug_output_folder){
            observations.writeToBinaryFile(debug_output_folder + "/observations" + std::to_string(time_step) + ".pltb");
            forecast_init.writeToBinaryFile(debug_output_folder + "/forecast_init_" + std::to_string(time_step) + ".pltb");
            forecast_result.writeToBinaryFile(debug_output_folder + "/forecast_result_" + std::to_string(time_step) + ".pltb");
            error_covariance.writeToBinaryFile(debug_output_folder + "/error_covariance_" + std::to_string(time_step) + ".pltb");
            model.writeToBinaryFile(debug_output_folder + "/model_" + std::to_string(time_step) + ".pltb");
        };
    };
#endif

}  // namespace data_assimilation
