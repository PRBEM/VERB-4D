// SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file DataAssimilationConvection.h
**/

#pragma once

#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <filesystem>
#include <exception>

#include "Matrix.h"
#include "UpdatableMatrix.h"
#include "Parameters.h"
#include "DataAssimilation.h"

namespace data_assimilation {

class Convection2DAnalysisCovariances {
   public:
    Convection2DAnalysisCovariances() {};
    Convection2DAnalysisCovariances(size_t P_size, size_t R_size, size_t V_size, size_t K_size);
    const std::vector<Matrix2D<double>>& operator[](size_t iV) const;
    std::vector<Matrix2D<double>>& operator[](size_t iV);

   private:
    std::vector<std::vector<Matrix2D<double>>> _data;
};

class DataAssimilationManagerConvection {
   public:
    DataAssimilationManagerConvection() {};
    DataAssimilationManagerConvection(const std::string& parametersFile,
                                      const Matrix2D<double>& P,
                                      const Matrix2D<double>& R,
                                      const Matrix2D<double>& V,
                                      const Matrix2D<double>& K,
                                      const std::string& debug_output_folder);

    void assimilate(double start_time, double end_time, Matrix4D<double>& PSD,
                    const Matrix4D<double>& VP, const Matrix4D<double>& VR,
                    const bool has_VX_updated,
                    const Matrix4D<double> &SaturationDensity,
                    const Matrix4D<double>& Loss, double dt);

    Convection2DAnalysisCovariances _analysisCovarianceConvection;

   private:
    std::string _debug_output_folder;
    Matrix2D<double> _P;
    Matrix2D<double> _R;
    Matrix2D<double> _V;
    Matrix2D<double> _K;
    std::vector<std::vector<Matrix2D<double>>> _model_matrix;

    bool _runDataAssimilation;
    Parameters _assimilationParameters;
    std::unique_ptr<DataSource> _dataSource;

};

struct DebugOutput2D {
    Matrix2D<double> observations;
    Matrix2D<double> forecast_init;
    Matrix2D<double> forecast_result;
    Matrix2D<double> a_priori_err_cov;
    Matrix2D<double> a_posteriori_err_cov;
    Matrix2D<double> model;
};

data_assimilation::DebugOutput2D runKalmanFilterConvection2D(
    Matrix2D<double>& forecast, Matrix2D<double>& analysisErrorCovariance, Matrix2D<double>& model_matrix,
    const Matrix2D<double>& P, const Matrix2D<double>& R,
    const Matrix2D<double>& VP, const Matrix2D<double>& VR,
    const bool has_VR_updated, const Matrix2D<double>& Loss, const Matrix2D<bool> &saturation_map,
    const Matrix2D<double>& observations, const double timeStep, const Parameters& parameters);

#ifdef DATA_ASSIMILATION_DEBUG
    struct DebugOutput4D {
        Matrix4D<double> observations;
        Matrix4D<double> forecast_init;
        Matrix4D<double> forecast_result;
        Matrix4D<double> a_priori_err_cov;
        Matrix4D<double> a_posteriori_err_cov;
        Matrix4D<double> model;

        DebugOutput4D(int P_size, int R_size, int V_size, int K_size) {
            observations = Matrix4D<double>(P_size, R_size, V_size, K_size);
            forecast_init = Matrix4D<double>(P_size, R_size, V_size, K_size);
            forecast_result = Matrix4D<double>(P_size, R_size, V_size, K_size);
            a_priori_err_cov = Matrix4D<double>(P_size * R_size, P_size * R_size, V_size, K_size);
            a_posteriori_err_cov = Matrix4D<double>(P_size * R_size, P_size * R_size, V_size, K_size);
            model = Matrix4D<double>(P_size * R_size, P_size * R_size, V_size, K_size);
        };

        void insert_output_2d(DebugOutput2D output_2d, int iV, int iK) {
            for (size_t iP = 0; iP < observations.size_w; iP++) {
                for (size_t iR = 0; iR < observations.size_x; iR++) {
                    observations[iP][iR][iV][iK] = output_2d.observations[iP][iR];
                    forecast_init[iP][iR][iV][iK] = output_2d.forecast_init[iP][iR];
                    forecast_result[iP][iR][iV][iK] = output_2d.forecast_result[iP][iR];
                }
            }
            for (size_t iP = 0; iP < a_priori_err_cov.size_w; iP++) {
                for (size_t iR = 0; iR < a_priori_err_cov.size_x; iR++) {
                    a_priori_err_cov[iP][iR][iV][iK] = output_2d.a_priori_err_cov[iP][iR];
                    a_posteriori_err_cov[iP][iR][iV][iK] = output_2d.a_posteriori_err_cov[iP][iR];
                    model[iP][iR][iV][iK] = output_2d.model[iP][iR];
                }
            }
        };

        void write_files(double time_step, std::string debug_output_folder){
            observations.writeToBinaryFile(debug_output_folder + "/observations" + std::to_string(time_step) + ".pltb");
            forecast_init.writeToBinaryFile(debug_output_folder + "/forecast_init_" + std::to_string(time_step) + ".pltb");
            forecast_result.writeToBinaryFile(debug_output_folder + "/forecast_result_" + std::to_string(time_step) + ".pltb");
            a_priori_err_cov.writeToBinaryFile(debug_output_folder + "/a_priori_err_cov" + std::to_string(time_step) + ".pltb");
            a_posteriori_err_cov.writeToBinaryFile(debug_output_folder + "/a_posteriori_err_cov" + std::to_string(time_step) + ".pltb");
            model.writeToBinaryFile(debug_output_folder + "/model_" + std::to_string(time_step) + ".pltb");
        };
    };
#endif

}  // namespace data_assimilation
