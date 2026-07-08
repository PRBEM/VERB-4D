// SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file DataAssimilationDiffusion.h
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
#include "DataAssimilation.h"
#include "BoundaryConditionType.hpp"

namespace data_assimilation {

class Diffusion1DAnalysisCovariances {
   public:
    Diffusion1DAnalysisCovariances() {};
    Diffusion1DAnalysisCovariances(size_t R_size, size_t V_size, size_t K_size);
    const std::vector<Matrix2D<double>>& operator[](size_t iV) const;
    std::vector<Matrix2D<double>>& operator[](size_t iV);

   private:
    std::vector<std::vector<Matrix2D<double>>> _data;
};

class DataAssimilationManagerDiffusion {
   public:
    DataAssimilationManagerDiffusion() {};
    DataAssimilationManagerDiffusion(const std::string& parametersFile,
                                     const Matrix1D<double>& R,
                                     const Matrix2D<double>& V,
                                     const Matrix2D<double>& K,
                                     const std::string& debug_output_folder);

    void assimilate(double start_time,
                    double end_time,
                    Matrix4D<double>& PSD,
                    const Matrix4D<double>& DLL,
                    const Matrix4D<double>& G,
                    Matrix3D<double>& x_LBC, Matrix3D<double>& x_UBC,
                    BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
                    double dt);

    Diffusion1DAnalysisCovariances _analysisCovariance;

   private:
    std::string _debug_output_folder;
    Matrix1D<double> _R;
    Matrix2D<double> _V;
    Matrix2D<double> _K;
    std::vector<std::vector<Matrix2D<double>>> _model_matrix;

    bool _runDataAssimilation;
    Parameters _assimilationParameters;
    std::unique_ptr<DataSource> _dataSource;

};

struct DebugOutput2DDiffusion {
    Matrix1D<double> observations;
    Matrix1D<double> forecast_init;
    Matrix1D<double> forecast_result;
    Matrix2D<double> a_priori_err_cov;
    Matrix2D<double> a_posteriori_err_cov;
    Matrix2D<double> model;
};

data_assimilation::DebugOutput2DDiffusion runKalmanFilterDiffusion1D(
    Matrix1D<double>& forecast,
    Matrix2D<double>& analysisErrorCovariance,
    Matrix2D<double>& model_matrix,
    const Matrix1D<double>& R,
    const Matrix1D<double>& DLL,
    const Matrix1D<double>& G,
	double x_LBC, double x_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
    const Matrix1D<double>& observations,
    const double timeStep,
    const Parameters& parameters);

#ifdef DATA_ASSIMILATION_DEBUG
    struct DebugOutput3D {
        Matrix3D<double> observations;
        Matrix3D<double> forecast_init;
        Matrix3D<double> forecast_result;
        Matrix4D<double> a_priori_err_cov;
        Matrix4D<double> a_posteriori_err_cov;
        Matrix4D<double> model;

        DebugOutput3D(int R_size, int V_size, int K_size) {
            observations = Matrix3D<double>(R_size, V_size, K_size);
            forecast_init = Matrix3D<double>(R_size, V_size, K_size);
            forecast_result = Matrix3D<double>(R_size, V_size, K_size);
            a_priori_err_cov = Matrix4D<double>(R_size, R_size, V_size, K_size);
            a_posteriori_err_cov = Matrix4D<double>(R_size, R_size, V_size, K_size);
            model = Matrix4D<double>(R_size, R_size, V_size, K_size);
        };

        void insert_output_2d(DebugOutput2DDiffusion output_2d, int iV, int iK) {
                        
            for (size_t iR = 0; iR < observations.size_q1; iR++) {
                observations[iR][iV][iK] = output_2d.observations[iR];
                forecast_init[iR][iV][iK] = output_2d.forecast_init[iR];
                forecast_result[iR][iV][iK] = output_2d.forecast_result[iR];
            }

            size_t R_size_2D = a_priori_err_cov.size_x;

            for (size_t iR = 0; iR < R_size_2D; iR++) {
                    for (size_t iC = 0; iC < R_size_2D; iC++) {
                        a_priori_err_cov[iR][iC][iV][iK] = output_2d.a_priori_err_cov[iR][iC];
                        a_posteriori_err_cov[iR][iC][iV][iK] = output_2d.a_posteriori_err_cov[iR][iC];
                        model[iR][iC][iV][iK] = output_2d.model[iR][iC];
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
