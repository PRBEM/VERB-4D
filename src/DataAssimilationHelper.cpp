// SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file DataAssimilationHelper.cpp
**/

#include "DataAssimilationHelper.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "CustomDate.h"
#include "DataAssimilation.h"
#include "MatrixOperations.h"
#include "PMF.h"

namespace internal = data_assimilation::internal;

Matrix2D<double> internal::getModelMatrixConvection2DNoStabilityCheck(
    const Matrix2D<double> &VP,
    const Matrix2D<double> &VR,
    const Matrix2D<double> &Loss,
    const Matrix2D<bool> &saturation_map,
    double timeStep,
    double dP,
    double dR,
    bool use_log,
    double correlation_factor) {

    int P_size = VP.size_q1;
    int R_size = VP.size_q2;

    Matrix2D<double> result(P_size * R_size, P_size * R_size);
    result = 0.;

    for (auto iP = 0; iP < P_size; ++iP) {
        for (auto iR = 0; iR < R_size; ++iR) {

            // just write boundary
            if (iR == 0 || iR == R_size - 1) {
                result[iR + iP * R_size][iR + iP * R_size] = correlation_factor;
                continue;
            }

            auto iP_c = iP;
            size_t iP_l, iP_r;
            if (iP == 0 || iP == P_size - 1) {
                iP_l = P_size - 2;
                iP_r = 1;
            } else {
                iP_l = iP - 1;
                iP_r = iP + 1;
            }

            auto iR_c = iR;
            auto iR_l = iR - 1;
            auto iR_r = iR + 1;

            auto VP_c = VP[iP_c][iR_c];
            auto VR_c = VR[iP_c][iR_c];
            auto Loss_c = Loss[iP_c][iR_c];
            bool sat_c = saturation_map[iP_c][iR_c];
            Loss_c = sat_c ? 0 : Loss_c;

            double m_cc = 0;

            if (use_log) {
                m_cc = 1. - std::abs(VP_c) * timeStep / dP - std::abs(VR_c) * timeStep / dR;
            } else {
                m_cc = 1. - std::abs(VP_c) * timeStep / dP - std::abs(VR_c) * timeStep / dR + timeStep * Loss_c;
            }
    
            if (m_cc < 0) {
                std::cout << "iR: " << iR << std::endl;
                std::cout << "sat_c: " << sat_c << std::endl;
                std::cout << "VP: " << std::abs(VP_c) * timeStep / dP << std::endl;
                std::cout << "VR: " << std::abs(VR_c) * timeStep / dR << std::endl;
                std::cout << "Loss: " << timeStep * Loss_c << std::endl;
                std::cout << m_cc << std::endl;
                throw std::exception();
            } 

            double m_lc, m_rc, m_cl, m_cr;
            if (VP_c > 0.) {
                m_lc = VP_c * timeStep / dP;
                m_rc = 0.;
            } else {
                m_lc = 0.;
                m_rc = -VP_c * timeStep / dP;
            }

            if (VR_c > 0.) {
                m_cl = VR_c * timeStep / dR;
                m_cr = 0.;
            } else {
                m_cl = 0.;
                m_cr = -VR_c * timeStep / dR;
            }

            if (m_cc < 0 || m_lc < 0 || m_rc < 0 || m_cl < 0 || m_cr < 0) {
                std::cout << "This model matrix is unstable. Probably due to a misscalculation in the CLF condition." << std::endl;
                exit(EXIT_FAILURE);
            }

            result[iR_c + iP_c * R_size][iR_c + iP_c * R_size] = m_cc;
            result[iR_c + iP_c * R_size][iR_c + iP_l * R_size] = m_lc;
            result[iR_c + iP_c * R_size][iR_c + iP_r * R_size] = m_rc;
            result[iR_c + iP_c * R_size][iR_l + iP_c * R_size] = m_cl;
            result[iR_c + iP_c * R_size][iR_r + iP_c * R_size] = m_cr;

        }
    }

    return result;
}

std::pair<int, double> internal::splitTimeStepCourantCondition(
    double maximumCourantNumber,
    double timeStep,
    const Matrix2D<double>& VP,
    const Matrix2D<double>& VR,
    double dP,
    double dR) {

    Matrix2D<double> c_p = VP * timeStep / dP;
    Matrix2D<double> c_r = VR * timeStep / dR;
    Matrix2D<double> c_pr = c_p.abs() + c_r.abs();

    double cour_num = c_pr.maxabs();

    int nt = (cour_num <= maximumCourantNumber)
                 ? 1
                 : std::ceil(cour_num / maximumCourantNumber);

    double timeStep_new = timeStep / static_cast<double>(nt);

    return {nt, timeStep_new};
}

Matrix2D<double> internal::getModelMatrixConvection2D(
    const Matrix2D<double> &VP,
    const Matrix2D<double> &VR,
    const Matrix2D<double> &Loss,
    const Matrix2D<bool> &saturation_map,
    double timeStep,
    double dP,
    double dR,
    bool use_log,
    double correlation_time) {

    const double maximumCourantNumber = 1;
    auto stepSplit = splitTimeStepCourantCondition(
        maximumCourantNumber, timeStep,
        VP, VR, dP, dR);
    int &nt = stepSplit.first;
    const double &timeStep_new = stepSplit.second;

    Matrix2D<double> result_1step = getModelMatrixConvection2DNoStabilityCheck(
        VP, VR, Loss, saturation_map, timeStep_new, dP, dR, use_log, exp(-timeStep_new/correlation_time));
    
    // std::cout << "NUMBER OF TIMESTEPS: " << nt << std::endl;

    Matrix2D<double> result = result_1step^nt;

    return result;
}

void internal::getModelMatrixDiffusion1D(
    Matrix2D<double> &F,
    const Matrix1D<double>& R,
    const Matrix1D<double>& DLL,
    const Matrix1D<double>& G,
    const double x_LBC,
    const double x_UBC,
    const BoundaryConditionType x_LBC_type,
    const BoundaryConditionType x_UBC_type,
    const double timeStep) {

    Matrix1D<double> zero_mat{R.size_q1};
    zero_mat = 0;

    std::tuple<CalculationMatrix, CalculationMatrix, CalculationMatrix> matrices = 
        getCalculationMatrix1DDiffusion(R, R.size_q1, x_LBC, x_UBC, x_LBC_type, x_UBC_type, DLL, G, zero_mat, zero_mat, timeStep);

    CalculationMatrix A = std::get<0>(matrices);

    Matrix2D<double> matr_A = std::get<0>(matrices).toMatrix2D();
    Matrix2D<double> matr_B = std::get<1>(matrices).toMatrix2D();

    F = trans_solve(matr_A, matr_B);

}

double internal::conditionalMean(
    const Matrix1D<double>& a,
    const Matrix1D<bool>& flag) {
    if (a.size_q1 != flag.size_q1) {
        std::cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        std::cout << "Matrix sizes are different. ";
        std::cout << "a.size_q1 = " << a.size_q1 << ", ";
        std::cout << "flag.size_q1 = " << flag.size_q1 << std::endl;
        ;
        exit(EXIT_FAILURE);
    }

    double result{0.};
    int counter{0};
    for (size_t i = 0; i < a.size_q1; ++i) {
        if (flag[i] && !std::isnan(a[i])) {
            result += a[i];
            ++counter;
        }
    }

    if (!counter) {
        result = NAN;
    } else {
        result /= counter;
    }

    return result;
}

data_assimilation::ObservationSpace internal::convertToObservationSpace(
    const Matrix2D<double>& data) {
    size_t size{0};
    for (size_t iP = 0; iP < data.size_q1; ++iP) {
        for (size_t iR = 0; iR < data.size_q2; ++iR) {
            if (std::isfinite(data[iP][iR])) {
                ++size;
            }
        }
    }

    if (size == 0) {
        ObservationSpace result{
            Matrix1D<double>{1},
            Matrix2D<double>{1, data.size_q1 * data.size_q2}};
        result.H = 0.0;
        result.data = 0.0;

        return result;
    }

    ObservationSpace result{
        Matrix1D<double>{size},
        Matrix2D<double>{size, data.size_q1 * data.size_q2}};

    result.H = 0.;

    int counter{0};
    for (size_t iP = 0; iP < data.size_q1; ++iP) {
        for (size_t iR = 0; iR < data.size_q2; ++iR) {
            int H_idx = iR + iP * data.size_q2;
            if (std::isfinite(data[iP][iR])) {
                result.data[counter] = data[iP][iR]; // +1 to account for the boundary 
                result.H[counter][H_idx] = 1.0;
                ++counter;
            }
        }
    }
    return result;
}

data_assimilation::ObservationSpace internal::convertToObservationSpace(
    const Matrix1D<double>& data) {
    size_t size{0};
    for (size_t iR = 0; iR < data.size_q1; ++iR) {
        if (std::isfinite(data[iR])) {
            ++size;
        }
    }

    if (size == 0) {
        ObservationSpace result{
            Matrix1D<double>{1},
            Matrix2D<double>{1, data.size_q1}};
        result.H = 0.0;
        result.data = 0.0;

        return result;
    }

    ObservationSpace result{
        Matrix1D<double>{size},
        Matrix2D<double>{size, data.size_q1}};

    result.H = 0.;

    int counter{0};
    for (size_t iR = 0; iR < data.size_q1; ++iR) {
        if (std::isfinite(data[iR])) {
            result.data[counter] = data[iR];
            result.H[counter][iR] = 1.0;
            ++counter;
        }
    }
    return result;
}
