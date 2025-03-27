// SPDX-FileCopyrightText: 2015 UCLA
// SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * \file Convection_2D.cpp
 *
 *
 * Using the Convection_1D_Ultimate_QUICKEST6() for calculating P and R.
 * All equations and formulas for these calculations can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf .
 * Mostly coming from 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme.
 * In the source code mathematical equations are numbered corresponding to the numbering of equations found in the paper
 *
 * \brief Calculates the convection in 2D given a 2D matrix of Phase Space Densities, P, R, boundary conditions, diffusion, sources and losses
 */

#include "Convection_2D.h"

using namespace std;

/**
 * Calculation of 2D convection
 *
 * Uses Convection_1D_Ultimate_QUICKEST6.h for calculating P and R if either are of size > 3.
 *
 * All equations and formulas for these calculations can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf .
 * Mostly coming from 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme.
 * In the source code mathematical equations are numbered corresponding to the numbering of equations found in the paper.
 * Refer to the equations/function in Convection_1D_Ultimate_QUICKEST6.h
 *
 * @param PSD_PR - Phase Space Density - P and R held constant
 * @param P - Time (magnetic local time) based on phi
 * @param R - radial distance
 * @param P_size - dimension of P
 * @param R_size - dimension of R
 * @param P_LBC - P lower boundary
 * @param P_UBC - P upper boundary
 * @param R_LBC - R lower boundary
 * @param R_UBC - R upper boundary
 * @param P_LBC_type - Type of boundary for p lower
 * @param P_UBC_type - Type of boundary for p upper
 * @param R_LBC_type - Type of boundary for r lower
 * @param R_UBC_type - Type of boundary for r upper
 * @param VP - Matrix of energy and time(phi)
 * @param VR - Matrix of energy and radial distance
 * @param Sources - Sources matrix
 * @param Losses - Losses matrix (loss cone)
 * @param dt_total - total time change
 * @param min_PSD - minimum value from the PSD matrix
 * @param min_V - minimum value for V from the
 */
#ifdef SAVE_PSD_LOST_CONV
bool Convection_2D(
    Matrix2D<double>& PSD_PR,
    Matrix2D<double>& PSD_lost_PR,
    const Matrix2D<double>& P, const Matrix2D<double>& R,
    int P_size, int R_size,
    const Matrix1D<double>& P_LBC, const Matrix1D<double>& P_UBC,
    const Matrix1D<double>& R_LBC, const Matrix1D<double>& R_UBC,
    BoundaryConditionType P_LBC_type, BoundaryConditionType P_UBC_type,
    BoundaryConditionType R_LBC_type, BoundaryConditionType R_UBC_type,
    const Matrix2D<double>& VP, const Matrix2D<double>& VR,
    const Matrix2D<double>& Sources, const Matrix2D<double>& Losses,
    const Matrix2D<double>& G_conv, double dt_total) {
#else
bool Convection_2D(
    Matrix2D<double>& PSD_PR,
    const Matrix2D<double>& P, const Matrix2D<double>& R,
    int P_size, int R_size,
    const Matrix1D<double>& P_LBC, const Matrix1D<double>& P_UBC,
    const Matrix1D<double>& R_LBC, const Matrix1D<double>& R_UBC,
    BoundaryConditionType P_LBC_type, BoundaryConditionType P_UBC_type,
    BoundaryConditionType R_LBC_type, BoundaryConditionType R_UBC_type,
    const Matrix2D<double>& VP, const Matrix2D<double>& VR,
    const Matrix2D<double>& Sources, const Matrix2D<double>& Losses,
    const Matrix2D<double>& G_conv, double dt_total) {
#endif
    // indexes
    int iR, iP;

    double dP, dR, num_steps_P, num_steps_R, num_steps, dt;

    // Maximum Courant number. Can't be > 1, but can be smaller (0.25 is fairly common value)
    double maxCourNum = 0.25;

    // Find number of sub-time steps required to satisfy the Courant condition for both directions: P and R
    if (P_size >= 3) {
        dP = (P[1][0] - P[0][0]);  // FIXME - will work only for regular grid
        double num_steps_double = (double)dt_total * VP.maxabs() / dP / maxCourNum;
        num_steps_P = num_steps_double <= 1 ? 1 : ceil(num_steps_double);
    } else {
        num_steps_P = 1;
    }
    if (R_size >= 3) {
        dR = (R[0][1] - R[0][0]);  // FIXME - will work only for regular grid
        double num_steps_double = (double)dt_total * VR.maxabs() / dR / maxCourNum;
        num_steps_R = num_steps_double <= 1 ? 1 : ceil(num_steps_double);
    } else {
        num_steps_R = 1;
    }

    // max_dt_P = dP / DP.maxabs();
    // max_dt_R = dP / DR.maxabs();
    // max_dt = (max_dt_P < max_dt_R)? max_dt_P : max_dt_R;
    // num_steps = (max_dt < dt_total) ? ceil((double)dt_total/max_dt) : 1;

    // Either use the smallest time step for both, or specify the maximum time step here and then different time steps will be used
    // (the Courant condition will be checked inside of Convection_1D one more time)
    // "max" is more accurate, but "min" is much-much faster and (hopefully?) still more accurate than completely unrelated time steps
#ifdef FAST_CONVECTION
    num_steps = std::min(num_steps_P, num_steps_R);
#else
    num_steps = std::max(num_steps_P, num_steps_R);
#endif

    dt = dt_total / num_steps;

    // 1d-slices of PSD to pass into Convection_1D
    Matrix1D<double> PSD_P(P_size);
    Matrix1D<double> P_P(P_size);
    Matrix1D<double> VP_P(P_size);

    Matrix1D<double> PSD_R(R_size);
    Matrix1D<double> R_R(R_size);
    Matrix1D<double> VR_R(R_size);
    Matrix1D<double> G_conv_R(R_size);

    Matrix2D<double> PSD_PR_G;

    // precompute losses
    Matrix2D<double> losses_exp = Losses.exp(dt);
    for (int it = 0; it < num_steps; it++) {
        if (P_size >= 3) {
            for (iR = 0; iR < R_size - 1; iR++) {
                // 1d slice
                PSD_PR.ySlice(PSD_P, iR);
                P.ySlice(P_P, iR);
                VP.ySlice(VP_P, iR);

                Convection_1D_ULTIMATE_QUICKEST6(
                    PSD_P, P_P, P_size,
                    P_LBC[iR], P_UBC[iR],
                    P_LBC_type, P_UBC_type,
                    VP_P, dt);

                // copy results back
                for (iP = 0; iP < P_size; iP++)
                    PSD_PR[iP][iR] = PSD_P[iP];
            }
        }

        PSD_PR_G = PSD_PR.times(G_conv);
        //PSD_PR_G = PSD_PR;

        if (R_size >= 3) {
            for (iP = 0; iP < P_size; iP++) {
                // 2d slice
                PSD_PR_G.xSlice(PSD_R, iP);
                R.xSlice(R_R, iP);
                VR.xSlice(VR_R, iP);
                G_conv.xSlice(G_conv_R, iP);

                Convection_1D_ULTIMATE_QUICKEST6(
                    PSD_R, R_R, R_size,
                    R_LBC[iP] * G_conv_R[0], 
                    R_UBC[iP] * G_conv_R[R_size-1],  
                    R_LBC_type, R_UBC_type,
                    VR_R, dt);

                // copy results back
                for (iR = 0; iR < R_size; iR++){
                    PSD_PR[iP][iR] = PSD_R[iR] / G_conv_R[iR];
                }
            } 
        }

        // Losses
        for (iP = 0; iP < P_size; iP++) {
            for (iR = 0; iR < R_size; iR++) {
#ifdef SAVE_PSD_LOST_CONV
                double PSD_tmp = PSD_PR[iP][iR];
                PSD_PR[iP][iR] = PSD_PR[iP][iR] * losses_exp[iP][iR];
                PSD_lost_PR[iP][iR] += (PSD_tmp - PSD_PR[iP][iR]);
#else
                PSD_PR[iP][iR] = PSD_PR[iP][iR] * losses_exp[iP][iR];
#endif
            }
        }

        // Sources
        for (iP = 0; iP < P_size; iP ++) {
            for (iR = 1; iR < R_size - 1; iR++) {
                PSD_PR[iP][iR] = PSD_PR[iP][iR] + Sources[iP][iR] * dt;
            }
        }
    }

    return true;
}
