// SPDX-FileCopyrightText: 2015 UCLA
// SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * \file Convection_3D.cpp
 *
 *
 * Using the Convection_1D_Ultimate_QUICKEST6() for calculating P and R.
 *  * Using the Convection_1D_2ndORDER_NONUNIFORM_GRID() for calculating V.
 * All equations and formulas for these calculations can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf .
 * Mostly coming from 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme.
 * In the source code mathematical equations are numbered corresponding to the numbering of equations found in the paper
 *
 * \brief Calculates the convection in 2D given a 2D matrix of Phase Space Densities, P, R, boundary conditions, diffusion, sources and losses
 */

#include "Convection_3D.h"

using namespace std;

/**
 * Calculation of 3D convection
 *
 * Uses Convection_1D_Ultimate_QUICKEST6.h for calculating P and R if either are of size > 3.
 * Uses Convection_1D_2ndORDER_NONUNIFORM_GRID.h for calculating V if size > 3.
 *
 * All equations and formulas for these calculations can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf .
 * Mostly coming from 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme.
 * In the source code mathematical equations are numbered corresponding to the numbering of equations found in the paper.
 * Refer to the equations/function in Convection_1D_Ultimate_QUICKEST6.h
 *
 * @param PSD_PRV - Phase Space Density - P, R and V held constant
 * @param P - Time (magnetic local time) based on phi
 * @param R - radial distance
 * @param V - V variable
 * @param P_size - dimension of P
 * @param R_size - dimension of R
 * @param V_size - dimension of V
 * @param P_LBC - P lower boundary
 * @param P_UBC - P upper boundary
 * @param R_LBC - R lower boundary
 * @param R_UBC - R upper boundary
 * @param V_LBC - V lower boundary
 * @param V_UBC - V upper boundary
 * @param P_LBC_type - Type of boundary for p lower
 * @param P_UBC_type - Type of boundary for p upper
 * @param R_LBC_type - Type of boundary for r lower
 * @param R_UBC_type - Type of boundary for r upper
 * @param VP - Matrix of energy and time(phi)
 * @param VR - Matrix of energy and radial distance
 * @param VV - Matrix of Coulomb velocity
 * @param Sources - Sources matrix
 * @param Losses - Losses matrix (loss cone)
 * @param dt_total - total time change
 * @param min_PSD - minimum value from the PSD matrix
 * @param min_V - minimum value for V from the
 */

bool Convection_3D(
    Matrix3D<double>& PSD_PRV,
    const Matrix3D<double>& P, const Matrix3D<double>& R, const Matrix3D<double>& V,
    int P_size, int R_size, int V_size,
    const Matrix2D<double>& P_LBC, const Matrix2D<double>& P_UBC,
    const Matrix2D<double>& R_LBC, const Matrix2D<double>& R_UBC,
	const Matrix2D<double>& V_LBC, const Matrix2D<double>& V_UBC,
    BoundaryConditionType P_LBC_type, BoundaryConditionType P_UBC_type,
    BoundaryConditionType R_LBC_type, BoundaryConditionType R_UBC_type,
	BoundaryConditionType V_LBC_type, BoundaryConditionType V_UBC_type,
    const Matrix3D<double>& VP, const Matrix3D<double>& VR, const Matrix3D<double>& VV,
    const Matrix3D<double>& Sources, const Matrix3D<double>& Losses,
    const Matrix3D<double>& G_conv, double dt_total) {

    // indexes
    int iR, iP, iV; 

	// Note: dV is not forgotten; the courant number will be checked in 1D_2ndORDER  
    double dP, dR, num_steps_P, num_steps_R, num_steps, dt;

    // Maximum Courant number. Can't be > 1, but can be smaller (0.25 is fairly common value)
    double maxCourNum = 0.25;

    // Find number of sub-time steps required to satisfy the Courant condition for both directions: P and R
    if (P_size >= 3) {
        dP = (P[1][0][0] - P[0][0][0]);  // FIXME - will work only for regular grid
        double num_steps_double = (double)dt_total * VP.maxabs() / dP / maxCourNum;
        num_steps_P = num_steps_double <= 1 ? 1 : ceil(num_steps_double);
    } else {
        num_steps_P = 1;
    }
    if (R_size >= 3) {
        dR = (R[0][1][0] - R[0][0][0]);  // FIXME - will work only for regular grid
        double num_steps_double = (double)dt_total * VR.maxabs() / dR / maxCourNum;
        num_steps_R = num_steps_double <= 1 ? 1 : ceil(num_steps_double);
    } else {
        num_steps_R = 1;
    }

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

	Matrix1D<double> PSD_V(V_size);
    Matrix1D<double> V_V(V_size);
    Matrix1D<double> VV_V(V_size);

    Matrix3D<double> PSD_PRV_GR;
	Matrix3D<double> PSD_PRV_GV;


    // precompute losses
    Matrix3D<double> losses_exp = Losses.exp(dt);
    for (int it = 0; it < num_steps; it++) {

    if (P_size >= 3) {
        for (iR = 0; iR < R_size - 1; iR++) {
			for (iV = V_size - 1; iV >=  0; iV--) {
				// 1d slice
				PSD_PRV.yzSlice(PSD_P, iR, iV);
				P.yzSlice(P_P, iR, iV);
				VP.yzSlice(VP_P, iR, iV);
				// Speedup: if PSD~=0 or V~=0, skip the thing
				//if (PSD_P.max() < min_PSD || VP_P.maxabs() < min_V)
				//if (PSD_P.max() < min_PSD)
				//	continue;

				Convection_1D_ULTIMATE_QUICKEST6(
					PSD_P, P_P, P_size,
					P_LBC[iR][iV], P_UBC[iR][iV],
					P_LBC_type, P_UBC_type,
					VP_P, dt);

				// copy results back
				for (iP = 0; iP < P_size; iP++) {
					PSD_PRV[iP][iR][iV] = PSD_P[iP];
				}
			}
        }
    }

    PSD_PRV_GR = PSD_PRV.times(G_conv);
    if (R_size >= 3) {
        for (iP = 0; iP < P_size; iP++) {
			for (iV = V_size - 1; iV >=  0; iV--) {
				// 2d slice
				PSD_PRV_GR.xzSlice(PSD_R, iP, iV);
				R.xzSlice(R_R, iP, iV);
				VR.xzSlice(VR_R, iP, iV);

				Convection_1D_ULTIMATE_QUICKEST6(
					PSD_R, R_R, R_size,
					R_LBC[iP][iV] * G_conv_R[0], 
					R_UBC[iP][iV] * G_conv_R[R_size-1],  
					R_LBC_type, R_UBC_type,
					VR_R, dt);

				// copy results back
				for (iR = 0; iR < R_size; iR++) {
					PSD_PRV[iP][iR][iV] = PSD_R[iR] / G_conv_R[iR];
				}
			}
        } 
    }
    

    PSD_PRV_GV = PSD_PRV.times(V.sqrt());
    if (V_size >= 3) {
         for (iP = 0; iP < P_size; iP++) {
		  	for (iR = 0; iR < R_size; iR++) {
				// 2d slice
				PSD_PRV_GV.xySlice(PSD_V, iP, iR);
				V.xySlice(V_V, iP, iR);
				VV.xySlice(VV_V, iP, iR);

				Convection_1D_2ndORDER_NONUNIFORM_GRID(
					PSD_V, V_V, V_size,
					V_LBC[iP][iR] * std::sqrt(V_V[0]), 
					V_UBC[iP][iR] * std::sqrt(V_V[V_size-1]),  
					V_LBC_type, V_UBC_type,
					VV_V, dt);

				// copy results back
				for (iV = 0; iV < V_size; iV++) {
					PSD_PRV[iP][iR][iV] = PSD_V[iV] / std::sqrt(V_V[iV]);
				}
	 		}
         } 
    } 
    
        // Losses
        for (iP = 0; iP < P_size; iP++) {
            for (iR = 0; iR < R_size; iR++) {
				for (iV = 0; iV < V_size; iV++) {
                	PSD_PRV[iP][iR][iV] = PSD_PRV[iP][iR][iV] * losses_exp[iP][iR][iV];
				}
            }
        }

        // Sources
        for (iP = 0; iP < P_size; iP ++) {
            for (iR = 1; iR < R_size - 1; iR++) {
                for (iV = 0; iV < V_size; iV++) {
                    PSD_PRV[iP][iR][iV] = PSD_PRV[iP][iR][iV] + Sources[iP][iR][iV] * dt;
                }
            }
        }

    }

    return true;
}
