/**
 * \file Convection_1D_ULTIMATE_QUICKEST6.cpp
 *
 * Leonard, 1991; Leonard and Niknafs, 1991;
 * 
 * Leonard BP (1988) Universal Limiter for transient interpolation modeling of the advective transport equations: the ULTIMATE conservative difference scheme, NASA technical Memorandum 100916 ICOMP-88-11
 *
 * It works somehow, edit with a great care! - edited by Hayley to correct and issue with the PSD update. Without this correction, the total PSD will gain and drop periodically with the drift
 *
 * All equations and formulas for these calculations can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf .
 * Mostly coming from 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme.
 * In the source code mathematical equations are numbered corresponding to the numbering of equations found in the paper
 * 
 * \brief Calculates the convection in 1D given a 1D matrix of Phase Space Densities, boundary conditions, diffusion, sources and losses
 */

#include "Convection_1D_ULTIMATE_QUICKEST6.h"
#include "BoundaryConditionType.hpp"

/** Using the Standard namespace
*/
using namespace std;

#define gst 5
#define uw_n 4
#define dw_n 5

bool Convection_1D_ULTIMATE_QUICKEST6( 
		Matrix1D<double> &PSD, 
		const Matrix1D<double>& x,
		int x_size,
		double x_LBC, double x_UBC,
		BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
		const Matrix1D<double>& Ux,
		const Matrix1D<double>& Sources, const Matrix1D<double>& Losses,
		double dt_total)
{
    int ix;
    [[maybe_unused]] double CURVAV, GRAD;
    double DEL, CURV, ADEL, ACURV, PSD_REF;  // temporary variables
    Matrix1D<double> PSD_f(x_size);                        // PSD on a 'face' (between nodes)

    Matrix1D<double> CourNum(x_size);  // Courant number

    double dx = (x[1] - x[0]);  // FIXME - will work only for regular (uniform) grid
    // double max_dt = dx / Ux.maxabs();
    // int num_steps = (max_dt < dt_total) ? ceil((double)dt_total/max_dt) : 1;

    double maxCourNum = 1;
    double ux_max = Ux.maxabs();
    double steps_double = (double)dt_total * ux_max / dx / maxCourNum;
    int num_steps = (steps_double <= 1) ? 1 : ceil(steps_double);

    double dt = dt_total / num_steps;

    CourNum = Ux / dx * dt;
    const double courant_max = std::abs(ux_max / dx * dt);
    if (courant_max > 1) {
        Logger::error << "max(CourNum) = " << courant_max << ", calculation can't be performed." << endl;
        exit(EXIT_FAILURE);
    }

    double d1, d3, d5, d7, d9;
    [[maybe_unused]] double D10;
    double D2, D4, D6, D8;
    [[maybe_unused]] double dif0;
    double dif1, dif2, dif3, dif4, dif5, dif6;  // Leonard and Niknfas, 1991
    [[maybe_unused]] double THC1 = 1e6, THC2 = 1e12, THG = 1e6;        // Leonard and Niknfas, 1991

    // PSD down-wind, center, and up-wind
    // double PSD_U[uw_n], PSD_D[dw_n], PSD_C;
    double PSD_U[uw_n + 1], PSD_D[dw_n + 1], PSD_C;

    if (x_size < gst * 2 + 1) {
        Logger::error << "Number of convection grid can't be smaller than " << (gst * 2 + 1) << " due to numerical method used." << endl;
        // exit(EXIT_FAILURE);
    }

    bool use_discreminator = true, use_limiting = true;
    [[maybe_unused]] bool never_3rd_order = false, always_3rd_order = false;
    [[maybe_unused]] bool never_5rd_order = false, always_5rd_order = false;
    [[maybe_unused]] bool never_7rd_order = false, always_7rd_order = false;

    constexpr double three_factorial = 3 * 2 * 1;
    constexpr double five_factorial  = 5 * 4 * three_factorial;
    constexpr double seven_factorial = 7 * 6 * five_factorial;
    constexpr double nine_factorial  = 9 * 8 * seven_factorial;

    // Temporary PSD to use during calculations, includes ghost points
    Matrix1D<double> PSD_t(x_size + gst + gst);

    int ig;  // ghost points iterator
    for (int it = 0; it < num_steps; it++) {
        // copying PSD into PSD_t offset by gst, with gst extra spots at the end
        // the offset (gst) is used in order to wrap values around from the end to the beginning and vice versa for the periodic case
        for (ix = 0; ix <= x_size - 1; ix++)
            PSD_t[gst + ix] = PSD[ix];

        // if periodic then wrap around the values (last gst values in the PSD go into the first gst spots of PSD_t)
        // add ghost points from boundary conditions
        switch(x_LBC_type)
        {
            case BoundaryConditionType::Periodic:
            // FYI: PSD_t[gst] = PSD[0] = PSD[x_size-1] = PSD_t[gst + x_size-1]; - this is by grid construction
                for (ig = 1; ig <= gst; ig++) {
                    PSD_t[gst - ig] = PSD[(x_size - 1) - ig];
                }
                break;
            // if constant value set the first gst spots to x_LBC
            case BoundaryConditionType::ConstantValue:
                for (ig = 1; ig <= gst; ig++) {
                    PSD_t[gst - ig] = x_LBC;
                }
                break;
            // take the derivative of each point with the previous one
            case BoundaryConditionType::ConstantDerivative:
                dx = (x[1] - x[0]);
                PSD_t[gst - 1] = PSD[0] - x_LBC * dx;
                for (ig = 1; ig <= gst; ig++) {
                    PSD_t[gst - ig] = PSD_t[gst - ig + 1] - x_LBC * dx;
                }
                break;
            default:
                Logger::error << "CONV_1D_BOUNDARY: unknown boundary type: " << x_LBC_type << std::endl;
        }

        // same thing for the upper boundary conditions

        switch(x_UBC_type)
        {
            case BoundaryConditionType::Periodic:
                // FYI: PSD_t[gst] = PSD[0] = PSD[x_size-1] = PSD_t[gst + x_size-1]; - this is by grid construction
                for (ig = 1; ig <= gst; ig++) {
                    PSD_t[gst + (x_size - 1) + ig] = PSD[ig];
                }
                break;
            case BoundaryConditionType::ConstantValue:
                for (ig = 1; ig <= gst; ig++) {
                    PSD_t[gst + (x_size - 1) + ig] = x_UBC;
                }
                break;
            case BoundaryConditionType::ConstantDerivative:
                dx = (x[x_size - 1] - x[x_size - 2]);
                PSD_t[gst + (x_size - 1) + 1] = PSD[x_size - 1] + x_UBC * dx;
                for (ig = 1; ig <= gst; ig++) {
                    PSD_t[gst + (x_size - 1) + ig] = PSD_t[gst + (x_size - 1) + ig - 1] + x_UBC * dx;
                }
                break;
            default:
                Logger::error << "CONV_1D_BOUNDARY: unknown boundary type: " << x_UBC_type << std::endl;
        }

        // All equations and formulas for these calculations can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf
        // Start calculation
        double psd_face_left, psd_face_right, courant_left, courant_right;
        for (ix = 0; ix <= x_size - 1; ix++) {
            // calculate at ix = 0 if the boundary condition is periodic
            if (ix == 0) {                           // special case
                if (x_LBC_type == BoundaryConditionType::Periodic) {  // Periodic
                    courant_right = (CourNum[x_size - 2] + CourNum[0]) / 2;
                    // CourNum_f = CourNum[0];
                } else {
                    continue;  // skip to the next ix - we need to calculate at ix = 1 only for periodic conditions
                }
            } else {
                courant_right = (CourNum[ix - 1] + CourNum[ix]) / 2;
                // CourNum_f = CourNum[ix];
            }
            
            const double courant_squared = courant_right * courant_right;
            const double courant_abs = fabs(courant_right);

            // Copy values into PSD_U, PSD_C, and PSD_D to use only these later
            if (courant_right >= 0) {  // Page 33
                // PSD_U = PSD[ix-2];
                // PSD_C = PSD[ix-1];
                // PSD_D = PSD[ix];

                for (ig = 1; ig <= uw_n; ig++)
                    PSD_U[ig] = PSD_t[gst + ix - 1 - ig];

                PSD_C = PSD_t[gst + ix - 1];

                for (ig = 1; ig <= dw_n; ig++)
                    PSD_D[ig] = PSD_t[gst + ix - 1 + ig];

            } else {
                // PSD_D = PSD[ix-1];
                // PSD_C = PSD[ix];
                // PSD_U = PSD[ix+1];

                for (ig = 1; ig <= dw_n; ig++)
                    PSD_D[ig] = PSD_t[gst + ix - ig];

                PSD_C = PSD_t[gst + ix];

                for (ig = 1; ig <= uw_n; ig++)
                    PSD_U[ig] = PSD_t[gst + ix + ig];
            }

            // From now on, only PSD_D, PSD_U, and PSD_C are in use

            // notes from Hayley: Because we only use PSD_D, PSD_U, and PSD_C, which we set based on the velocity direction - we now no longer need to consider the sign of the velocity - the absolute value is now used.

            // Using 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme
            // (2) Calculate some useful parameters
            DEL = PSD_D[1] - PSD_U[1];
            CURV = PSD_D[1] - 2 * PSD_C + PSD_U[1];
            ADEL = fabs(DEL);
            ACURV = fabs(CURV);

            CURVAV = 0.5 * fabs(PSD_U[1] - PSD_C - PSD_D[1] + PSD_D[2]);  // Leonard and Niknafs, 1991, (21)
            GRAD = fabs(PSD_C - PSD_D[1]);

            // Construct the upwind scheme
            //			if (ACURV >= ADEL) { // (4)
            //				// not monotionic
            //				PSD_f[ix] = PSD_C;
            //			} else
            // if (true) {
            // desired order scheme

            // PSD_f[ix] = 0.5 * (PSD_D[1] + PSD_C) - fabs(CourNum_f)/2 * (PSD_D[1] - PSD_C) - (1 - pow(CourNum_f,2))/6 * CURV; // Leonard and Niknafs, 1991, (85)
            //  note: CURV = PSD_D - 2*PSD_C + PSD_U;
            // continue;

            psd_face_right = 0.5 * (PSD_C + PSD_D[1]);  // A.1

            d1 = PSD_C - PSD_D[1];  // A.2
            // 2 order
            psd_face_right += 0.5 * courant_abs * d1;

            // Using table A1 to get the coefficients for the odd differences and table A2 for the even differences
            // Even is divided by 2
            D2 = (PSD_U[1] - PSD_C - PSD_D[1] + PSD_D[2]) / 2;    // A.3
            d3 = PSD_U[1] - 3 * PSD_C + 3 * PSD_D[1] - PSD_D[2];  //

            if (always_3rd_order) {  // !never_3rd_order && (always_3rd_order || (CURVAV < THC1 && GRAD < THG))) {
                // 3 order - QUICKEST method
                psd_face_right += (courant_squared - 1) / three_factorial * (D2 + d3 / 2);  // A.16

            } else {
                // 4 order - PSD_f is still second order
                psd_face_right += (courant_squared - 1) / three_factorial * (D2 + courant_abs / 4 * d3);

                D4 = (PSD_U[2] - 3 * PSD_U[1] + 2 * PSD_C + 2 * PSD_D[1] - 3 * PSD_D[2] + PSD_D[3]) / 2;
                d5 = PSD_U[2] - 5 * PSD_U[1] + 10 * PSD_C - 10 * PSD_D[1] + 5 * PSD_D[2] - PSD_D[3];

                if (always_5rd_order) {  // !never_5rd_order && always_5rd_order) {
                    // 5 order
                    psd_face_right += (courant_squared - 1) * (courant_squared - 4) / five_factorial * (D4 + d5 / 2);

                } else {
                    // 6 order - PSD_f is still 4th order.
                    psd_face_right += (courant_squared - 1) * (courant_squared - 4) / five_factorial * (D4 + courant_abs / 6 * d5);

                    D6 = (PSD_U[3] - 5 * PSD_U[2] + 9 * PSD_U[1] - 5 * PSD_C - 5 * PSD_D[1] + 9 * PSD_D[2] - 5 * PSD_D[3] + PSD_D[4]) / 2;
                    d7 = PSD_U[3] - 7 * PSD_U[2] + 21 * PSD_U[1] - 35 * PSD_C + 35 * PSD_D[1] - 21 * PSD_D[2] + 7 * PSD_D[3] - PSD_D[4];

                    if (always_7rd_order) {  // !never_7rd_order && (always_7rd_order || (CURVAV < THC2 && GRAD < THG))) {
                        // 7 order
                        psd_face_right += (courant_squared - 1) * (courant_squared - 4) * (courant_squared - 9) / seven_factorial * (D6 + d7 / 2);
                    } else {
                        // I think there is a bug somewhere here for 8-9 orders :( <- this comment is not from hayley
                        // 8 order
                        psd_face_right += (courant_squared - 1) * (courant_squared - 4) * (courant_squared - 9) / seven_factorial * (D6 + courant_abs / 8 * d7);

                        // 8 even:  +1 -7 +20 -28 +14  +14  -28 +20 -7 +1
                        // 9 odd: +1 -9 +36 -84 +126 -126 +84 -36 +9 -1
                        D8 = (PSD_U[4] - 7 * PSD_U[3] + 20 * PSD_U[2] - 28 * PSD_U[1] + 14 * PSD_C + 14 * PSD_D[1] - 28 * PSD_D[2] + 20 * PSD_D[3] - 7 * PSD_D[4] + PSD_D[5]) / 2;
                        d9 = PSD_U[4] - 9 * PSD_U[3] + 36 * PSD_U[2] - 84 * PSD_U[1] + 126 * PSD_C - 126 * PSD_D[1] + 84 * PSD_D[2] - 36 * PSD_D[3] + 9 * PSD_D[4] - PSD_D[5];
                        // 9 order
                        psd_face_right += (courant_squared - 1) * (courant_squared - 4) * (courant_squared - 9) * (courant_squared - 16) / nine_factorial * (D8 + d9 / 2);
                    }
                }
            }
            //}

            // Limit the result, so there is no 'waves' next to sharp gradients
            // also considering peaks - the peaks should remain to be peaks as much as possible
            // Leonidov and Niknafs, 1991, page 151
            //	dif1 = PSD_t[gst + ix - 2] - PSD_t[gst + ix - 3];
            //	dif2 = PSD_t[gst + ix - 1] - PSD_t[gst + ix - 2];
            //	dif3 = PSD_t[gst + ix - 0] - PSD_t[gst + ix - 1];
            //	dif4 = PSD_t[gst + ix + 1] - PSD_t[gst + ix + 0];
            //	dif5 = PSD_t[gst + ix + 2] - PSD_t[gst + ix + 1];
            //	dif6 = PSD_t[gst + ix + 3] - PSD_t[gst + ix + 2];

            // dif0 = PSD_D[3] - PSD_D[4];

            dif1 = PSD_D[2] - PSD_D[3];
            dif2 = PSD_D[1] - PSD_D[2];
            dif3 = PSD_C    - PSD_D[1];  //
            dif4 = PSD_U[1] - PSD_C;
            dif5 = PSD_U[2] - PSD_U[1];
            dif6 = PSD_U[3] - PSD_D[2];

            constexpr double treshhold = 0;                          // 1e-99;
            if (use_limiting && (!use_discreminator || !(   // if don't use the descriminator - always go inside
            // (((dif0 > treshhold && dif1 > treshhold) && (dif2 < -treshhold && dif3 < -treshhold)) && ((fabs(dif1) < fabs(dif0) - treshhold) && (fabs(dif2) < fabs(dif3) - treshhold))) ||
                (dif1 > treshhold && dif2 > treshhold && dif3 < -treshhold && dif4 < -treshhold && dif2 < dif1 - treshhold && dif3 > dif4 + treshhold) ||    
                (dif2 > treshhold && dif3 > treshhold && dif4 < -treshhold && dif5 < -treshhold && dif3 < dif2 - treshhold && dif4 > dif5 + treshhold) ||                    
                (dif3 > treshhold && dif4 > treshhold && dif5 < -treshhold && dif6 < -treshhold && dif4 < dif3 - treshhold && dif5 > dif6 + treshhold)
            ))) {
                //			if (!(
                //					(((dif1 > 0 && dif2 > 0) && (dif3 < 0 && dif4 < 0)) && ((fabs(dif2) < fabs(dif1)) && (fabs(dif3) < fabs(dif4)))) ||
                //					(((dif2 > 0 && dif3 > 0) && (dif4 < 0 && dif5 < 0)) && ((fabs(dif3) < fabs(dif2)) && (fabs(dif4) < fabs(dif5)))) ||
                //					(((dif3 > 0 && dif4 > 0) && (dif5 < 0 && dif6 < 0)) && ((fabs(dif4) < fabs(dif3)) && (fabs(dif5) < fabs(dif6))))
                //				)) {

                // limiting
                if (ACURV >= ADEL) {  // (4)
                    // not monotionic, or something
                    psd_face_right = PSD_C;
                } else if (courant_abs > 1e-31) {  // this will not work if velocity is zero
                    // cont. from page 48
                    // !!! PSD_REF = PSD_U + (PSD_C - PSD_U)/CourNum[ix]; // (3)
                    // !!! PSD_REF = PSD_U + (PSD_C - PSD_U)/CourNum_f; // (3)
                    PSD_REF = PSD_U[1] + (PSD_C - PSD_U[1]) / courant_abs;  // (3)

                    if (DEL > 0) {                                             // (5) limit
                        psd_face_right = fmax(psd_face_right, PSD_C);                    // max
                        psd_face_right = fmin(psd_face_right, fmin(PSD_REF, PSD_D[1]));  // min, min
                    } else if (DEL < 0) {                                      // (6)
                        psd_face_right = fmin(psd_face_right, PSD_C);                    // min
                        psd_face_right = fmax(psd_face_right, fmax(PSD_REF, PSD_D[1]));  // max, max
                    }
                }
            }

            // update PSD
            // courant_left, courant_right are values on left and right faces (between (i and i-1), and (i and i+1))
            if (ix > 0)  // special case
            {
                // (7)
                // PSD[ix] = PSD[ix] - CourNum[ix] * (PSD_r - PSD_l); - old version that makes the total PSD oscillate
                PSD[ix - 1] -= (courant_right * psd_face_right - courant_left * psd_face_left);
            }
            psd_face_left = psd_face_right;
            courant_left = courant_right;
        }
        if (x_LBC_type == BoundaryConditionType::Periodic) {  // special case
            PSD[x_size - 1] = PSD[0];        // need to update ix == x_size-1 point for periodic;
        }
        else if (x_LBC_type == BoundaryConditionType::ConstantValue) {
            PSD[0] = x_LBC;
        }
        else if (x_LBC_type == BoundaryConditionType::ConstantDerivative) {
            PSD[0] = PSD[1];
        }

        if (x_UBC_type == BoundaryConditionType::ConstantValue) {
            PSD[x_size - 1] = x_UBC;
        }
        else if (x_UBC_type == BoundaryConditionType::ConstantDerivative) {
            PSD[x_size - 1] = PSD[x_size - 2];
        }
    }

    return true;
}
