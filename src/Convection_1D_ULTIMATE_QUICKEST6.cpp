// SPDX-FileCopyrightText: 2015 UCLA
// SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * \file Convection_1D_ULTIMATE_QUICKEST6.cpp
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
		double dt_total)
{
    constexpr bool use_limiting      = true;
    constexpr bool use_discriminator = true;
    constexpr bool always_3rd_order  = false;
    constexpr bool always_5rd_order  = false;
    constexpr bool always_7rd_order  = false;
    
    Matrix1D<double> CourNum(x_size);  // Courant number
    double dx = (x[1] - x[0]);  // FIXME - will work only for regular (uniform) grid

    double maxCourNum   = 1;
    double ux_max       = Ux.maxabs();
    double steps_double = (double)dt_total * ux_max / dx / maxCourNum;
    int num_steps       = (steps_double <= 1) ? 1 : (int)ceil(steps_double);

    double dt = dt_total / num_steps;

    CourNum = Ux / dx * dt;
    const double courant_max = std::abs(ux_max / dx * dt);
    if (courant_max > 1) {
        Logger::error << "max(CourNum) = " << courant_max << ", calculation can't be performed." << endl;
        exit(EXIT_FAILURE);
    }

    // double THC1 = 1e6, THC2 = 1e12, THG = 1e6;        // Leonard and Niknfas, 1991

    if (x_size < gst * 2 + 1) {
        Logger::error << "Number of convection grid can't be smaller than " << (gst * 2 + 1) << " due to numerical method used." << endl;
        // exit(EXIT_FAILURE);
    }

    // bool never_3rd_order = false
    // bool never_5rd_order = false
    // bool never_7rd_order = false

    constexpr double three_factorial = 3 * 2 * 1;
    constexpr double five_factorial  = 5 * 4 * three_factorial;
    constexpr double seven_factorial = 7 * 6 * five_factorial;
    constexpr double nine_factorial  = 9 * 8 * seven_factorial;

    // Temporary PSD to use during calculations, includes ghost points
    Matrix1D<double> PSD_t(x_size + gst + gst);
    
    Matrix1D<double> PSD_unlimited;
    Matrix1D<double> PSD_unlimited_t;
    Matrix1D<double> PSD_faces;
    if(use_limiting)
    {
        PSD_unlimited.AllocateMemory(x_size);
        PSD_unlimited_t.AllocateMemory(x_size + gst + gst);
        PSD_faces.AllocateMemory(x_size);
    }
    // if no limiter is used, write directly to PSD
    Matrix1D<double>& PSD_out = use_limiting ? PSD_unlimited : PSD;

    int ig;  // ghost points iterator
    for (int it = 0; it < num_steps; it++) {
        // copying PSD into PSD_t offset by gst, with gst extra spots at the end
        // the offset (gst) is used in order to wrap values around from the end to the beginning and vice versa for the periodic case
        for (int ix = 0; ix <= x_size - 1; ix++)
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
        double psd_face_left = 0.0, psd_face_right = 0.0, courant_left = 0.0, courant_right = 0.0;
        for (int ix = 0; ix <= x_size - 1; ix++) {
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
            const double courant_abs = std::fabs(courant_right);

            // PSD down-wind, center, and up-wind
            double PSD_U[uw_n + 1], PSD_D[dw_n + 1], PSD_C;
            
            // Copy values into PSD_U, PSD_C, and PSD_D to use only these later
            if (courant_right >= 0) {  // Page 33
                for (ig = 1; ig <= uw_n; ig++)
                    PSD_U[ig] = PSD_t[gst + ix - 1 - ig];

                PSD_C = PSD_t[gst + ix - 1];

                for (ig = 1; ig <= dw_n; ig++)
                    PSD_D[ig] = PSD_t[gst + ix - 1 + ig];

            } else {
                for (ig = 1; ig <= dw_n; ig++)
                    PSD_D[ig] = PSD_t[gst + ix - ig];

                PSD_C = PSD_t[gst + ix];

                for (ig = 1; ig <= uw_n; ig++)
                    PSD_U[ig] = PSD_t[gst + ix + ig];
            }

            // From now on, only PSD_D, PSD_U, and PSD_C are in use
            // notes from Hayley: Because we only use PSD_D, PSD_U, and PSD_C, which we set based on the velocity direction
            // we now no longer need to consider the sign of the velocity - the absolute value is now used.

            psd_face_right = 0.5 * (PSD_C + PSD_D[1]);  // A.1

            const double d1 = PSD_C - PSD_D[1];  // A.2
            // 2 order
            psd_face_right += 0.5 * courant_abs * d1;

            // Using table A1 to get the coefficients for the odd differences and table A2 for the even differences
            // Even is divided by 2
            const double D2 = (PSD_U[1] - PSD_C - PSD_D[1] + PSD_D[2]) / 2;    // A.3
            const double d3 = PSD_U[1] - 3 * PSD_C + 3 * PSD_D[1] - PSD_D[2];  //

            if (always_3rd_order) {  // !never_3rd_order && (always_3rd_order || (CURVAV < THC1 && GRAD < THG))) {
                // 3 order - QUICKEST method
                psd_face_right += (courant_squared - 1) / three_factorial * (D2 + d3 / 2);  // A.16

            } else {
                // 4 order - PSD_f is still second order
                psd_face_right += (courant_squared - 1) / three_factorial * (D2 + courant_abs / 4 * d3);

                const double D4 = (PSD_U[2] - 3 * PSD_U[1] + 2 * PSD_C + 2 * PSD_D[1] - 3 * PSD_D[2] + PSD_D[3]) / 2;
                const double d5 = PSD_U[2] - 5 * PSD_U[1] + 10 * PSD_C - 10 * PSD_D[1] + 5 * PSD_D[2] - PSD_D[3];

                if (always_5rd_order) {  // !never_5rd_order && always_5rd_order) {
                    // 5 order
                    psd_face_right += (courant_squared - 1) * (courant_squared - 4) / five_factorial * (D4 + d5 / 2);

                } else {
                    // 6 order - PSD_f is still 4th order.
                    psd_face_right += (courant_squared - 1) * (courant_squared - 4) / five_factorial * (D4 + courant_abs / 6 * d5);

                    const double D6 = (PSD_U[3] - 5 * PSD_U[2] + 9 * PSD_U[1] - 5 * PSD_C - 5 * PSD_D[1] + 9 * PSD_D[2] - 5 * PSD_D[3] + PSD_D[4]) / 2;
                    const double d7 = PSD_U[3] - 7 * PSD_U[2] + 21 * PSD_U[1] - 35 * PSD_C + 35 * PSD_D[1] - 21 * PSD_D[2] + 7 * PSD_D[3] - PSD_D[4];

                    if (always_7rd_order) {  // !never_7rd_order && (always_7rd_order || (CURVAV < THC2 && GRAD < THG))) {
                        // 7 order
                        psd_face_right += (courant_squared - 1) * (courant_squared - 4) * (courant_squared - 9) / seven_factorial * (D6 + d7 / 2);
                    } else {
                        // I think there is a bug somewhere here for 8-9 orders :( <- this comment is not from hayley
                        // 8 order
                        psd_face_right += (courant_squared - 1) * (courant_squared - 4) * (courant_squared - 9) / seven_factorial * (D6 + courant_abs / 8 * d7);

                        // 8 even:  +1 -7 +20 -28 +14  +14  -28 +20 -7 +1
                        // 9 odd: +1 -9 +36 -84 +126 -126 +84 -36 +9 -1
                        const double D8 = (PSD_U[4] - 7 * PSD_U[3] + 20 * PSD_U[2] - 28 * PSD_U[1] + 14 * PSD_C + 14 * PSD_D[1] - 28 * PSD_D[2] + 20 * PSD_D[3] - 7 * PSD_D[4] + PSD_D[5]) / 2;
                        const double d9 = PSD_U[4] - 9 * PSD_U[3] + 36 * PSD_U[2] - 84 * PSD_U[1] + 126 * PSD_C - 126 * PSD_D[1] + 84 * PSD_D[2] - 36 * PSD_D[3] + 9 * PSD_D[4] - PSD_D[5];
                        // 9 order
                        psd_face_right += (courant_squared - 1) * (courant_squared - 4) * (courant_squared - 9) * (courant_squared - 16) / nine_factorial * (D8 + d9 / 2);
                    }
                }
            }

            // update PSD_out
            // courant_left, courant_right are values on left and right faces (between (i and i-1), and (i and i+1))
            if (ix > 0)  // special case
            {
                // (7)
                PSD_out[ix-1] = PSD[ix - 1] - (courant_right * psd_face_right - courant_left * psd_face_left);
            }
            if(use_limiting)
            {
                PSD_faces[ix] = psd_face_right;
            }
            psd_face_left = psd_face_right;
            courant_left = courant_right;
        }
        if (x_LBC_type == BoundaryConditionType::Periodic) {  // special case
            PSD_out[x_size - 1] = PSD_out[0];        // need to update ix == x_size-1 point for periodic;
        }
        else if (x_LBC_type == BoundaryConditionType::ConstantValue) {
            PSD_out[0] = x_LBC;
        }
        else if (x_LBC_type == BoundaryConditionType::ConstantDerivative) {
            PSD_out[0] = PSD_out[1];
        }

        if (x_UBC_type == BoundaryConditionType::ConstantValue) {
            PSD_out[x_size - 1] = x_UBC;
        }
        else if (x_UBC_type == BoundaryConditionType::ConstantDerivative) {
            PSD_out[x_size - 1] = PSD_out[x_size - 2];
        }
        
        if(use_limiting)
        {   
            // copy the unlimited result of the upwind method into PSD_t and set the ghost point
            for (int ix = 0; ix <= x_size - 1; ix++)
            {
                PSD_unlimited_t[gst + ix] = PSD_unlimited[ix];
            }
            switch(x_LBC_type)
            {
                case BoundaryConditionType::Periodic:
                // FYI: PSD_t[gst] = PSD[0] = PSD[x_size-1] = PSD_t[gst + x_size-1]; - this is by grid construction
                    for (ig = 1; ig <= gst; ig++) {
                        PSD_unlimited_t[gst - ig] = PSD[(x_size - 1) - ig];
                    }
                    break;
                // if constant value set the first gst spots to x_LBC
                case BoundaryConditionType::ConstantValue:
                    for (ig = 1; ig <= gst; ig++) {
                        PSD_unlimited_t[gst - ig] = x_LBC;
                    }
                    break;
                // take the derivative of each point with the previous one
                case BoundaryConditionType::ConstantDerivative:
                    dx = (x[1] - x[0]);
                    PSD_unlimited_t[gst - 1] = PSD[0] - x_LBC * dx;
                    for (ig = 1; ig <= gst; ig++) {
                        PSD_unlimited_t[gst - ig] = PSD_unlimited_t[gst - ig + 1] - x_LBC * dx;
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
                        PSD_unlimited_t[gst + (x_size - 1) + ig] = PSD_unlimited[ig];
                    }
                    break;
                case BoundaryConditionType::ConstantValue:
                    for (ig = 1; ig <= gst; ig++) {
                        PSD_unlimited_t[gst + (x_size - 1) + ig] = x_UBC;
                    }
                    break;
                case BoundaryConditionType::ConstantDerivative:
                    dx = (x[x_size - 1] - x[x_size - 2]);
                    PSD_unlimited_t[gst + (x_size - 1) + 1] = PSD_unlimited[x_size - 1] + x_UBC * dx;
                    for (ig = 1; ig <= gst; ig++) {
                        PSD_unlimited_t[gst + (x_size - 1) + ig] = PSD_unlimited_t[gst + (x_size - 1) + ig - 1] + x_UBC * dx;
                    }
                    break;
                default:
                    Logger::error << "CONV_1D_BOUNDARY: unknown boundary type: " << x_UBC_type << std::endl;
            }
            
            //double psd_face_left = 0.0, courant_left = 0.0, courant_right = 0.0; // Alreayd defined
            psd_face_left = 0.0, courant_left = 0.0, courant_right = 0.0;
            for (int ix = 0; ix <= x_size - 1; ix++) 
            {
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
                
                const double courant_abs = std::fabs(courant_right);
                
                // PSD down-wind, center, and up-wind
                double orig_upwind;
                double orig_downwind;
                double orig_centre;
                double PSD_unlimited_U[uw_n + 1];
                double PSD_unlimited_D[dw_n + 1];
                double PSD_unlimited_C;
                
                // Copy the original values of PSD into upwind, centre and downwind for the limiter
                // and the updated unlimited values of PSD for the discriminator 
                if (courant_right >= 0) {  // Page 33
                    for (ig = 1; ig <= uw_n; ig++)
                    {
                        PSD_unlimited_U[ig] = PSD_unlimited_t[gst + ix - 1 - ig];
                    }
                    orig_upwind = PSD_t[gst + ix - 2];

                    PSD_unlimited_C = PSD_unlimited_t[gst + ix - 1];
                    orig_centre = PSD_t[gst + ix - 1];

                    for (ig = 1; ig <= dw_n; ig++)
                    {
                        PSD_unlimited_D[ig] = PSD_unlimited_t[gst + ix - 1 + ig];
                    }
                    orig_downwind = PSD_t[gst + ix];

                } else {
                    for (ig = 1; ig <= dw_n; ig++)
                    {
                        PSD_unlimited_D[ig] = PSD_unlimited_t[gst + ix - ig];
                    }
                    orig_downwind = PSD_t[gst + ix - 1];

                    PSD_unlimited_C = PSD_unlimited_t[gst + ix];
                    orig_centre = PSD_t[gst + ix];

                    for (ig = 1; ig <= uw_n; ig++)
                    {
                        PSD_unlimited_U[ig] = PSD_unlimited_t[gst + ix + ig];
                    }
                    orig_upwind = PSD_t[gst + ix + 1];
                }
            
                // Using 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme
                // (2) Calculate some useful parameters
                const double DEL   = orig_downwind - orig_upwind;
                const double CURV  = orig_downwind - 2 * orig_centre + orig_upwind;
                const double ADEL  = std::fabs(DEL);
                const double ACURV = std::fabs(CURV);
                
                // Leonard and Niknafs, 1991, (21)
                // const double CURVAV = 0.5 * std::fabs(PSD_U[1] - PSD_C - PSD_D[1] + PSD_D[2]);
                // const double GRAD   = std::fabs(PSD_C - PSD_D[1]);

                // the updated unlimited values are used to decide wether or not an unphysical oscillation occured
                const double dif1 = PSD_unlimited_D[2] - PSD_unlimited_D[3];
                const double dif2 = PSD_unlimited_D[1] - PSD_unlimited_D[2];
                const double dif3 = PSD_unlimited_C    - PSD_unlimited_D[1];
                const double dif4 = PSD_unlimited_U[1] - PSD_unlimited_C;
                const double dif5 = PSD_unlimited_U[2] - PSD_unlimited_U[1];
                const double dif6 = PSD_unlimited_U[3] - PSD_unlimited_U[2];
                
                constexpr double treshhold = 0; // 1e-99;

                // the discriminator - is there a realistic peak at one of the three points upwind, centre, downwind?
                const bool peak_downwind = 
                    dif1 > treshhold && dif2 > treshhold && dif3 < -treshhold &&
                    dif4 < -treshhold && dif2 < dif1 - treshhold && dif3 > dif4 + treshhold;

                const bool peak_centre =
                    dif2 > treshhold && dif3 > treshhold && dif4 < -treshhold &&
                    dif5 < -treshhold && dif3 < dif2 - treshhold && dif4 > dif5 + treshhold;

                const bool peak_upwind =
                    dif3 > treshhold && dif4 > treshhold && dif5 < -treshhold &&
                    dif6 < -treshhold && dif4 < dif3 - treshhold && dif5 > dif6 + treshhold;

                // Already defined
                // double psd_face_right = PSD_faces[ix];
                psd_face_right = PSD_faces[ix];
                
                // use the limiter only if there is no realistic peak or if the discriminator is switched off
                if (!use_discriminator || !(peak_downwind || peak_centre || peak_upwind))
                {
                    // limiting
                    if (ACURV >= ADEL) {  // (4)
                        // non monotionic case 
                        psd_face_right = orig_centre;
                    } else if (courant_abs > 1e-31) {  // this will not work if velocity is zero
                        // cont. from page 48
                        const double PSD_REF = orig_upwind + (orig_centre - orig_upwind) / courant_abs;  // (3)

                        if (DEL > 0) {                                             // (5) limit
                            psd_face_right = std::fmax(psd_face_right, orig_centre);                    // max
                            psd_face_right = std::fmin(psd_face_right, std::fmin(PSD_REF, orig_downwind));  // min, min
                        } else if (DEL < 0) {                                      // (6)
                            psd_face_right = std::fmin(psd_face_right, orig_centre);                    // min
                            psd_face_right = std::fmax(psd_face_right, std::fmax(PSD_REF, orig_downwind));  // max, max
                        }
                    }
                }
                if (ix > 0)
                {
                    // use the limited psd_face_right to compute the final output
                    PSD[ix - 1] -= (courant_right * psd_face_right - courant_left * psd_face_left);
                }
                psd_face_left = psd_face_right;
                courant_left  = courant_right;
            }
            // apply boundary conditions on the results of the limiter
            if (x_LBC_type == BoundaryConditionType::Periodic) 
            {
                PSD[x_size - 1] = PSD[0];        // need to update ix == x_size-1 point for periodic;
            }
            else if (x_LBC_type == BoundaryConditionType::ConstantValue)
            {
                PSD[0] = x_LBC;
            }
            else if (x_LBC_type == BoundaryConditionType::ConstantDerivative)
            {
                PSD[0] = PSD[1];
            }

            if (x_UBC_type == BoundaryConditionType::ConstantValue)
            {
                PSD[x_size - 1] = x_UBC;
            }
            else if (x_UBC_type == BoundaryConditionType::ConstantDerivative)
            {
                PSD[x_size - 1] = PSD[x_size - 2];
            }
        }
    }

    return true;
}
