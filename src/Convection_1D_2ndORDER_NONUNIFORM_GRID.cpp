/**
 * \file Convection_1D_2ndODER_NONUNIFORM_GRID.cpp
 *
 * Leonard, 1997
 *
 * Leonard BP (1988) Universal Limiter for transient interpolation modeling of the advective transport equations: the ULTIMATE conservative difference scheme, NASA technical Memorandum 100916 ICOMP-88-11
 *
 */

#include "Convection_1D_2ndORDER_NONUNIFORM_GRID.h"
#include "BoundaryConditionType.hpp"

/** Using the Standard namespace
 */
using namespace std;

#define gst 5
#define uw_n 4
#define dw_n 4

bool Convection_1D_2ndORDER_NONUNIFORM_GRID(
        Matrix1D<double> &PSD, 
        const Matrix1D<double>& x,
        int x_size,
        double x_LBC, double x_UBC,
        BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
        const Matrix1D<double>& Ux,
        double dt_total)
{
    constexpr bool use_limiting      = false;
    constexpr bool use_discriminator = false;

    Matrix1D<double> CourNum(x_size);   // Courant number on faces
    Matrix1D<double> dx(x_size);  	    // different spacing in log-grid
    Matrix1D<double> UxOverDx(x_size); 	// velocity(iV)/dx(iV) on faces 

    // Face values
    for (int ix = 0; ix <= x_size-1; ix++) {
       if (ix == x_size-1) {
	       dx[ix] = (x[ix] - x[ix-1]);
       } else {
	       dx[ix] = (x[ix+1] - x[ix]);
       }
    }

    double maxCourNum = 1; 
    UxOverDx = Ux.divide(dx);   
    double ux_max = UxOverDx.maxabs();
    double steps_double = (double)dt_total * ux_max / maxCourNum;
    int num_steps = (steps_double <= 1) ? 1 : ceil(steps_double);
    double dt = dt_total / num_steps;
    CourNum = UxOverDx * dt;
    const double courant_max = CourNum.maxabs();
    if (courant_max > 1) {
        Logger::error << "max(CourNum) = " << courant_max << ", calculation can't be performed." << endl;
        exit(EXIT_FAILURE);
    }

    if (x_size < gst * 2 + 1) {
        Logger::error << "Number of convection grid can't be smaller than " << (gst * 2 + 1) << " due to numerical method used." << endl;
        // exit(EXIT_FAILURE);
    }

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
                Logger::error << "Constant derivative boundaries are not implemented yet" << std::endl;
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
                Logger::error << "Constant derivative boundaries are not implemented yet" << std::endl;
                break;
            default:
                Logger::error << "CONV_1D_BOUNDARY: unknown boundary type: " << x_UBC_type << std::endl;
        }

        // All equations and formulas for these calculations can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf
        // Start calculation
        double psd_face_left = 0.0, psd_face_right = 0.0, courant_right = 0.0;
        double velocity_right = 0.0, velocity_left = 0.0;
        double grad_FR = 0.0, grad_R = 0.0, curv_R = 0.0;  
        double dx_C = 0.0, dx_R = 0.0, dx_r = 0.0, dx_fr = 0.0; 

        grad_FR = 1.0 / dx[0] * (PSD_t[gst] - PSD_t[gst-1]); 
        dx_r  = dx[0]; 
        dx_fr = dx[1];
        dx_R  = 0.5*(dx_fr+dx_r);

        for (int ix = 0; ix <= (x_size-1); ix++) {
            // calculate at ix = 0 if the boundary condition is periodic
            if (ix == 0) {                           // special case
                if (x_LBC_type == BoundaryConditionType::Periodic) {  // Periodic
                    courant_right = (CourNum[x_size - 2] + CourNum[0]) / 2;
                    velocity_right = (Ux[x_size - 2] + Ux[0]) / 2;
                } 
                else {
                    courant_right  = (CourNum[ix+1] + CourNum[ix]) / 2;
                    velocity_right = (Ux[ix+1] + Ux[ix]) / 2; 
                //    continue;  // skip to the next ix - we need to calculate at ix = 1 only for periodic conditions
                }
            } else if ( ix == x_size-1 ) {
                if (x_LBC_type == BoundaryConditionType::Periodic) {  
                    courant_right  = (CourNum[x_size - 2] + CourNum[0]) / 2;
                    velocity_right = (Ux[x_size - 2] + Ux[0]) / 2;
                } else {
                    courant_right  = (CourNum[ix - 1] + CourNum[ix]) / 2;
                    velocity_right = (Ux[ix - 1] + Ux[ix]) / 2; 
                }
            } else {
                courant_right = (CourNum[ix + 1] + CourNum[ix]) / 2;
                velocity_right= (Ux[ix + 1] + Ux[ix]) / 2; 
            }
            
            // PSD down-wind, center, and up-wind
            double PSD_U[uw_n + 1], PSD_D[dw_n + 1], PSD_C;
            
            // Copy values into PSD_U, PSD_C, and PSD_D to use only these later
            if (courant_right >= 0) {  // Page 33
                for (ig = 1; ig <= uw_n; ig++)
                    PSD_U[ig] = PSD_t[gst + ix - 1 - ig];

                PSD_C = PSD_t[gst + ix - 1];

                for (ig = 1; ig <= dw_n; ig++)
                    PSD_D[ig] = PSD_t[gst + ix - 1 + ig];

                curv_R  = 1/dx_R * (grad_FR - grad_R);

                dx_C = dx_R;
                dx_r = dx_fr;
                dx_fr = (ix == x_size-1) ? dx[ix] : dx[ix + 1]; 
                dx_R = 0.5*(dx_fr + dx_r);
                grad_R = grad_FR;
                grad_FR = 1/dx_fr * (PSD_D[2] - PSD_D[1]);
            } else {
                for (ig = 1; ig <= dw_n; ig++)
                    PSD_D[ig] = PSD_t[gst + ix - ig];

                PSD_C = PSD_t[gst + ix];

                for (ig = 1; ig <= uw_n; ig++)
                    PSD_U[ig] = PSD_t[gst + ix + ig];

                dx_C = dx_R;
                dx_r = dx_fr;
                dx_fr = (ix == x_size-1) ? dx[ix] : dx[ix + 1]; 
                dx_R = 0.5*(dx_fr + dx_r);

                grad_R = grad_FR;
                grad_FR = 1/dx_fr * (PSD_U[1] - PSD_C);

                curv_R  = 1/dx_R * (grad_FR - grad_R);
            }

            // ~ 1st order upwind scheme
            // psd_face_right = PSD_C;
            // ~ ULTIMATE QUICKEST
            psd_face_right = 0.5*(PSD_C + PSD_D[1]) - 0.5*courant_right*grad_R*dx_r - 0.5 * pow(dx_r,2)*1.0/3*(1-pow(courant_right,2))*curv_R; 

            // update PSD_out
            // courant_left, courant_right are values on left and right faces (between (i and i-1), and (i and i+1))
            if (ix > 0)  // special case
            {
                PSD_out[ix - 1] = PSD[ix - 1] - dt / dx_C * (velocity_right * psd_face_right - velocity_left * psd_face_left);
            }
            if(use_limiting)
            {
                PSD_faces[ix] = psd_face_right;
            }
            psd_face_left = psd_face_right;
            velocity_left= velocity_right;

        }


        if (x_LBC_type == BoundaryConditionType::Periodic) {  // special case
            PSD_out[x_size - 1] = PSD_out[0];        // need to update ix == x_size-1 point for periodic;
        }
        else if (x_LBC_type == BoundaryConditionType::ConstantValue) {
            PSD_out[0] = x_LBC;
        }
        else if (x_LBC_type == BoundaryConditionType::ConstantDerivative) {
            Logger::error << "Constant derivative boundaries are not implemented yet" << std::endl;
            break;
        }

        if (x_UBC_type == BoundaryConditionType::ConstantValue) {
            PSD_out[x_size - 1] = x_UBC;
        }
        else if (x_UBC_type == BoundaryConditionType::ConstantDerivative) {
            Logger::error << "Constant derivative boundaries are not implemented yet" << std::endl;
            break;
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
                        //////PSD_unlimited_t[gst - ig] = PSD_unlimited[(x_size - 1) - ig];
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
                    Logger::error << "Constant derivative boundaries are not implemented yet" << std::endl;
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
                    Logger::error << "Constant derivative boundaries are not implemented yet" << std::endl;
                    break;
                default:
                    Logger::error << "CONV_1D_BOUNDARY: unknown boundary type: " << x_UBC_type << std::endl;
            }
            double psd_face_left = 0.0, psd_face_right = 0.0, courant_right = 0.0;
            double velocity_right = 0.0, velocity_left = 0.0;
            double grad_FR = 0.0, grad_R = 0.0, curv_R = 0.0;  
            double dx_C = 0.0, dx_R = 0.0, dx_l, dx_r = 0.0, dx_fr = 0.0; 
            double dx_ffl = 0.0, dx_fl = 0.0, dx_ffr = 0.0; 

            grad_FR = 1.0 / dx[0] * (PSD_unlimited_t[gst] - PSD_unlimited_t[gst-1]); 
            dx_ffl= dx[0]; 
            dx_fl = dx[0]; 
            dx_l  = dx[0]; 
            dx_r  = dx[0]; 
            dx_fr = dx[1];
            dx_ffr= dx[2];
            dx_R  = 0.5*(dx_fr+dx_r);

            for (int ix = 0; ix <= x_size - 1; ix++)
            {
                if (ix == 0) {                           // special case
                    if (x_LBC_type == BoundaryConditionType::Periodic) {  // Periodic
                        courant_right = (CourNum[x_size - 2] + CourNum[0]) / 2;
                        velocity_right= (Ux[x_size - 2] + Ux[0]) / 2;
                    } else {
                        continue;  // skip to the next ix - we need to calculate at ix = 1 only for periodic conditions
                    }
                } else if ( ix == x_size-1 ) {
                    if (x_LBC_type == BoundaryConditionType::Periodic) {  
                        courant_right  = (CourNum[x_size - 2] + CourNum[0]) / 2;
                        velocity_right = (Ux[x_size - 2] + Ux[0]) / 2;
                    } else {
                        courant_right  = (CourNum[ix - 1] + CourNum[ix]) / 2;
                        velocity_right = (Ux[ix - 1] + Ux[ix]) / 2; 
                    }
                } else {
                    courant_right = (CourNum[ix + 1] + CourNum[ix]) / 2;
                    velocity_right= (Ux[ix + 1] + Ux[ix]) / 2; 
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

                    curv_R  = 1/dx_R * (grad_FR - grad_R);

                    dx_C = dx_R;
                    dx_ffl= dx_fl;
                    dx_fl = dx_l; 
                    dx_l  = dx_r; 
                    dx_r  = dx_fr;
                    dx_fr = (ix == x_size-1) ? dx[ix] : dx[ix + 1]; 
                    dx_ffr= (ix <= x_size-3) ? dx[ix] : dx[ix + 2]; 
                    dx_R = 0.5*(dx_fr + dx_r);
                    grad_R = grad_FR;
                    grad_FR = 1/dx_fr * (PSD_unlimited_D[2] - PSD_unlimited_D[1]);
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

                    dx_C = dx_R;
                    dx_ffl= dx_fl;
                    dx_fl = dx_l; 
                    dx_l = dx_r;
                    dx_r = dx_fr;
                    dx_fr = (ix == x_size-1) ? dx[ix] : dx[ix + 1]; 
                    dx_ffr= (ix <= x_size-3) ? dx[ix] : dx[ix + 2]; 
                    dx_R = 0.5*(dx_fr + dx_r);

                    grad_R = grad_FR;
                    grad_FR = 1/dx_fr * (PSD_unlimited_U[1] - PSD_unlimited_C);

                    curv_R  = 1/dx_R * (grad_FR - grad_R);
                }
            
                // Using 3.5 Simplified Ultimate Quickest strategy from B.P Leonard the Ultimate conservative difference scheme
                // (2) Calculate some useful parameters
                const double DEL   = (orig_downwind - orig_upwind) / (dx_l + dx_r); 
                const double CURV  = curv_R; 
                const double ADEL  = std::fabs(DEL);
                const double ACURV = std::fabs(CURV);
                
                // Leonard and Niknafs, 1991, (21)
                // const double CURVAV = 0.5 * std::fabs(PSD_U[1] - PSD_C - PSD_D[1] + PSD_D[2]);
                // const double GRAD   = std::fabs(PSD_C - PSD_D[1]);

                // the updated unlimited values are used to decide wether or not an unphysical oscillation occured
                const double dif1 = (PSD_unlimited_D[2] - PSD_unlimited_D[3]) / dx_ffr;
                const double dif2 = (PSD_unlimited_D[1] - PSD_unlimited_D[2]) / dx_fr;
                const double dif3 = (PSD_unlimited_C    - PSD_unlimited_D[1]) / dx_r;  
                const double dif4 = (PSD_unlimited_U[1] - PSD_unlimited_C) / dx_l;
                const double dif5 = (PSD_unlimited_U[2] - PSD_unlimited_U[1]) / dx_fl;
                const double dif6 = (PSD_unlimited_U[3] - PSD_unlimited_U[2]) / dx_ffl;
                
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
                    PSD[ix - 1] -= dt / dx_C * (velocity_right * psd_face_right - velocity_left * psd_face_left);
                }
                psd_face_left = psd_face_right;
                velocity_left = velocity_right;
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
                Logger::error << "Constant derivative boundaries are not implemented yet" << std::endl;
                break;
            }

            if (x_UBC_type == BoundaryConditionType::ConstantValue)
            {
                PSD[x_size - 1] = x_UBC;
            }
            else if (x_UBC_type == BoundaryConditionType::ConstantDerivative)
            {
                Logger::error << "Constant derivative boundaries are not implemented yet" << std::endl;
                break;
            }
        }
    }

    return true;
}
