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

/** Using the Standard namespace
*/
using namespace std;

#define gst 5
#define uw_n 5
#define dw_n 4

bool Convection_1D_ULTIMATE_QUICKEST6( 
		Matrix1D<double> &PSD, 
		const Matrix1D<double>& x,
		int x_size,
		double x_LBC, double x_UBC,
		string x_LBC_type, string x_UBC_type,
		const Matrix1D<double>& Ux,
		const Matrix1D<double>& Sources, const Matrix1D<double>& Losses,
		double dt_total) 
{
	int ix;

	double DEL, CURV, ADEL, ACURV, PSD_REF, CURVAV, GRAD; // temporary variables
	Matrix1D < double > PSD_f(x_size); // PSD on a 'face' (between nodes)

	Matrix1D < double > CourNum(x_size); // Courant number

	double dx = (x[1] - x[0]); // FIXME - will work only for regular (uniform) grid
	//double max_dt = dx / Ux.maxabs();
	//int num_steps = (max_dt < dt_total) ? ceil((double)dt_total/max_dt) : 1;

	double maxCourNum = 1;
	int num_steps = (((double) dt_total * Ux.maxabs() / dx / maxCourNum) <= 1) ? 1 : ceil((double) dt_total * Ux.maxabs() / dx / maxCourNum);

	double dt = dt_total / num_steps;

	CourNum = Ux/dx*dt;
	if (CourNum.maxabs() > 1) {
		Logger::error << "max(CourNum) = " << CourNum.maxabs() << ", calculation can't be performed." << endl;
		exit(EXIT_FAILURE);
	}

	//double CourNum_f; // Courant number of a face
	//added by hayley for testing
	Matrix1D < double > CourNum_f(x_size);


	double d1, d3, d5, d7, d9;
	double D2, D4, D6, D8, D10;
	double dif0, dif1, dif2, dif3, dif4, dif5, dif6; // Leonard and Niknfas, 1991
	double THC1 = 1e6, THC2 = 1e12, THG = 1e6; // Leonard and Niknfas, 1991

	// PSD down-wind, center, and up-wind
	// double PSD_U[uw_n], PSD_D[dw_n], PSD_C;
	double PSD_U[uw_n+1], PSD_D[dw_n+1], PSD_C;

	if (x_size < gst*2+1) {
		Logger::error << "Number of convection grid can't be smaller than " << (gst * 2 + 1) << " due to numerical method used." << endl;
		//exit(EXIT_FAILURE);
	}

	bool use_discreminator = true, use_limiting = true;
	bool never_3rd_order = false, always_3rd_order = false;
	bool never_5rd_order = false, always_5rd_order = false;
	bool never_7rd_order = false, always_7rd_order = false;

	// Temporary PSD to use during calculations, includes ghost points
	Matrix1D < double > PSD_t(x_size + gst + gst);

	int ig; // ghost points iterator
	for (int it = 0; it < num_steps; it++) {
		// copying PSD into PSD_t offset by gst, with gst extra spots at the end
		// the offset (gst) is used in order to wrap values around from the end to the beginning and vice versa for the periodic case
		for (ix = 0; ix <= x_size-1; ix++)
			PSD_t[gst + ix] = PSD[ix];


		// if periodic then wrap around the values (last gst values in the PSD go into the first gst spots of PSD_t)
		// add ghost points from boundary conditions
		if (x_LBC_type == "BCT_PERIODIC") { // Periodic
			//FYI: PSD_t[gst] = PSD[0] = PSD[x_size-1] = PSD_t[gst + x_size-1]; - this is by grid construction
			for (ig = 1; ig <= gst; ig++)
				PSD_t[gst - ig] = PSD[(x_size-1) - ig];

		// if constant value set the first gst spots to x_LBC
		} else if (x_LBC_type == "BCT_CONSTANT_VALUE") { // for condition on value
			for (ig = 1; ig <= gst; ig++)
				PSD_t[gst - ig] = x_LBC;

		// take the derivative of each point with the previous one
		} else if (x_LBC_type == "BCT_CONSTANT_DERIVATIVE") { // for condition on derivative
			dx = (x[1] - x[0]);
			PSD_t[gst - 1] = PSD[0] - x_LBC * dx;
			for (ig = 1; ig <= gst; ig++)
				PSD_t[gst - ig] = PSD_t[gst - ig + 1] - x_LBC * dx;

		} else {
			Logger::error << "2D_DIFF_BOUNDARY: unknown boundary type: " << x_LBC_type.c_str() << endl;
			//exit(EXIT_FAILURE);
		}

		// same thing for the upper boundary conditions

		if (x_UBC_type == "BCT_PERIODIC") { // Periodic
			//FYI: PSD_t[gst] = PSD[0] = PSD[x_size-1] = PSD_t[gst + x_size-1]; - this is by grid construction
			for (ig = 1; ig <= gst; ig++)
				PSD_t[gst + (x_size-1) + ig] = PSD[ig];

		} else if (x_UBC_type == "BCT_CONSTANT_VALUE") { // for condition on value
			for (ig = 1; ig <= gst; ig++)
				PSD_t[gst + (x_size-1) + ig] = x_UBC;

		} else if (x_UBC_type == "BCT_CONSTANT_DERIVATIVE") { // for condition on derivative
			dx = (x[x_size-1] - x[x_size-2]);
			PSD_t[gst + (x_size-1) + 1] = PSD[x_size-1] + x_UBC * dx;
			for (ig = 1; ig <= gst; ig++)
				PSD_t[gst + (x_size-1) + ig] = PSD_t[gst + (x_size-1) + ig - 1] + x_UBC * dx;

		} else {
			Logger::error << "2D_DIFF_BOUNDARY: unknown boundary type: " << x_UBC_type.c_str() << endl;
			//exit(EXIT_FAILURE);
		}

		// All equations and formulas for these calculations can be found at http://www.hadian.ir/teaching/CompHydr/3.pdf
		// Start calculation
		for (ix = 0; ix <= x_size-1; ix++) {

			// calculate at ix = 0 if the boundary condition is periodic
			if (ix == 0) { // special case
				if (x_LBC_type == "BCT_PERIODIC") { // Periodic
					CourNum_f[ix] = (CourNum[x_size-2] + CourNum[0])/2;
					//CourNum_f = CourNum[0];
				} else {
					continue; // skip to the next ix - we need to calculate at ix = 1 only for periodic conditions
				}
			} else {
				CourNum_f[ix] = (CourNum[ix-1] + CourNum[ix])/2;
				//CourNum_f = CourNum[ix];
			}

			// Copy values into PSD_U, PSD_C, and PSD_D to use only these later
			if (CourNum_f[ix] >= 0) { // Page 33
				//PSD_U = PSD[ix-2];
				//PSD_C = PSD[ix-1];
				//PSD_D = PSD[ix];

				for (ig = 1; ig <= uw_n; ig++)
					PSD_U[ig] = PSD_t[gst + ix - 1 - ig];

				PSD_C = PSD_t[gst + ix - 1];

				for (ig = 1; ig <= dw_n; ig++)
					PSD_D[ig] = PSD_t[gst + ix - 1 + ig];

			} else {
				//PSD_D = PSD[ix-1];
				//PSD_C = PSD[ix];
				//PSD_U = PSD[ix+1];

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
			CURV = PSD_D[1] - 2*PSD_C + PSD_U[1];
			ADEL = fabs(DEL);
			ACURV = fabs(CURV);

			CURVAV = 0.5 * fabs(PSD_U[1] - PSD_C - PSD_D[1] + PSD_D[2]); // Leonard and Niknafs, 1991, (21)
			GRAD = fabs(PSD_C - PSD_D[1]);

			// Construct the upwind scheme
//			if (ACURV >= ADEL) { // (4)
//				// not monotionic
//				PSD_f[ix] = PSD_C;
//			} else
			//if (true) {
				// desired order scheme
				
				//PSD_f[ix] = 0.5 * (PSD_D[1] + PSD_C) - fabs(CourNum_f)/2 * (PSD_D[1] - PSD_C) - (1 - pow(CourNum_f,2))/6 * CURV; // Leonard and Niknafs, 1991, (85)
				// note: CURV = PSD_D - 2*PSD_C + PSD_U;
				//continue;

				PSD_f[ix] = 0.5 * (PSD_C + PSD_D[1]); // A.1

				d1 = PSD_C - PSD_D[1]; //A.2
				// 2 order
				PSD_f[ix] = PSD_f[ix] + 0.5 * fabs(CourNum_f[ix]) * d1;


				// Using table A1 to get the coefficients for the odd differences and table A2 for the even differences
				// Even is divided by 2
				D2 = (PSD_U[1] - PSD_C - PSD_D[1] + PSD_D[2])/2; // A.3
				d3 = PSD_U[1] - 3*PSD_C + 3*PSD_D[1] - PSD_D[2]; //

				if (always_3rd_order) { // !never_3rd_order && (always_3rd_order || (CURVAV < THC1 && GRAD < THG))) {
					// 3 order - QUICKEST method
					PSD_f[ix] = PSD_f[ix] + (pow(CourNum_f[ix],2) - pow(1.,2.))/(3*2*1)
														* (D2 + d3 / 2); // A.16

				} else {
					// 4 order - PSD_f is still second order
					PSD_f[ix] = PSD_f[ix] + (pow(CourNum_f[ix],2) - pow(1.,2.))/(3*2*1)
													* (D2 + fabs(CourNum_f[ix])/4 * d3);

					D4 = (PSD_U[2] - 3*PSD_U[1] + 2*PSD_C + 2*PSD_D[1] - 3*PSD_D[2] + PSD_D[3])/2;
					d5 = PSD_U[2] - 5*PSD_U[1] + 10*PSD_C - 10*PSD_D[1] + 5*PSD_D[2] - PSD_D[3];

					if (always_5rd_order) { // !never_5rd_order && always_5rd_order) {
						// 5 order
						PSD_f[ix] = PSD_f[ix] + (pow(CourNum_f[ix],2) - pow(1.,2.)) * (pow(CourNum_f[ix],2) - pow(2.,2.))/(5*4*3*2*1)
														* (D4 + d5 / 2);

					} else {
						// 6 order - PSD_f is still 4th order. 
						PSD_f[ix] = PSD_f[ix] + (pow(CourNum_f[ix],2) - pow(1.,2.)) * (pow(CourNum_f[ix],2) - pow(2.,2.))/(5*4*3*2*1)
														* (D4 + fabs(CourNum_f[ix])/6 * d5);

						D6 = (PSD_U[3] - 5*PSD_U[2] + 9*PSD_U[1] - 5*PSD_C - 5*PSD_D[1] + 9*PSD_D[2] - 5*PSD_D[3] + PSD_D[4])/2;
						d7 = PSD_U[3] - 7*PSD_U[2] + 21*PSD_U[1] - 35*PSD_C + 35*PSD_D[1] - 21*PSD_D[2] + 7*PSD_D[3] - PSD_D[4];

						if (always_7rd_order) { // !never_7rd_order && (always_7rd_order || (CURVAV < THC2 && GRAD < THG))) {
							// 7 order
							PSD_f[ix] = PSD_f[ix] + (pow(CourNum_f[ix],2) - pow(1.,2.)) * (pow(CourNum_f[ix],2) - pow(2.,2.)) * (pow(CourNum_f[ix],2) - pow(3.,2.))/(7*6*5*4*3*2*1)
															* (D6 + d7 / 2);
						} else {
							// I think there is a bug somewhere here for 8-9 orders :( <- this comment is not from hayley
							// 8 order
							PSD_f[ix] = PSD_f[ix] + (pow(CourNum_f[ix],2) - pow(1.,2.)) * (pow(CourNum_f[ix],2) - pow(2.,2.)) * (pow(CourNum_f[ix],2) - pow(3.,2.))/(7*6*5*4*3*2*1)
															* (D6 + fabs(CourNum_f[ix])/8 * d7);

							// 8 even:  +1 -7 +20 -28 +14  +14  -28 +20 -7 +1
							// 9 odd: +1 -9 +36 -84 +126 -126 +84 -36 +9 -1
							D8 = (PSD_U[4] - 7*PSD_U[3] + 20*PSD_U[2] - 28*PSD_U[1] + 14*PSD_C + 14*PSD_D[1] - 28*PSD_D[2] + 20*PSD_D[3] - 7*PSD_D[4] + PSD_D[5])/2;
							d9 = PSD_U[4] - 9*PSD_U[3] + 36*PSD_U[2] - 84*PSD_U[1] + 126*PSD_C - 126*PSD_D[1] + 84*PSD_D[2] - 36*PSD_D[3] + 9*PSD_D[4] - PSD_D[5];
							// 9 order
							PSD_f[ix] = PSD_f[ix] + (pow(CourNum_f[ix],2) - pow(1.,2.)) * (pow(CourNum_f[ix],2) - pow(2.,2.)) * (pow(CourNum_f[ix],2) - pow(3.,2.)) * (pow(CourNum_f[ix],2) - pow(4.,2.))/(9*8*7*6*5*4*3*2*1)
															* (D8 + d9 / 2);
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
			dif3 = PSD_C    - PSD_D[1]; //
			dif4 = PSD_U[1] - PSD_C;
			dif5 = PSD_U[2] - PSD_U[1];
			dif6 = PSD_U[3] - PSD_D[2];

			double treshhold = 0;//1e-99;
			if (use_limiting && (!use_discreminator || !( // if don't use the descriminator - always go inside
//					   (((dif0 > treshhold && dif1 > treshhold) && (dif2 < -treshhold && dif3 < -treshhold)) && ((fabs(dif1) < fabs(dif0) - treshhold) && (fabs(dif2) < fabs(dif3) - treshhold))) ||
					   (((dif1 > treshhold && dif2 > treshhold) && (dif3 < -treshhold && dif4 < -treshhold)) && ((fabs(dif2) < fabs(dif1) - treshhold) && (fabs(dif3) < fabs(dif4) - treshhold)))
					|| (((dif2 > treshhold && dif3 > treshhold) && (dif4 < -treshhold && dif5 < -treshhold)) && ((fabs(dif3) < fabs(dif2) - treshhold) && (fabs(dif4) < fabs(dif5) - treshhold)))
					|| (((dif3 > treshhold && dif4 > treshhold) && (dif5 < -treshhold && dif6 < -treshhold)) && ((fabs(dif4) < fabs(dif3) - treshhold) && (fabs(dif5) < fabs(dif6) - treshhold)))
			))) {


//			if (!(
//					(((dif1 > 0 && dif2 > 0) && (dif3 < 0 && dif4 < 0)) && ((fabs(dif2) < fabs(dif1)) && (fabs(dif3) < fabs(dif4)))) ||
//					(((dif2 > 0 && dif3 > 0) && (dif4 < 0 && dif5 < 0)) && ((fabs(dif3) < fabs(dif2)) && (fabs(dif4) < fabs(dif5)))) ||
//					(((dif3 > 0 && dif4 > 0) && (dif5 < 0 && dif6 < 0)) && ((fabs(dif4) < fabs(dif3)) && (fabs(dif5) < fabs(dif6))))
//				)) {

				// limiting
				if (ACURV >= ADEL) { // (4)
					// not monotionic, or something
					PSD_f[ix] = PSD_C;
				} else
				if (fabs(CourNum_f[ix]) > 1e-31	) { // this will not work if velocity is zero
					// cont. from page 48
					// !!! PSD_REF = PSD_U + (PSD_C - PSD_U)/CourNum[ix]; // (3)
					// !!! PSD_REF = PSD_U + (PSD_C - PSD_U)/CourNum_f; // (3)
					PSD_REF = PSD_U[1] + (PSD_C - PSD_U[1]) / fabs(CourNum_f[ix]); // (3)

					if (DEL > 0) { // (5) limit
						PSD_f[ix] = fmax(PSD_f[ix], PSD_C); // max
						PSD_f[ix] = fmin(PSD_f[ix], fmin(PSD_REF, PSD_D[1])); // min, min
					} else if (DEL < 0) { // (6)
						PSD_f[ix] = fmin(PSD_f[ix], PSD_C); // min
						PSD_f[ix] = fmax(PSD_f[ix], fmax(PSD_REF, PSD_D[1])); // max, max
					}
				}
			}


		}


		// update PSD
		double PSD_r, PSD_l, CourNum_r, CourNum_l; // values on left and right faces (between (i and i-1), and (i and i+1))
		for (ix = 0; ix < x_size-1; ix++) {

			if (ix == 0 && x_LBC_type != "BCT_PERIODIC") { // special case
				continue; // need to calculate ix == 0 point only for periodic;
			}

			// (7)
			PSD_r     = PSD_f[ix+1];
			CourNum_r = CourNum_f[ix+1];
			PSD_l     = PSD_f[ix];
			CourNum_l = CourNum_f[ix];

			//PSD[ix] = PSD[ix] - CourNum[ix] * (PSD_r - PSD_l); - old version that makes the total PSD oscillate
			PSD[ix] = PSD[ix] - (CourNum_r*PSD_r - CourNum_l*PSD_l);
		}
		if (x_LBC_type == "BCT_PERIODIC") { // special case
			PSD[x_size-1] = PSD[0]; // need to update ix == x_size-1 point for periodic;
		}
		if (x_UBC_type == "BCT_CONSTANT_VALUE") {
            PSD[x_size-1] = x_UBC;
		}

	}

	return true;
}
