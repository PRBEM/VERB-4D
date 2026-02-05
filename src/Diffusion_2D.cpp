// SPDX-FileCopyrightText: 2015 UCLA
// SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * \file Diffusion_2D.cpp
 *
 * \brief Getting the diffusion in 2D for the given a 2D matrix of Phase Space Densities, x, y, boundary conditions, jacobian, sources, losses and change in time
 *
 *  Created on: May 28, 2011
 *      Author: dimath
 */

#include "Diffusion_2D.h"
#include <iostream>


#ifdef LU_CACHING
bool Diffusion_2D(
		Matrix2D<double> &psd,
		const Matrix2D<double>& x, const Matrix2D<double>& y,
		int x_size, int y_size,
		const Matrix1D<double>& x_LBC, const Matrix1D<double>& x_UBC,
		const Matrix1D<double>& y_LBC, const Matrix1D<double>& y_UBC,
		BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
		BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
		const Matrix2D<double>& Dxx, const Matrix2D<double>& Dyy, const Matrix2D<double>& Dxy, const Matrix2D<double>& Dyx,
		const Matrix2D<double>& G, const Matrix2D<double>& Sources, const Matrix2D<double>& Losses, double dt, double sub_dt,
		double* lu_cache, long* index_cache, bool recompute_lu) 
#else
bool Diffusion_2D(
		Matrix2D<double> &psd,
		const Matrix2D<double>& x, const Matrix2D<double>& y,
		int x_size, int y_size,
		const Matrix1D<double>& x_LBC, const Matrix1D<double>& x_UBC,
		const Matrix1D<double>& y_LBC, const Matrix1D<double>& y_UBC,
		BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
		BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
		const Matrix2D<double>& Dxx, const Matrix2D<double>& Dyy, const Matrix2D<double>& Dxy, const Matrix2D<double>& Dyx,
		const Matrix2D<double>& G, const Matrix2D<double>& Sources, const Matrix2D<double>& Losses, double dt, double sub_dt) 
#endif
{

	const int num_substeps = std::max(1, int(dt / sub_dt));
	dt /= num_substeps;

	// diagonals are zero-initialized on construction
	CalculationMatrix
		matr_A(x_size, y_size, 1, 1),
		matr_B(x_size, y_size, 1, 0),
		matr_C(x_size, y_size, 1, 0);

	DiagMatrix::iterator it;
	for (it = matr_A.begin(); it != matr_A.end(); it++)	it->second = 0;
	for (it = matr_B.begin(); it != matr_B.end(); it++)	it->second = 0;
	for (it = matr_C.begin(); it != matr_C.end(); it++)	it->second = 0;

	//
	// Create matrix form of the Fokker-Planck equation: matr_A * PSD(t+1) = matr_B * PSD(t) + matr_C
	//
	// Calculation matr_A, matr_B, and matr_C, i.e. numerical approximation of the derivatives
	//
	// (f^{t+1} - f^{t})/dt = L1(f^{t+1}) + L2(f^{t+1}) + L3(f^{t+1})[main equation, losses are calculated separately]
	// L1, L2, L3 - diffusion operators
	//  + need to add there the equations for boundary conditions also
	int ix, iy, in;
	for (ix = 0; ix < x_size; ix++) {
		for (iy = 0; iy < y_size; iy++) {
			// calculating current line number (in)
			in = matr_A.index1d(ix, iy);

			if((ix == 0 && x_size >= 3)
				|| (ix == x_size - 1 && x_size >= 3)
				|| (iy == 0 		 && y_size >= 3)
				|| (iy == y_size - 1 && y_size >= 3)) {

				// if at the boundary
				// add boundary conditions
				AddBoundaries_2D(
					matr_A, matr_C,
					x, y, x_size, y_size,
					x_LBC, x_UBC, y_LBC, y_UBC,
					x_LBC_type, x_UBC_type, y_LBC_type, y_UBC_type,
					ix, iy, in);

			} else {

				// now we are sure we are not on a boundary, can do the Fokker-Planck equation approximation in the inner area

				// f^{t+1}/dt
				matr_A[0][in] += 1.0 / dt;
				// f^{t}/dt
				matr_B[0][in] += 1.0 / dt;

				// Sources and losses
				// + Sources - Losses * f(t+1)
				matr_C[0][in] += Sources[ix][iy];
				matr_A[0][in] -= Losses[ix][iy];

				// Dxx
				if (x_size >= 3) {
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "x_left", "x_right", x, y, Dxx, G, -0.5);
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "x_right", "x_left", x, y, Dxx, G, -0.5);
				}

				// Dyy
				if (y_size >= 3) {
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "y_left", "y_right", x, y, Dyy, G, -0.5);
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "y_right", "y_left", x, y, Dyy, G, -0.5);
				}

				// mixed
				if (x_size >= 3 && y_size >= 3) {
					// d/dx Dxy d/dy f
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "y_left", "x_right", x, y, Dxy, G, -0.25);
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "y_right", "x_left", x, y, Dxy, G, -0.25);
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "y_left", "x_left", x, y, Dxy, G, -0.25);
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "y_right", "x_right", x, y, Dxy, G, -0.25);

					// d/dy Dyx d/dx f
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "x_left", "y_right", x, y, Dyx, G, -0.25);
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "x_right", "y_left", x, y, Dyx, G, -0.25);
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "x_left", "y_left", x, y, Dyx, G, -0.25);
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "x_right", "y_right", x, y, Dyx, G, -0.25);
				}

			}
		}
	}
	// Output::echo("recalculated.\n");

	// save the time of matrix change
	matr_A.change_ind = clock();
	matr_B.change_ind = clock();
	matr_C.change_ind = clock();


	//matr_A.writeToFile("./Debug_output/matr_A.dat");
	//matr_B.writeToFile("./Debug_output/matr_B.dat");
	//matr_C.writeToFile("./Debug_output/matr_C.dat");

	//RHS.writeToFile("RHS.dat");
#ifdef LU_CACHING
	Lapack(matr_A, matr_B, matr_C, psd, num_substeps, lu_cache, index_cache, recompute_lu);
#else
	Lapack(matr_A, matr_B, matr_C, psd, num_substeps);
#endif
	//RHS.writeToFile("result.dat");

	return true;

}
