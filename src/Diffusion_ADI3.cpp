/**
 * \file Diffusion_ADI3.cpp
 *
 * \brief Calculating the diffusion in 2D using method from Jihye Shin and Sungsoo Kim (2008)
 *
 *  Created on: May 28, 2011
 *      Author: dimath
 */

#include "Diffusion_ADI3.h"
#include <iostream>

#include <ctime>

/** Create model matrices and solve the system,
 * The method is from
 * Jihye Shin and Sungsoo Kim (2008), ALTERNATING DIRECTION IMPLICIT METHOD FOR	TWO-DIMENSIONAL FOKKER-PLANCK EQUATION OF DENSE SPHERICAL STELLAR SYSTEMS
 *
 * http://arxiv.org/pdf/0805.0054v1.pdf
 *
 * \todo Fix this function - currently results in negative PSD values - using Lapack instead
 *
 * Method:
 *
 * For the x direction and then for the y direction
 *
 * 1. Add boundary conditions AddBoundary()
 *
 * 2. Add sources and Losses (only for the y direction)
 *
 * 3. get the seond derivative approximation with diffusion coeficient SecondDerivativeApproximation_2D()
 *
 * 4. Multiply matrix B by f (PSD_1D_x) and add matrix C - this will be the RHS for tridag()
 *
 * 5. solve matrix with tridag()
 *
 * @param psd - phase space density
 * @param x - one dimensional slice
 * @param y - one dimensional slice 
 * @param x_size - size of x slice
 * @param y_size - size of y slice
 * @param x_LBC - lower boundary condition for param x
 * @param x_UBC - upper boundary condition for param x
 * @param y_LBC - lower boundary condition for param y
 * @param y_UBC - upper boundary condition for param y
 * @param x_LBC_type - type of lower boundary condition for param x
 * @param x_UBC_type - type of upper boundary condition for param x
 * @param y_LBC_type - type of lower boundary condition for param y
 * @param y_UBC_type - type of upper boundary condition for param y
 * @param Dxx - 2D Diffusion matrix
 * @param Dyy - 2D Diffusion matrix
 * @param Dxy - 2D Diffusion matrix
 * @param Dyx - 2D Diffusion matrix
 * @param G - 2D used for Jacobian to normalize matrix
 * @param Sources - matrix used for Sources
 * @param Losses - Matrix used for Losses (loss cone)
 * @param dt - change in time of single time step
 */
bool Diffusion_2D_ADI3(
		Matrix2D<double> &psd,
		Matrix2D<double> x, Matrix2D<double> y,
		int x_size, int y_size,
		Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
		Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
		BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
		BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
		Matrix2D<double> Dxx, Matrix2D<double> Dyy, Matrix2D<double> Dxy, Matrix2D<double> Dyx,
		Matrix2D<double> G, Matrix2D<double> Sources, Matrix2D<double> Losses, double dt) 
{

	DiagMatrix::iterator it_B, it_A;
	Matrix1D<double> RHS(x_size * y_size);
	int ix, iy, in, id;
	DiagMatrix::iterator it;
	Matrix2D<double> psd_prev(x_size, y_size);
	Matrix1D<double> psd_1d_x(x_size * y_size); // Rearranged PSD into one vector of unknown variables
	Matrix1D<double> psd_1d_y(x_size * y_size); // Rearranged PSD into one vector of unknown variables
	Matrix1D<double> psd_1d_prev(x_size * y_size); // Rearranged PSD into one vector of unknown variables
	double dh;

	CalculationMatrix
		matr_A_x(x_size, y_size, 1, 1),
		matr_B_x(x_size, y_size, 1, 1),
		matr_C_x(x_size, y_size, 1, 0);


	for (ix = 0; ix < x_size; ix++) {
		for (iy = 0; iy < y_size; iy++) {
			in = matr_B_x.index1d(ix, iy);
			psd_1d_x[in] = psd[ix][iy];
		}
	}

/*	MakeModelMatrix_2D_ADI3_x(
					  matr_A_x, matr_B_x, matr_C_x,
					  x, y,
					  x_size, y_size,
					  x_LBC, x_UBC,
					  psd.ySlice(0), psd.ySlice(y.size_y-1), //y_LBC, y_UBC,
					  x_LBC_type, x_UBC_type,
					  "BCT_CONSTANT_VALUE", "BCT_CONSTANT_VALUE", //y_LBC_type, y_UBC_type,
					  Dxx, Dyy, Dxy, Dyx,
					  G, dt);*/

	// Make diagonals to be equal to zero
	for (it = matr_A_x.begin(); it != matr_A_x.end(); it++)	it->second = 0;
	for (it = matr_B_x.begin(); it != matr_B_x.end(); it++)	it->second = 0;
	for (it = matr_C_x.begin(); it != matr_C_x.end(); it++)	it->second = 0;

	// create a new model matrix
	// (f^{t+1} - f^{t})/dt = L1(f^{t+1}) + L2(f^{t+1}) + L3(f^{t+1})[main equation, losses are calculated separately]
	// L1, L2, L3 - diffusion operators
	//  + need to add there the equations for boundary conditions also
	for (ix = 0; ix < x_size; ix++) {
		for (iy = 0; iy < y_size; iy++) {
			// calculating current line number (in)
			in = matr_A_x.index1d(ix, iy);

			// Boundary conditions
			if (ix == 0 && x_size >= 3) {

				matr_C_x[0][in] = x_LBC[iy];
				id = matr_A_x.index1d(ix + 1, iy) - in;
				dh = x[ix + 1][iy] - x[ix][iy];
				AddBoundary(matr_A_x, x_LBC_type, in, id, dh);

			} else if (ix == x_size - 1 && x_size >= 3) {

				matr_C_x[0][in] = x_UBC[iy];
				id = matr_A_x.index1d(ix - 1, iy) - in;
				dh = x[ix][iy] - x[ix - 1][iy];
				AddBoundary(matr_A_x, x_UBC_type, in, id, dh);

			} else if (iy == 0 && y_size >= 3) {

				// Keep PSD as it was on that boundary, since we don't calculate diffusion in this direction
				//matr_C_x[0][in] = y_LBC[ix];
				matr_C_x[0][in] = psd[ix][0];
				id = matr_A_x.index1d(ix, iy + 1) - in;
				dh = y[ix][iy + 1] - y[ix][iy];
				//AddBoundary(matr_A_x, y_LBC_type, in, id, dh);
				AddBoundary(matr_A_x, BoundaryConditionType::ConstantValue, in, id, dh);

			} else if (iy == y_size - 1 && y_size >= 3) {

				// Keep PSD as it was on that boundary, since we don't calculate diffusion in this direction
				//matr_C_x[0][in] = y_UBC[ix];
				matr_C_x[0][in] = psd[ix][y.size_q2-1];
				id = matr_A_x.index1d(ix, iy - 1) - in;
				dh = y[ix][iy] - y[ix][iy - 1];
				//AddBoundary(matr_A_x, y_UBC_type, in, id, dh);
				AddBoundary(matr_A_x, BoundaryConditionType::ConstantValue, in, id, dh);

			} else {

				// now we are sure we are not on a boundary, can do the Fokker-Planck equation approximation in the inner area

				// f^{t+1}/dt
				matr_A_x[0][in] += 1.0 / dt;
				// f^{t}/dt
				matr_B_x[0][in] += 1.0 / dt;

				// Sources and losses
				// + Sources - Losses * f(t+1)
				//matr_C_x[0][in] += Sources[ix][iy];
				//matr_A_x[0][in] -= Losses[ix][iy];

				// Dxx
				if (x_size >= 3) {
					SecondDerivativeApproximation_2D(matr_A_x, ix, iy, "x_left", "x_right", x, y, Dxx, G, -0.25);
					SecondDerivativeApproximation_2D(matr_A_x, ix, iy, "x_right", "x_left", x, y, Dxx, G, -0.25);

					SecondDerivativeApproximation_2D(matr_B_x, ix, iy, "x_left", "x_right", x, y, Dxx, G, 0.25);
					SecondDerivativeApproximation_2D(matr_B_x, ix, iy, "x_right", "x_left", x, y, Dxx, G, 0.25);
				}

				// Dyy
				if (y_size >= 3) {
					SecondDerivativeApproximation_2D(matr_B_x, ix, iy, "y_left", "y_right", x, y, Dyy, G, 0.5);
					SecondDerivativeApproximation_2D(matr_B_x, ix, iy, "y_right", "y_left", x, y, Dyy, G, 0.5);
				}

				// mixed
				if (x_size >= 3 && y_size >= 3) {

					// d/dx Dxy d/dy f
					SecondDerivativeApproximation_2D(matr_B_x, ix, iy, "y_left", "x_right", x, y, Dxy, G, 0.25);
					SecondDerivativeApproximation_2D(matr_B_x, ix, iy, "y_right", "x_left", x, y, Dxy, G, 0.25);
					SecondDerivativeApproximation_2D(matr_B_x, ix, iy, "y_left", "x_left", x, y, Dxy, G, 0.25);
					SecondDerivativeApproximation_2D(matr_B_x, ix, iy, "y_right", "x_right", x, y, Dxy, G, 0.25);

					// d/dy Dyx d/dx f
					SecondDerivativeApproximation_2D(matr_B_x, ix, iy, "x_left", "y_right", x, y, Dyx, G, 0.25);
					SecondDerivativeApproximation_2D(matr_B_x, ix, iy, "x_right", "y_left", x, y, Dyx, G, 0.25);
					SecondDerivativeApproximation_2D(matr_B_x, ix, iy, "x_left", "y_left", x, y, Dyx, G, 0.25);
					SecondDerivativeApproximation_2D(matr_B_x, ix, iy, "x_right", "y_right", x, y, Dyx, G, 0.25);

				}

			}
		}
	}
	// Output::echo("recalculated.\n");

	// save the time of matrix change
	matr_A_x.change_ind = clock();
	matr_B_x.change_ind = clock();
	matr_C_x.change_ind = clock();


	//matr_A_x.writeToFile("matr_A_x.dat");
	//matr_B_x.writeToFile("matr_B_x.dat");
	//matr_C_x.writeToFile("matr_C_x.dat");

	// make RHS = B*f + C
	for (ix = 0; ix < x_size; ix++) {
		for (iy = 0; iy < y_size; iy++) {
			in = matr_B_x.index1d(ix, iy);
			RHS[in] = matr_C_x[0][in];
			for (it_B = matr_B_x.begin(); it_B != matr_B_x.end(); it_B++)
				// multiplication B * f
				if (in + it_B->first >= 0 && in + it_B->first < matr_B_x.total_size) {
					RHS[in] += it_B->second[in] * psd_1d_x[in + it_B->first];
				}
		}
	}

	tridag(
			&matr_A_x[-1][0], // just a way to get a pointer to first element of an array
			&matr_A_x[0][0],
			&matr_A_x[+1][0],
			&RHS[0],
			&psd_1d_x[0],
			psd_1d_x.size_q1);




	// Rearrange psd_1d for y-direction
	psd_prev = psd;

	// copy back
	for (ix = 0; ix < x_size; ix++) {
		for (iy = 0; iy < y_size; iy++) {
			in = matr_B_x.index1d(ix, iy);
			psd[ix][iy] = psd_1d_x[in];
		}
	}

	// y direction
	CalculationMatrix
		matr_A_y(y_size, x_size, 1, 1),
		matr_B_y(y_size, x_size, 1, 1),
		matr_C_y(y_size, x_size, 1, 0);

	for (ix = 0; ix < x_size; ix++) {
		for (iy = 0; iy < y_size; iy++) {
			in = matr_B_y.index1d(iy, ix);
			psd_1d_y[in] = psd[ix][iy];
			psd_1d_prev[in] = psd_prev[ix][iy];
		}
	}

/*
	MakeModelMatrix_2D_ADI3_y(
					  matr_A_y, matr_B_y, matr_C_y,
					  x, y,
					  x_size, y_size,
					  psd.xSlice(0), psd.xSlice(x.size_x-1),	//x_LBC, x_UBC,
					  y_LBC, y_UBC,
					  "BCT_CONSTANT_VALUE", "BCT_CONSTANT_VALUE", //x_LBC_type, x_UBC_type,
					  y_LBC_type, y_UBC_type,
					  Dxx, Dyy, Dxy, Dyx,
					  G, dt);*/

	// Make diagonals to be equal to zero
	for (it = matr_A_y.begin(); it != matr_A_y.end(); it++)	it->second = 0;
	for (it = matr_B_y.begin(); it != matr_B_y.end(); it++)	it->second = 0;
	for (it = matr_C_y.begin(); it != matr_C_y.end(); it++)	it->second = 0;

	// create a new model matrix
	// (f^{t+1} - f^{t})/dt = L1(f^{t+1}) + L2(f^{t+1}) + L3(f^{t+1})[main equation, losses are calculated separately]
	// L1, L2, L3 - diffusion operators
	//  + need to add there the equations for boundary conditions also
	for (ix = 0; ix < x_size; ix++) {
		for (iy = 0; iy < y_size; iy++) {
			// calculating current line number (in)
			in = matr_A_y.index1d(iy, ix);

			// Bboundary conditions
			if (ix == 0 && x_size >= 3) {

				// Keep PSD as it was on that boundary, since we don't calculate diffusion in this direction
				// matr_C_y[0][in] = x_LBC[iy];
				matr_C_y[0][in] = psd[0][iy];
				id = matr_A_y.index1d(iy, ix + 1) - in;
				dh = x[ix + 1][iy] - x[ix][iy];
				//AddBoundary(matr_A_y, x_LBC_type, in, id, dh);
				AddBoundary(matr_A_y, BoundaryConditionType::ConstantValue, in, id, dh);

			} else if (ix == x_size - 1 && x_size >= 3) {

				// Keep PSD as it was on that boundary, since we don't calculate diffusion in this direction
				// matr_C_y[0][in] = x_UBC[iy];
				matr_C_y[0][in] = psd[x.size_q1-1][iy];
				id = matr_A_y.index1d(iy, ix - 1) - in;
				dh = x[ix][iy] - x[ix - 1][iy];
				//AddBoundary(matr_A_y, x_UBC_type, in, id, dh);
				AddBoundary(matr_A_y, BoundaryConditionType::ConstantValue, in, id, dh);

			} else if (iy == 0 && y_size >= 3) {

				matr_C_y[0][in] = y_LBC[ix];
				id = matr_A_y.index1d(iy + 1, ix) - in;
				dh = y[ix][iy + 1] - y[ix][iy];
				AddBoundary(matr_A_y, y_LBC_type, in, id, dh);

			} else if (iy == y_size - 1 && y_size >= 3) {

				matr_C_y[0][in] = y_UBC[ix];
				id = matr_A_y.index1d(iy - 1, ix) - in;
				dh = y[ix][iy] - y[ix][iy - 1];
				AddBoundary(matr_A_y, y_UBC_type, in, id, dh);

			} else {

				// now we are sure we are not on a boundary, can do the Fokker-Planck equation approximation in the inner area

				// f^{t+1}/dt
				matr_A_y[0][in] += 1.0 / dt;
				// f^{t}/dt
				// !!!!!!!!!!!!!!! matr_B_y[0][in] += 1.0 / dt;

				// Sources and losses
				// + Sources - Losses * f(t+1)
				matr_C_y[0][in] += Sources[ix][iy];
				matr_A_y[0][in] -= Losses[ix][iy];

				// Dxx
				if (x_size >= 3) {
					//SecondDerivativeApproximation_2D_y(matr_B_y, ix, iy, "x_left", "x_right", x, y, Dxx, G, 0.5);
					//SecondDerivativeApproximation_2D_y(matr_B_y, ix, iy, "x_right", "x_left", x, y, Dxx, G, 0.5);
				}

				// Dyy
				if (y_size >= 3) {
					SecondDerivativeApproximation_2D_y(matr_A_y, ix, iy, "y_left", "y_right", x, y, Dyy, G, -0.25);
					SecondDerivativeApproximation_2D_y(matr_A_y, ix, iy, "y_right", "y_left", x, y, Dyy, G, -0.25);

					SecondDerivativeApproximation_2D_y(matr_B_y, ix, iy, "y_left", "y_right", x, y, Dyy, G, -0.25);
					SecondDerivativeApproximation_2D_y(matr_B_y, ix, iy, "y_right", "y_left", x, y, Dyy, G, -0.25);
				}

				// mixed
				if (x_size >= 3 && y_size >= 3) {

				}

			}
		}
	}
	// Output::echo("recalculated.\n");

	// save the time of matrix change
	matr_A_y.change_ind = clock();
	matr_B_y.change_ind = clock();
	matr_C_y.change_ind = clock();

	//matr_A_y.writeToFile("matr_A_y.dat");
	//matr_B_y.writeToFile("matr_B_y.dat");
	//matr_C_y.writeToFile("matr_C_y.dat");

	// make RHS = B*f + C
	for (iy = 0; iy < y_size; iy++) {
		for (ix = 0; ix < x_size; ix++) {
			in = matr_B_y.index1d(iy, ix);
			RHS[in] = matr_C_y[0][in];
			// This is from B-matrix, but B is multiplied to psd-prev, whil this one should be multiplied to psd
			if (iy != 0 && ix != 0 && ix != x_size-1 && iy != y_size -1)
				RHS[in] += 1.0 / dt * psd_1d_y[in];
			//cout << RHS[in] << endl;
			for (it_B = matr_B_y.begin(); it_B != matr_B_y.end(); it_B++)
				// multiplication B * f
				if (in + it_B->first >= 0 && in + it_B->first < matr_B_y.total_size) {
					//cout << "RHS["<< in <<"] = " << RHS[in] << " + " << it_B->second[in] << " * " << psd_1d[in + it_B->first] << " (" << it_B->first << ")";
					// !!!!!!!!! Multiplying to previous PSD
					RHS[in] += it_B->second[in] * psd_1d_prev[in + it_B->first];
					//cout << " = " << RHS[in] << endl;
				}
		}
	}

	//RHS.writeToFile("RHS.dat");

	tridag(
			&matr_A_y[-1][0], // just a way to get a pointer to first element of an array
			&matr_A_y[0][0],
			&matr_A_y[+1][0],
			&RHS[0],
			&psd_1d_y[0],
			psd_1d_y.size_q1);


	// copy back
	for (ix = 0; ix < x_size; ix++) {
		for (iy = 0; iy < y_size; iy++) {
			in = matr_B_y.index1d(iy, ix);
			//psd[ix][iy] = RHS[in];
			psd[ix][iy] = psd_1d_y[in];
		}
	}

	return true;
}
