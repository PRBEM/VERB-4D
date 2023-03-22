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

/**
 * 2D Diffusion - returns true on successful completion
 *
 * Make 2D arrays from 3D, create model matrices and solve them
 *
 * Method:
 *
 * 1. Add boundary conditions AddBoundaries_2D()
 *
 * 2. Add sources and losses
 *
 * 3. get the seond derivative approximation with diffusion coeficient SecondDerivativeApproximation_2D()
 *
 * 4. solve matrix with Lapack()
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
 * @param Dxx - 2D Diffusion matrix - derivative with respect to x twice
 * @param Dyy - 2D Diffusion matrix - derivative with respect to y twice
 * @param Dxy - 2D Diffusion matrix - derivative with respect to x then y
 * @param Dyx - 2D Diffusion matrix - derivative with respect to y then x
 * @param G - 2D used for Jacobian to normalize matrix
 * @param Sources - matrix used for Sources
 * @param Losses - Matrix used for Losses (loss cone)
 * @param dt - change in time of single time step
 */
bool Diffusion_2D(
		Matrix2D<double> &psd,
		const Matrix2D<double>& x, const Matrix2D<double>& y,
		int x_size, int y_size,
		const Matrix1D<double>& x_LBC, const Matrix1D<double>& x_UBC,
		const Matrix1D<double>& y_LBC, const Matrix1D<double>& y_UBC,
		BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
		BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
		const Matrix2D<double>& Dxx, const Matrix2D<double>& Dyy, const Matrix2D<double>& Dxy, const Matrix2D<double>& Dyx,
		const Matrix2D<double>& G, const Matrix2D<double>& Sources, const Matrix2D<double>& Losses, double dt) 
{
	
	// diagonals are zero-initialized on construction
	CalculationMatrix matr_A(x_size, y_size, 1, 1);

 	// only main diagonal non-zero, can be replaced with Matrix1D(x_size * y_size)
	CalculationMatrix matr_B(x_size, y_size, 1, 0);
	CalculationMatrix matr_C(x_size, y_size, 1, 0);

	// Create matrix form of the Fokker-Planck equation: matr_A * PSD(t+1) = matr_B * PSD(t) + matr_C
	//
	// Calculation matr_A, matr_B, and matr_C, i.e. numerical approximation of the derivatives
	//
	// (f^{t+1} - f^{t})/dt = L1(f^{t+1}) + L2(f^{t+1}) + L3(f^{t+1})[main equation, losses are calculated separately]
	// L1, L2, L3 - diffusion operators
	//  + need to add there the equations for boundary conditions also
	for (int ix = 0; ix < x_size; ix++) {
		for (int iy = 0; iy < y_size; iy++) {
			// calculating current line number (in)
			int in = matr_A.index1d(ix, iy);

			if((ix == 0 && x_size >= 3)	|| (ix == x_size - 1 && x_size >= 3)
			|| (iy == 0 && y_size >= 3) || (iy == y_size - 1 && y_size >= 3)) 
			{
				// if at the boundary add boundary conditions
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
	// matr_A.change_ind = clock();
	// matr_B.change_ind = clock();
	// matr_C.change_ind = clock();


	//matr_A.writeToFile("./Debug_output/matr_A.dat");
	//matr_B.writeToFile("./Debug_output/matr_B.dat");
	//matr_C.writeToFile("./Debug_output/matr_C.dat");


	// rearranged PSD into one 1d vector in column major form
	Matrix1D<double> rhs(matr_A.total_size);
	for (int ix = 0; ix < x_size; ix++) {
		for (int iy = 0; iy < y_size; iy++) {
			int in = matr_B.index1d(ix, iy);
			rhs[in] = psd[ix][iy];
		}
	}

	// make rhs = B*f + C
	rhs.times_equal(matr_B[0]);
	rhs += matr_C[0];

	//rhs.writeToFile("rhs.dat");
	Lapack(matr_A, rhs);
	//rhs.writeToFile("result.dat");


	// copy back from column major form
	for (int ix = 0; ix < x_size; ix++) {
		for (int iy = 0; iy < y_size; iy++) {
			int in = matr_B.index1d(ix, iy);
			psd[ix][iy] = rhs[in];
		}
	}

	return true;
}

