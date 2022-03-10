/**
 * \file Diffusion_1D.cpp
 *
 * \brief Getting the diffusion in 1D for the given a 1D matrix of Phase Space Densities, x, boundary conditions, jacobian, sources, losses and change in time
 *
 */

#include "Diffusion_1D.h"
#include <iostream>


/**
 * 1D Diffusion - returns true on successful completion
 *
 * Method:
 *
 * 1. Add boundary conditions
 *
 * 2. Add sources and losses
 *
 * 3. get the seond derivative approximation with diffusion coeficient SecondDerivativeApproximation_1D()
 *
 * 4. solve matrix with tridag()
 *
 * @param psd - phase space density
 * @param x - one dimensional slice 
 * @param x_size - size of slice
 * @param x_LBC - lower boundary condition for param x
 * @param x_UBC - upper boundary condition for param x
 * @param x_LBC_type - type of lower boundary condition for param x
 * @param x_UBC_type - type of upper boundary condition for param x
 * @param Dxx - 1D Diffusion matrix
 * @param G - 1D used for Jacobian to normalize matrix
 * @param Sources - matrix used for Sources
 * @param Losses - Matrix used for Losses (loss cone)
 * @param dt - change in time of single time step
 */
using std::string;
bool Diffusion_1D(
		Matrix1D<double> &psd,
		Matrix1D<double> x, int x_size,
		double x_LBC, double x_UBC,
		string x_LBC_type, string x_UBC_type,
		Matrix1D<double> Dxx, Matrix1D<double> G,
		Matrix1D<double> Sources, Matrix1D<double> Losses,
		double dt) {

	CalculationMatrix matr_A(x_size, 1, 1, 1), matr_B(x_size, 1, 1, 0), matr_C(x_size, 1, 1, 0);

	// Set all diagonals (1-D matrices) to 0
	DiagMatrix::iterator it;
	for (it = matr_A.begin(); it != matr_A.end(); it++)	it->second = 0;
	for (it = matr_B.begin(); it != matr_B.end(); it++)	it->second = 0;
	for (it = matr_C.begin(); it != matr_C.end(); it++)	it->second = 0;

	// create a new model matrix
	int ix, in, id1;
	double dh;
	for (ix = 0; ix < x_size; ix++) {
		// calculating current line number (in)
		in = ix;


		
		// Boundary conditions
		if (ix == 0 && x_size >= 3) {

			matr_C[0][in] = x_LBC;
			id1 = matr_A.index1d(ix + 1) - in;
			dh = x[ix + 1] - x[ix];
			//AddBoundary(matr_A, x_LBC_type, in, id, dh);
			int id0 = 0;
			if (x_LBC_type == "BCT_CONSTANT_VALUE") { // for condition on value
				matr_A[id0][in] = 1;
				matr_A[id1][in] = 0;
			} else if (x_LBC_type == "BCT_CONSTANT_DERIVATIVE") { // for condition on derivative
				// !!! Works incorrectly for a derivative != 0
				matr_A[id0][in] = 1 / dh;
				matr_A[id1][in] = -1 / dh;
			} else {
				printf("2D_DIFF_BOUNDARY: unknown boundary type: %s", x_LBC_type.c_str());
				exit(EXIT_FAILURE);
			}

		} else if (ix == x_size - 1 && x_size >= 3) {

			matr_C[0][in] = x_UBC;
			id1 = matr_A.index1d(ix - 1) - in;
			dh = x[ix] - x[ix - 1];
			//AddBoundary(matr_A, x_UBC_type, in, id, dh);
			int id0 = 0;
			if (x_UBC_type == "BCT_CONSTANT_VALUE") { // for condition on value
				matr_A[id0][in] = 1;
				matr_A[id1][in] = 0;
			} else if (x_UBC_type == "BCT_CONSTANT_DERIVATIVE") { // for condition on derivative
				// !!! Works incorrectly for a derivative != 0
				matr_A[id0][in] = 1 / dh;
				matr_A[id1][in] = -1 / dh;
			} else {
				printf("2D_DIFF_BOUNDARY: unknown boundary type: %s", x_UBC_type.c_str());
				exit(EXIT_FAILURE);
			}

		} else {

			// f^{t+1}/dt
			matr_A[0][in] += 1.0 / dt;
			// f^{t}/dt
			matr_B[0][in] += 1.0 / dt;

			// Sources and losses
			// + Sources - Losses * f(t+1)
			matr_C[0][in] += Sources[ix];
			matr_A[0][in] -= Losses[ix];

			// Dxx
			if (x_size >= 3) {
				// Writing the finite difference approximation of the diffusion operator into matrix A
				SecondDerivativeApproximation_1D(matr_A, ix, "x_left", "x_right", x, Dxx, G, -0.5);
				SecondDerivativeApproximation_1D(matr_A, ix, "x_right", "x_left", x, Dxx, G, -0.5);
			}
		}
	}

	// matr_A.writeToFile("matr_A.dat");
	// matr_B.writeToFile("matr_B.dat");
	// matr_C.writeToFile("matr_C.dat");

	// Add reanalysis here?

	// make RHS = B*f + C
	DiagMatrix::iterator it_B;
	Matrix1D<double> RHS(matr_A.total_size);
	for (ix = 0; ix < x_size; ix++) {
		in = matr_B.index1d(ix);
		RHS[in] = matr_C[0][in];
		for (it_B = matr_B.begin(); it_B != matr_B.end(); it_B++)
			// multiplication B * f
			if (in + it_B->first >= 0 && in + it_B->first < matr_B.total_size)
				RHS[in] += it_B->second[in] * psd[in + it_B->first];
	}

	//matr_A.writeToFile("./Debug_output/matr_A.dat");
	//matr_B.writeToFile("./Debug_output/matr_B.dat");
	//matr_C.writeToFile("./Debug_output/matr_C.dat");
	//RHS.writeToFile("./Debug_output/RHS.dat");
	//psd.writeToFile("./Debug_output/PSD.dat");

	tridag(
		&matr_A[-1][0], // just a way to get a pointer to first element of an array
		&matr_A[0][0],
		&matr_A[+1][0],
		&RHS[0],
		&psd[0],
		x_size
	);

	return true;
}
// ADDS BOUNDARIES FROM SECOND DERIVATIVE APPROXIMATION AND THEN SOLVES THE MATRIX USING TRIDAG
