/**
 * \file MatrixSolver.cpp
 *
 * Making model matrixes, solving model matrixes.
 *
 * Matrix form of linear equations: A*X[t+1] = B*X[t],
 * where A - model matrix, B - RHS, X[t] - known values of function (PSD), X[t+1] - unknown values of function.
 *
 * In that file there are procedures for making model matrix for 1d-diffusion, 2d-diffusion, some  ideas of 3d-diffusion and mixed terms.
 * Solver for tridiagonal matrix, solver by gauss method and iteration method (upper relaxation).
 *
 * This file is under development, and has a lot of commented code, unfinished etc code.
 *
 * There is a checked version for 1d diffusion (or split method of 2d, 3d diffusions) - 1d_universal_solver.
 *
 * \author Developed under supervision of the PI Yuri Shprits
 */

//#define CPPL_DEBUG
//#define CPPL_VERBOSE

// Memory leaks
#if defined(_WIN32) || defined(_WIN64)
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

#include "MatrixSolver.h"
#include <math.h>
#include <malloc.h>
#include <string>
#include <ctime>
#include <iostream>
#include <stdlib.h>

using namespace std;

/** Supportive sub-function to add boundary conditions to model matrix
 */
void AddBoundary(DiagMatrix &matr_A, string type, int in, int id1, double dh) {
	int id0 = 0;
	if (type == "BCT_CONSTANT_VALUE") { // for condition on value
		matr_A[id0][in] = 1;
		matr_A[id1][in] = 0;
	} else if (type == "BCT_CONSTANT_DERIVATIVE") { // for condition on derivative
		// !!! Works incorrectly for a derivative != 0
		matr_A[id0][in] = 1 / dh;
		matr_A[id1][in] = -1 / dh;
	} else {
		printf("2D_DIFF_BOUNDARY: unknown boundary type: %s", type.c_str());
		abort();
	}
}

bool AddBoundaries_1D(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	Matrix1D<double> &x,
	int x_size,
	double x_LBC, double x_UBC,
	string x_LBC_type, string x_UBC_type,
	int ix) {

	int in, id1;
	double dh;

	in = ix;

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
			abort();
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
			abort();
		}
	}
	return true;
}


bool AddBoundaries_2D(
		CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
		Matrix2D<double> &x, Matrix2D<double> &y,
		int x_size, int y_size,
		Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
		Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
		string &x_LBC_type, string &x_UBC_type,
		string &y_LBC_type, string &y_UBC_type,
		int ix, int iy, int in) {

	int id1;
	double dh;

	// Bboundary conditions
	if (ix == 0 && x_size >= 3) {

		matr_C[0][in] = x_LBC[iy];
		id1 = matr_A.index1d(ix + 1, iy) - in;
		dh = x[ix + 1][iy] - x[ix][iy];
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
			abort();
		}

	} else if (ix == x_size - 1 && x_size >= 3) {

		matr_C[0][in] = x_UBC[iy];
		id1 = matr_A.index1d(ix - 1, iy) - in;
		dh = x[ix][iy] - x[ix - 1][iy];
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
			abort();
		}

	} else if (iy == 0 && y_size >= 3) {

		matr_C[0][in] = y_LBC[ix];
		id1 = matr_A.index1d(ix, iy + 1) - in;
		dh = y[ix][iy + 1] - y[ix][iy];
		//AddBoundary(matr_A, y_LBC_type, in, id, dh);
		int id0 = 0;
		if (y_LBC_type == "BCT_CONSTANT_VALUE") { // for condition on value
			matr_A[id0][in] = 1;
			matr_A[id1][in] = 0;
		} else if (y_LBC_type == "BCT_CONSTANT_DERIVATIVE") { // for condition on derivative
			// !!! Works incorrectly for a derivative != 0
			matr_A[id0][in] = 1 / dh;
			matr_A[id1][in] = -1 / dh;
		} else {
			printf("2D_DIFF_BOUNDARY: unknown boundary type: %s", y_LBC_type.c_str());
			abort();
		}

	} else if (iy == y_size - 1 && y_size >= 3) {

		matr_C[0][in] = y_UBC[ix];
		id1 = matr_A.index1d(ix, iy - 1) - in;
		dh = y[ix][iy] - y[ix][iy - 1];
		//AddBoundary(matr_A, y_UBC_type, in, id, dh);
		int id0 = 0;
		if (y_UBC_type == "BCT_CONSTANT_VALUE") { // for condition on value
			matr_A[id0][in] = 1;
			matr_A[id1][in] = 0;
		} else if (y_UBC_type == "BCT_CONSTANT_DERIVATIVE") { // for condition on derivative
			// !!! Works incorrectly for a derivative != 0
			matr_A[id0][in] = 1 / dh;
			matr_A[id1][in] = -1 / dh;
		} else {
			printf("2D_DIFF_BOUNDARY: unknown boundary type: %s", y_UBC_type.c_str());
			abort();
		}

	} else {
		// error - not a boundary
	}

	return true;
}


// Implicit diagonals, explicit mixed
bool MakeModelMatrix_2D_ADI1_x(CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
		Matrix2D<double> &x, Matrix2D<double> &y,
		int x_size, int y_size,
		Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
		Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
		string x_LBC_type, string x_UBC_type,
		string y_LBC_type, string y_UBC_type,
		Matrix2D<double> &Dxx, Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
		Matrix2D<double> &G, double dt) {

	// Make diagonals to be equal to zero
	DiagMatrix::iterator it;
	for (it = matr_A.begin(); it != matr_A.end(); it++)	it->second = 0;
	for (it = matr_B.begin(); it != matr_B.end(); it++)	it->second = 0;
	for (it = matr_C.begin(); it != matr_C.end(); it++)	it->second = 0;

	// create a new model matrix
	// (f^{t+1} - f^{t})/dt = L1(f^{t+1}) + L2(f^{t+1}) + L3(f^{t+1})[main equation, losses are calculated separately]
	// L1, L2, L3 - diffusion operators
	//  + need to add there the equations for boundary conditions also
	double dh;
	int ix, iy, in, id;
	for (ix = 0; ix < x_size; ix++) {
		for (iy = 0; iy < y_size; iy++) {
			// calculating current line number (in)
			in = matr_A.index1d(ix, iy);

			// Boundary conditions
			if (ix == 0 && x_size >= 3) {

				matr_C[0][in] = x_LBC[iy];
				id = matr_A.index1d(ix + 1, iy) - in;
				dh = x[ix + 1][iy] - x[ix][iy];
				AddBoundary(matr_A, x_LBC_type, in, id, dh);

			} else if (ix == x_size - 1 && x_size >= 3) {

				matr_C[0][in] = x_UBC[iy];
				id = matr_A.index1d(ix - 1, iy) - in;
				dh = x[ix][iy] - x[ix - 1][iy];
				AddBoundary(matr_A, x_UBC_type, in, id, dh);

			} else if (iy == 0 && y_size >= 3) {

				matr_C[0][in] = y_LBC[ix];
				id = matr_A.index1d(ix, iy + 1) - in;
				dh = y[ix][iy + 1] - y[ix][iy];
				AddBoundary(matr_A, y_LBC_type, in, id, dh);

			} else if (iy == y_size - 1 && y_size >= 3) {

				matr_C[0][in] = y_UBC[ix];
				id = matr_A.index1d(ix, iy - 1) - in;
				dh = y[ix][iy] - y[ix][iy - 1];
				AddBoundary(matr_A, y_UBC_type, in, id, dh);

			} else {

				// now we are sure we are not on a boundary, can do the Fokker-Planck equation approximation in the inner area

				// f^{t+1}/dt
				matr_A[0][in] += 1.0 / dt;
				// f^{t}/dt
				matr_B[0][in] += 1.0 / dt;

				// Dxx
				if (x_size >= 3) {
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "x_left", "x_right", x, y, Dxx, G, -0.5);
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "x_right", "x_left", x, y, Dxx, G, -0.5);
				}

				// Dyy
				if (y_size >= 3) {
					//SecondDerivativeApproximation_2D(matr_B, ix, iy, "y_left", "y_right", x, y, Dyy, G, 0.5);
					//SecondDerivativeApproximation_2D(matr_B, ix, iy, "y_right", "y_left", x, y, Dyy, G, 0.5);
				}

				// mixed
				if (x_size >= 3 && y_size >= 3) {

					// d/dx Dxy d/dy f
					// Multiplicatior = 0.125 because we approximate each term 4 times (0.25) and do it for each direction (/2)
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "y_left", "x_right", x, y, Dxy, G, 0.125);
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "y_right", "x_left", x, y, Dxy, G, 0.125);
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "y_left", "x_left", x, y, Dxy, G, 0.125);
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "y_right", "x_right", x, y, Dxy, G, 0.125);

					// d/dy Dyx d/dx f
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "x_left", "y_right", x, y, Dyx, G, 0.125);
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "x_right", "y_left", x, y, Dyx, G, 0.125);
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "x_left", "y_left", x, y, Dyx, G, 0.125);
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "x_right", "y_right", x, y, Dyx, G, 0.125);

				}

			}
		}
	}
	// Output::echo("recalculated.\n");

	// save the time of matrix change
	matr_A.change_ind = clock();
	matr_B.change_ind = clock();
	matr_C.change_ind = clock();

	return true;
}

bool MakeModelMatrix_2D_ADI1_y(CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
		Matrix2D<double> &x, Matrix2D<double> &y,
		int x_size, int y_size,
		Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
		Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
		string x_LBC_type, string x_UBC_type,
		string y_LBC_type, string y_UBC_type,
		Matrix2D<double> &Dxx, Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
		Matrix2D<double> &G, double dt) {

	// Make diagonals to be equal to zero
	DiagMatrix::iterator it;
	for (it = matr_A.begin(); it != matr_A.end(); it++)	it->second = 0;
	for (it = matr_B.begin(); it != matr_B.end(); it++)	it->second = 0;
	for (it = matr_C.begin(); it != matr_C.end(); it++)	it->second = 0;

	// create a new model matrix
	// (f^{t+1} - f^{t})/dt = L1(f^{t+1}) + L2(f^{t+1}) + L3(f^{t+1})[main equation, losses are calculated separately]
	// L1, L2, L3 - diffusion operators
	//  + need to add there the equations for boundary conditions also
	double dh;
	int ix, iy, in, id;
	for (ix = 0; ix < x_size; ix++) {
		for (iy = 0; iy < y_size; iy++) {
			// calculating current line number (in)
			in = matr_A.index1d(iy, ix);

			// Bboundary conditions
			if (ix == 0 && x_size >= 3) {

				matr_C[0][in] = x_LBC[iy];
				id = matr_A.index1d(iy, ix + 1) - in;
				dh = x[ix + 1][iy] - x[ix][iy];
				AddBoundary(matr_A, x_LBC_type, in, id, dh);

			} else if (ix == x_size - 1 && x_size >= 3) {

				matr_C[0][in] = x_UBC[iy];
				id = matr_A.index1d(iy, ix - 1) - in;
				dh = x[ix][iy] - x[ix - 1][iy];
				AddBoundary(matr_A, x_UBC_type, in, id, dh);

			} else if (iy == 0 && y_size >= 3) {

				matr_C[0][in] = y_LBC[ix];
				id = matr_A.index1d(iy + 1, ix) - in;
				dh = y[ix][iy + 1] - y[ix][iy];
				AddBoundary(matr_A, y_LBC_type, in, id, dh);

			} else if (iy == y_size - 1 && y_size >= 3) {

				matr_C[0][in] = y_UBC[ix];
				id = matr_A.index1d(iy - 1, ix) - in;
				dh = y[ix][iy] - y[ix][iy - 1];
				AddBoundary(matr_A, y_UBC_type, in, id, dh);

			} else {

				// now we are sure we are not on a boundary, can do the Fokker-Planck equation approximation in the inner area

				// f^{t+1}/dt
				matr_A[0][in] += 1.0 / dt;
				// f^{t}/dt
				matr_B[0][in] += 1.0 / dt;

				// Dxx
				if (x_size >= 3) {
					//SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "x_left", "x_right", x, y, Dxx, G, 0.5);
					//SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "x_right", "x_left", x, y, Dxx, G, 0.5);
				}

				// Dyy
				if (y_size >= 3) {
					SecondDerivativeApproximation_2D_y(matr_A, ix, iy, "y_left", "y_right", x, y, Dyy, G, -0.5);
					SecondDerivativeApproximation_2D_y(matr_A, ix, iy, "y_right", "y_left", x, y, Dyy, G, -0.5);
				}

				// mixed
				if (x_size >= 3 && y_size >= 3) {


					// Multiplicatior = 0.125 because we approximate each term 4 times (0.25) and do it for each direction (/2)
					// d/dx Dxy d/dy f
					SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "y_left", "x_right", x, y, Dxy, G, 0.125);
					SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "y_right", "x_left", x, y, Dxy, G, 0.125);
					SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "y_left", "x_left", x, y, Dxy, G, 0.125);
					SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "y_right", "x_right", x, y, Dxy, G, 0.125);
					// d/dy Dyx d/dx f
					SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "x_left", "y_right", x, y, Dyx, G, 0.125);
					SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "x_right", "y_left", x, y, Dyx, G, 0.125);
					SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "x_left", "y_left", x, y, Dyx, G, 0.125);
					SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "x_right", "y_right", x, y, Dyx, G, 0.125);

				}

			}
		}
	}
	// Output::echo("recalculated.\n");

	// save the time of matrix change
	matr_A.change_ind = clock();
	matr_B.change_ind = clock();
	matr_C.change_ind = clock();

	return true;
}

bool MakeModelMatrix_2D_ADI2_x(CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
		Matrix2D<double> &x, Matrix2D<double> &y,
		int x_size, int y_size,
		Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
		Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
		string x_LBC_type, string x_UBC_type,
		string y_LBC_type, string y_UBC_type,
		Matrix2D<double> &Dxx, Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
		Matrix2D<double> &G, double dt) {

	// Make diagonals to be equal to zero
	DiagMatrix::iterator it;
	for (it = matr_A.begin(); it != matr_A.end(); it++)	it->second = 0;
	for (it = matr_B.begin(); it != matr_B.end(); it++)	it->second = 0;
	for (it = matr_C.begin(); it != matr_C.end(); it++)	it->second = 0;

	Matrix2D<double> Coef1(x_size, y_size), Coef2(x_size, y_size);
	Coef1 = 1; Coef2 = 1;

	// create a new model matrix
	// (f^{t+1} - f^{t})/dt = L1(f^{t+1}) + L2(f^{t+1}) + L3(f^{t+1})[main equation, losses are calculated separately]
	// L1, L2, L3 - diffusion operators
	//  + need to add there the equations for boundary conditions also
	double dh;
	int ix, iy, in, id;
	double coef;
	for (ix = 0; ix < x_size; ix++) {
		for (iy = 0; iy < y_size; iy++) {
			// calculating current line number (in)
			in = matr_A.index1d(ix, iy);

			// Bboundary conditions
			if (ix == 0 && x_size >= 3) {

				matr_C[0][in] = x_LBC[iy];
				id = matr_A.index1d(ix + 1, iy) - in;
				dh = x[ix + 1][iy] - x[ix][iy];
				AddBoundary(matr_A, x_LBC_type, in, id, dh);

			} else if (ix == x_size - 1 && x_size >= 3) {

				matr_C[0][in] = x_UBC[iy];
				id = matr_A.index1d(ix - 1, iy) - in;
				dh = x[ix][iy] - x[ix - 1][iy];
				AddBoundary(matr_A, x_UBC_type, in, id, dh);

			} else if (iy == 0 && y_size >= 3) {

				matr_C[0][in] = y_LBC[ix];
				id = matr_A.index1d(ix, iy + 1) - in;
				dh = y[ix][iy + 1] - y[ix][iy];
				AddBoundary(matr_A, y_LBC_type, in, id, dh);

			} else if (iy == y_size - 1 && y_size >= 3) {

				matr_C[0][in] = y_UBC[ix];
				id = matr_A.index1d(ix, iy - 1) - in;
				dh = y[ix][iy] - y[ix][iy - 1];
				AddBoundary(matr_A, y_UBC_type, in, id, dh);

			} else {

				// now we are sure we are not on a boundary, can do the Fokker-Planck equation approximation in the inner area

				// f^{t+1}/dt
				matr_A[0][in] += 1.0 / dt;
				// f^{t}/dt
				matr_B[0][in] += 1.0 / dt;

				// Dxx
				if (x_size >= 3) {
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "x_left", "x_right", x, y, Dxx, G, -0.5);
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "x_right", "x_left", x, y, Dxx, G, -0.5);
				}

				// Dyy
				if (y_size >= 3) {
					//SecondDerivativeApproximation_2D(matr_B, ix, iy, "y_left", "y_right", x, y, Dyy, G, 0.5);
					//SecondDerivativeApproximation_2D(matr_B, ix, iy, "y_right", "y_left", x, y, Dyy, G, 0.5);
				}

				// mixed
				if (x_size >= 3 && y_size >= 3) {


					// mixed terms half implicit-half explicit
					//explicit part
					// Dxy d2f/dxy
					AnySecondDerivativeApproximation_2D(matr_B, ix, iy, "y_left", "x_right", x, y, Dxy, Coef2, 0.125);
					AnySecondDerivativeApproximation_2D(matr_B, ix, iy, "y_right", "x_left", x, y, Dxy, Coef2, 0.125);
					AnySecondDerivativeApproximation_2D(matr_B, ix, iy, "y_left", "x_left", x, y, Dxy, Coef2, 0.125);
					AnySecondDerivativeApproximation_2D(matr_B, ix, iy, "y_right", "x_right", x, y, Dxy, Coef2, 0.125);
					// Dyx d2f/dyx
					AnySecondDerivativeApproximation_2D(matr_B, ix, iy, "x_left", "y_right", x, y, Dyx, Coef2, 0.125);
					AnySecondDerivativeApproximation_2D(matr_B, ix, iy, "x_right", "y_left", x, y, Dyx, Coef2, 0.125);
					AnySecondDerivativeApproximation_2D(matr_B, ix, iy, "x_left", "y_left", x, y, Dyx, Coef2, 0.125);
					AnySecondDerivativeApproximation_2D(matr_B, ix, iy, "x_right", "y_right", x, y, Dyx, Coef2, 0.125);


					// Implicit part (first derivatives)
					// -1/G * d(G Dxy)/dx * df/dy
					coef = -1 / G[ix][iy] * (
							(Dxy[ix+1][iy] * G[ix+1][iy] - Dxy[ix][iy] * G[ix][iy]) / (x[ix+1][iy] - x[ix][iy]) +
							(Dxy[ix][iy] * G[ix][iy] - Dxy[ix-1][iy] * G[ix-1][iy]) / (x[ix][iy] - x[ix-1][iy]))/2;

					coef = coef / 2;
					//if (coef < 0) {
						id = matr_A.index1d(ix, iy+1) - in;
						matr_A[id][in] += +coef / (y[ix][iy+1] - y[ix][iy]);

						id = matr_A.index1d(ix, iy) - in;
						matr_A[id][in] += -coef / (y[ix][iy+1] - y[ix][iy]);
					//} else {
						id = matr_A.index1d(ix, iy) - in;
						matr_A[id][in] += +coef / (y[ix][iy] - y[ix][iy-1]);

						id = matr_A.index1d(ix, iy-1) - in;
						matr_A[id][in] += -coef / (y[ix][iy] - y[ix][iy-1]);
					//}

					// -1/G * d(G Dyx)/dy * df/dx
					coef = -1 / G[ix][iy] * (
							(Dyx[ix][iy+1] * G[ix][iy+1] - Dyx[ix][iy] * G[ix][iy]) / (y[ix][iy+1] - y[ix][iy]) +
							(Dyx[ix][iy] * G[ix][iy] - Dyx[ix][iy-1] * G[ix][iy-1]) / (y[ix][iy] - y[ix][iy-1]))/2;

					coef = coef / 2;
					//if (coef < 0) {
						id = matr_A.index1d(ix+1, iy) - in;
						matr_A[id][in] += +coef / (x[ix+1][iy] - x[ix][iy]);

						id = matr_A.index1d(ix, iy) - in;
						matr_A[id][in] += -coef / (x[ix+1][iy] - x[ix][iy]);
					//} else {
						id = matr_A.index1d(ix, iy) - in;
						matr_A[id][in] += +coef / (x[ix][iy] - x[ix-1][iy]);

						id = matr_A.index1d(ix-1, iy) - in;
						matr_A[id][in] += -coef / (x[ix][iy] - x[ix-1][iy]);
					//}

				}
			}
		}
	}
	// Output::echo("recalculated.\n");

	// save the time of matrix change
	matr_A.change_ind = clock();
	matr_B.change_ind = clock();
	matr_C.change_ind = clock();

	return true;
}

bool MakeModelMatrix_2D_ADI2_y(CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
		Matrix2D<double> &x, Matrix2D<double> &y,
		int x_size, int y_size,
		Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
		Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
		string x_LBC_type, string x_UBC_type,
		string y_LBC_type, string y_UBC_type,
		Matrix2D<double> &Dxx, Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
		Matrix2D<double> &G, double dt) {

	// Make diagonals to be equal to zero
	DiagMatrix::iterator it;
	for (it = matr_A.begin(); it != matr_A.end(); it++)	it->second = 0;
	for (it = matr_B.begin(); it != matr_B.end(); it++)	it->second = 0;
	for (it = matr_C.begin(); it != matr_C.end(); it++)	it->second = 0;

	Matrix2D<double> Coef1(x_size, y_size), Coef2(x_size, y_size);
	Coef1 = 1; Coef2 = 1;

	// create a new model matrix
	// (f^{t+1} - f^{t})/dt = L1(f^{t+1}) + L2(f^{t+1}) + L3(f^{t+1})[main equation, losses are calculated separately]
	// L1, L2, L3 - diffusion operators
	//  + need to add there the equations for boundary conditions also
	double dh;
	int ix, iy, in, id;
	double coef;
	for (ix = 0; ix < x_size; ix++) {
		for (iy = 0; iy < y_size; iy++) {
			// calculating current line number (in)
			in = matr_A.index1d(iy, ix);

			// Bboundary conditions
			if (ix == 0 && x_size >= 3) {

				matr_C[0][in] = x_LBC[iy];
				id = matr_A.index1d(iy, ix + 1) - in;
				dh = x[ix + 1][iy] - x[ix][iy];
				AddBoundary(matr_A, x_LBC_type, in, id, dh);

			} else if (ix == x_size - 1 && x_size >= 3) {

				matr_C[0][in] = x_UBC[iy];
				id = matr_A.index1d(iy, ix - 1) - in;
				dh = x[ix][iy] - x[ix - 1][iy];
				AddBoundary(matr_A, x_UBC_type, in, id, dh);

			} else if (iy == 0 && y_size >= 3) {

				matr_C[0][in] = y_LBC[ix];
				id = matr_A.index1d(iy + 1, ix) - in;
				dh = y[ix][iy + 1] - y[ix][iy];
				AddBoundary(matr_A, y_LBC_type, in, id, dh);

			} else if (iy == y_size - 1 && y_size >= 3) {

				matr_C[0][in] = y_UBC[ix];
				id = matr_A.index1d(iy - 1, ix) - in;
				dh = y[ix][iy] - y[ix][iy - 1];
				AddBoundary(matr_A, y_UBC_type, in, id, dh);

			} else {

				// now we are sure we are not on a boundary, can do the Fokker-Planck equation approximation in the inner area

				// f^{t+1}/dt
				matr_A[0][in] += 1.0 / dt;
				// f^{t}/dt
				matr_B[0][in] += 1.0 / dt;

				// Dxx
				if (x_size >= 3) {
					//SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "x_left", "x_right", x, y, Dxx, G, 0.5);
					//SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "x_right", "x_left", x, y, Dxx, G, 0.5);
				}

				// Dyy
				if (y_size >= 3) {
					SecondDerivativeApproximation_2D_y(matr_A, ix, iy, "y_left", "y_right", x, y, Dyy, G, -0.5);
					SecondDerivativeApproximation_2D_y(matr_A, ix, iy, "y_right", "y_left", x, y, Dyy, G, -0.5);
				}

				// mixed
				if (x_size >= 3 && y_size >= 3) {

					// mixed terms half implicit-half explicit
					//explicit part
					// Dxy d2f/dxy
					AnySecondDerivativeApproximation_2D_y(matr_B, ix, iy, "y_left", "x_right", x, y, Dxy, Coef2, 0.125);
					AnySecondDerivativeApproximation_2D_y(matr_B, ix, iy, "y_right", "x_left", x, y, Dxy, Coef2, 0.125);
					AnySecondDerivativeApproximation_2D_y(matr_B, ix, iy, "y_left", "x_left", x, y, Dxy, Coef2, 0.125);
					AnySecondDerivativeApproximation_2D_y(matr_B, ix, iy, "y_right", "x_right", x, y, Dxy, Coef2, 0.125);
					// Dyx d2f/dyx
					AnySecondDerivativeApproximation_2D_y(matr_B, ix, iy, "x_left", "y_right", x, y, Dyx, Coef2, 0.125);
					AnySecondDerivativeApproximation_2D_y(matr_B, ix, iy, "x_right", "y_left", x, y, Dyx, Coef2, 0.125);
					AnySecondDerivativeApproximation_2D_y(matr_B, ix, iy, "x_left", "y_left", x, y, Dyx, Coef2, 0.125);
					AnySecondDerivativeApproximation_2D_y(matr_B, ix, iy, "x_right", "y_right", x, y, Dyx, Coef2, 0.125);


					// Implicit part (first derivatives)
					// -1/G * d(G Dxy)/dx * df/dy
					coef = -1 / G[ix][iy] * (
							(Dxy[ix+1][iy] * G[ix+1][iy] - Dxy[ix][iy] * G[ix][iy]) / (x[ix+1][iy] - x[ix][iy]) +
							(Dxy[ix][iy] * G[ix][iy] - Dxy[ix-1][iy] * G[ix-1][iy]) / (x[ix][iy] - x[ix-1][iy]))/2;

					coef = coef / 2;
					//if (coef < 0) {
						id = matr_A.index1d(iy+1, ix) - in;
						matr_A[id][in] += +coef / (y[ix][iy+1] - y[ix][iy]);

						id = matr_A.index1d(iy, ix) - in;
						matr_A[id][in] += -coef / (y[ix][iy+1] - y[ix][iy]);
					//} else {
						id = matr_A.index1d(iy, ix) - in;
						matr_A[id][in] += +coef / (y[ix][iy] - y[ix][iy-1]);

						id = matr_A.index1d(iy-1, ix) - in;
						matr_A[id][in] += -coef / (y[ix][iy] - y[ix][iy-1]);
					//}

					// -1/G * d(G Dyx)/dy * df/dx
					coef = -1 / G[ix][iy] * (
							(Dyx[ix][iy+1] * G[ix][iy+1] - Dyx[ix][iy] * G[ix][iy]) / (y[ix][iy+1] - y[ix][iy]) +
							(Dyx[ix][iy] * G[ix][iy] - Dyx[ix][iy-1] * G[ix][iy-1]) / (y[ix][iy] - y[ix][iy-1]))/2;

					coef = coef / 2;
					//if (coef < 0) {
						id = matr_A.index1d(iy, ix+1) - in;
						matr_A[id][in] += +coef / (x[ix+1][iy] - x[ix][iy]);

						id = matr_A.index1d(iy, ix) - in;
						matr_A[id][in] += -coef / (x[ix+1][iy] - x[ix][iy]);
					//} else {
						id = matr_A.index1d(iy, ix) - in;
						matr_A[id][in] += +coef / (x[ix][iy] - x[ix-1][iy]);

						id = matr_A.index1d(iy, ix-1) - in;
						matr_A[id][in] += -coef / (x[ix][iy] - x[ix-1][iy]);
					//}

				}
			}
		}
	}
	// Output::echo("recalculated.\n");

	// save the time of matrix change
	matr_A.change_ind = clock();
	matr_B.change_ind = clock();
	matr_C.change_ind = clock();

	return true;
}


bool MakeModelMatrix_2D_ADI3_x(CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
		Matrix2D<double> &x, Matrix2D<double> &y,
		int x_size, int y_size,
		Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
		Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
		string x_LBC_type, string x_UBC_type,
		string y_LBC_type, string y_UBC_type,
		Matrix2D<double> &Dxx, Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
		Matrix2D<double> &G, double dt) {

	// Make diagonals to be equal to zero
	DiagMatrix::iterator it;
	for (it = matr_A.begin(); it != matr_A.end(); it++)	it->second = 0;
	for (it = matr_B.begin(); it != matr_B.end(); it++)	it->second = 0;
	for (it = matr_C.begin(); it != matr_C.end(); it++)	it->second = 0;

	// create a new model matrix
	// (f^{t+1} - f^{t})/dt = L1(f^{t+1}) + L2(f^{t+1}) + L3(f^{t+1})[main equation, losses are calculated separately]
	// L1, L2, L3 - diffusion operators
	//  + need to add there the equations for boundary conditions also
	double dh;
	int ix, iy, in, id;
	for (ix = 0; ix < x_size; ix++) {
		for (iy = 0; iy < y_size; iy++) {
			// calculating current line number (in)
			in = matr_A.index1d(ix, iy);

			// Bboundary conditions
/*			if (ix == 0 && x_size >= 3) {

				matr_C[0][in] = x_LBC[iy];
				id = matr_A.index1d(ix + 1, iy) - in;
				dh = x[ix + 1][iy] - x[ix][iy];
				AddBoundary(matr_A, x_LBC_type, in, id, dh);

			} else if (ix == x_size - 1 && x_size >= 3) {

				matr_C[0][in] = x_UBC[iy];
				id = matr_A.index1d(ix - 1, iy) - in;
				dh = x[ix][iy] - x[ix - 1][iy];
				AddBoundary(matr_A, x_UBC_type, in, id, dh);

			} else if (iy == 0 && y_size >= 3) {

				matr_C[0][in] = y_LBC[ix];
				id = matr_A.index1d(ix, iy + 1) - in;
				dh = y[ix][iy + 1] - y[ix][iy];
				AddBoundary(matr_A, y_LBC_type, in, id, dh);

			} else if (iy == y_size - 1 && y_size >= 3) {

				matr_C[0][in] = y_UBC[ix];
				id = matr_A.index1d(ix, iy - 1) - in;
				dh = y[ix][iy] - y[ix][iy - 1];
				AddBoundary(matr_A, y_UBC_type, in, id, dh);
*/
			if ( (ix == 0 && x_size >= 3)
					|| (ix == x_size - 1 && x_size >= 3)
					|| (iy == 0 && y_size >= 3)
					|| (iy == y_size - 1 && y_size >= 3)) {

				// if at the boundary
				// add boundary conditions
				AddBoundaries_2D(
					matr_A, matr_B, matr_C,
					x, y,
					x_size, y_size,
					x_LBC, x_UBC,
					//psd.ySlice(0), psd.ySlice(y.size_y-1), //y_LBC, y_UBC,
					y_LBC, y_UBC,
					x_LBC_type, x_UBC_type,
					//"BCT_CONSTANT_VALUE", "BCT_CONSTANT_VALUE", //y_LBC_type, y_UBC_type,
					y_LBC_type, y_UBC_type,
					ix, iy, in);

			} else {

				// now we are sure we are not on a boundary, can do the Fokker-Planck equation approximation in the inner area

				// f^{t+1}/dt
				matr_A[0][in] += 1.0 / dt;
				// f^{t}/dt
				matr_B[0][in] += 1.0 / dt;

				// Dxx
				if (x_size >= 3) {
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "x_left", "x_right", x, y, Dxx, G, -0.25);
					SecondDerivativeApproximation_2D(matr_A, ix, iy, "x_right", "x_left", x, y, Dxx, G, -0.25);

					SecondDerivativeApproximation_2D(matr_B, ix, iy, "x_left", "x_right", x, y, Dxx, G, 0.25);
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "x_right", "x_left", x, y, Dxx, G, 0.25);
				}

				// Dyy
				if (y_size >= 3) {
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "y_left", "y_right", x, y, Dyy, G, 0.5);
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "y_right", "y_left", x, y, Dyy, G, 0.5);
				}

				// mixed
				if (x_size >= 3 && y_size >= 3) {

					// d/dx Dxy d/dy f
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "y_left", "x_right", x, y, Dxy, G, 0.25);
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "y_right", "x_left", x, y, Dxy, G, 0.25);
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "y_left", "x_left", x, y, Dxy, G, 0.25);
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "y_right", "x_right", x, y, Dxy, G, 0.25);

					// d/dy Dyx d/dx f
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "x_left", "y_right", x, y, Dyx, G, 0.25);
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "x_right", "y_left", x, y, Dyx, G, 0.25);
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "x_left", "y_left", x, y, Dyx, G, 0.25);
					SecondDerivativeApproximation_2D(matr_B, ix, iy, "x_right", "y_right", x, y, Dyx, G, 0.25);

				}

			}
		}
	}
	// Output::echo("recalculated.\n");

	// save the time of matrix change
	matr_A.change_ind = clock();
	matr_B.change_ind = clock();
	matr_C.change_ind = clock();

	return true;
}

bool MakeModelMatrix_2D_ADI3_y(CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
		Matrix2D<double> &x, Matrix2D<double> &y,
		int x_size, int y_size,
		Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
		Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
		string x_LBC_type, string x_UBC_type,
		string y_LBC_type, string y_UBC_type,
		Matrix2D<double> &Dxx, Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
		Matrix2D<double> &G, double dt) {

	// Make diagonals to be equal to zero
	DiagMatrix::iterator it;
	for (it = matr_A.begin(); it != matr_A.end(); it++)	it->second = 0;
	for (it = matr_B.begin(); it != matr_B.end(); it++)	it->second = 0;
	for (it = matr_C.begin(); it != matr_C.end(); it++)	it->second = 0;

	// create a new model matrix
	// (f^{t+1} - f^{t})/dt = L1(f^{t+1}) + L2(f^{t+1}) + L3(f^{t+1})[main equation, losses are calculated separately]
	// L1, L2, L3 - diffusion operators
	//  + need to add there the equations for boundary conditions also
	double dh;
	int ix, iy, in, id;
	for (ix = 0; ix < x_size; ix++) {
		for (iy = 0; iy < y_size; iy++) {
			// calculating current line number (in)
			in = matr_A.index1d(iy, ix);

			// Bboundary conditions
/*			if (ix == 0 && x_size >= 3) {

				matr_C[0][in] = x_LBC[iy];
				id = matr_A.index1d(iy, ix + 1) - in;
				dh = x[ix + 1][iy] - x[ix][iy];
				AddBoundary(matr_A, x_LBC_type, in, id, dh);

			} else if (ix == x_size - 1 && x_size >= 3) {

				matr_C[0][in] = x_UBC[iy];
				id = matr_A.index1d(iy, ix - 1) - in;
				dh = x[ix][iy] - x[ix - 1][iy];
				AddBoundary(matr_A, x_UBC_type, in, id, dh);

			} else if (iy == 0 && y_size >= 3) {

				matr_C[0][in] = y_LBC[ix];
				id = matr_A.index1d(iy + 1, ix) - in;
				dh = y[ix][iy + 1] - y[ix][iy];
				AddBoundary(matr_A, y_LBC_type, in, id, dh);

			} else if (iy == y_size - 1 && y_size >= 3) {

				matr_C[0][in] = y_UBC[ix];
				id = matr_A.index1d(iy - 1, ix) - in;
				dh = y[ix][iy] - y[ix][iy - 1];
				AddBoundary(matr_A, y_UBC_type, in, id, dh);
				*/

			if ( (ix == 0 && x_size >= 3)
					|| (ix == x_size - 1 && x_size >= 3)
					|| (iy == 0 && y_size >= 3)
					|| (iy == y_size - 1 && y_size >= 3)) {

				// if at the boundary
				// add boundary conditions
				AddBoundaries_2D(
					matr_A, matr_B, matr_C,
					x, y,
					x_size, y_size,
					//psd.xSlice(0), psd.xSlice(x.size_x-1),	//x_LBC, x_UBC,
					x_LBC, x_UBC,
					y_LBC, y_UBC,
					//"BCT_CONSTANT_VALUE", "BCT_CONSTANT_VALUE", //x_LBC_type, x_UBC_type,
					x_LBC_type, x_UBC_type,
					y_LBC_type, y_UBC_type,
					ix, iy, in);


			} else {

				// now we are sure we are not on a boundary, can do the Fokker-Planck equation approximation in the inner area

				// f^{t+1}/dt
				matr_A[0][in] += 1.0 / dt;
				// f^{t}/dt
				// !!!!!!!!!!!!!!! matr_B[0][in] += 1.0 / dt;

				// Dxx
				if (x_size >= 3) {
					//SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "x_left", "x_right", x, y, Dxx, G, 0.5);
					//SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "x_right", "x_left", x, y, Dxx, G, 0.5);
				}

				// Dyy
				if (y_size >= 3) {
					SecondDerivativeApproximation_2D_y(matr_A, ix, iy, "y_left", "y_right", x, y, Dyy, G, -0.25);
					SecondDerivativeApproximation_2D_y(matr_A, ix, iy, "y_right", "y_left", x, y, Dyy, G, -0.25);

					SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "y_left", "y_right", x, y, Dyy, G, -0.25);
					SecondDerivativeApproximation_2D_y(matr_B, ix, iy, "y_right", "y_left", x, y, Dyy, G, -0.25);
				}

				// mixed
				if (x_size >= 3 && y_size >= 3) {

				}

			}
		}
	}
	// Output::echo("recalculated.\n");

	// save the time of matrix change
	matr_A.change_ind = clock();
	matr_B.change_ind = clock();
	matr_C.change_ind = clock();

	return true;
}


/**
 * Lapack inversion.
 *
 * A * X = B - equation
 */
void Lapack(DiagMatrix &A, Matrix1D<double> &B, Matrix1D<double> &X) {

	// Save A and B to check the solution at the end
	Matrix1D<double> B_res;
	B_res = B;

	// iterator for diagonals of the diagonal matrix
	DiagMatrix::iterator it;

	long m_size = A[0].size_q1;
	it = A.begin();
	long kl = -it->first; // first diagonal
	it = A.end();
	it--;
	long ku = it->first; // last diagonal
	long int NRHS(1), LDAB(2 * kl + ku + 1), *IPIV(new long[m_size]), LDB(B.size_q1), INFO(1);

	double *Array = new double[(kl+ku+kl+1)*m_size];
	double **newmat = new double*[m_size];
	for(int i=0; i<m_size; i++){ newmat[i] = &Array[i*(kl+ku+kl+1)]; }
	for (int i = 0; i < (kl+ku+kl+1)*m_size; i++) Array[i] = 0;

	int in, j;
	for (in = 0; in < m_size; in++) {
		for (it = A.begin(); it != A.end(); it++) {
			// Check if the element at line (in) and diagonal (it->first) is inside the matrix.
			j = in + it->first; // column number
			//if (in + it->first >= 0 && in + it->first < m_size) {
			if (j >= 0 && j < m_size) {
				// converting matrix, stored as diagonals, into lapack matrix (also diagonal)
				newmat[j][ku+kl-j+in] = it->second[in];
			}
		}
	}

	/*ofstream output("Lap_mat.dat");
	for (i = 0; i < m_size; i++) {
		for (j = 0; j < (kl+ku+kl+1); j++) {
			output << newmat[i][j] << "\t";
		}
		output << endl;
	}
	output.close();*/

	dgbsv_(&m_size, &kl, &ku, &NRHS, Array, &LDAB, IPIV, &B[0], &LDB, &INFO);

	if (INFO != 0) {
		printf("Lapack inversion Error!!! INFO = %ld.\n", INFO);
		abort();
	}

	// check
	// "A*X - B" should be zero
	// X is stored in B after Lapack solution
	// and we stored B in B_res previousely
	// So, it's "B_res - A*B" should be zero
	for (in = 0; in < m_size; in++) {
		for (it = A.begin(); it != A.end(); it++) {
			if (in + it->first >= 0 && in + it->first < m_size) {
				//cout << B_res[in] << "-=" << it->second[in] << "*" << B[in] << endl;
				B_res[in] -= it->second[in] * B[in + it->first];
			}
		}
	}

	// calculate max error
	double max = 0;
	for (int i = 0; i < m_size; i++) {
		max = (max > fabs(B_res(i))) ? max : B_res(i);
	}

	if(max>1) {
		printf(" Max error: %e.\n", max);
		abort();
	}

	delete IPIV;
	delete Array;
	delete newmat;
}

void SecondDerivativeApproximation_1D(CalculationMatrix &matr_A, int ix,
		string FirstDerivative, string SecondDerivative,
		Matrix1D<double> &x, ///< Coordinate x
		Matrix1D<double> &D, ///< Diffusion coefficient
		Matrix1D<double> &G, ///< Jacobian
		double multiplicator) {

	int dx1 = 0, dx2 = 0;
	if (FirstDerivative == "x_left") dx1 = -1;
	else dx1 = +1;
	if (SecondDerivative == "x_left") dx2 = -1;
	else dx2 = +1;

	// matr_A line number
	int in = matr_A.index1d(ix);

	double dh1, dh21, dh22;
	// The following trick work only for orthogonal grid
	// first and second dh in derivatove calculation

	// Grid steps.
	//dh1  = (x[ix][iy]             + y[ix][iy])             - (x[ix + dx1][iy + dy1]             + y[ix + dx1][iy + dy1]);
	//dh21 = (x[ix][iy]             + y[ix][iy])             - (x[ix + dx2][iy + dy2]             + y[ix + dx2][iy + dy2]);
	//dh22 = (x[ix + dx1][iy + dy1] + y[ix + dx1][iy + dy1]) - (x[ix + dx1 + dx2][iy + dy1 + dy2] + y[ix + dx1 + dx2][iy + dy1 + dy2]);
	dh1  = (x[ix])       - (x[ix + dx1]);
	dh21 = (x[ix])       - (x[ix + dx2]);
	dh22 = (x[ix + dx1]) - (x[ix + dx1 + dx2]);


	int id;

	// The same for all four coefficients, what's stending befire derivatives pretty much
	double common_part = multiplicator / G[ix] / dh1;

	// getting model matrix diagonal number according to derivative
	id = matr_A.index1d(ix) - in;
	// add corresponding coefficient for corresponding matrix cell
	matr_A[id][in] += +common_part / dh21 * D[ix] * G[ix];

	id = matr_A.index1d(ix + dx2) - in;
	matr_A[id][in] += -common_part / dh21 * D[ix] * G[ix];

	id = matr_A.index1d(ix + dx1) - in;
	matr_A[id][in] += -common_part / dh22 * D[ix + dx1] * G[ix + dx1];

	id = matr_A.index1d(ix + dx1 + dx2) - in;
	matr_A[id][in] += +common_part / dh22 * D[ix + dx1] * G[ix + dx1];

}


/**
 * Get change in indexes according to the derivatives direction.
 *
 * Used in approximation.
 * If it gets derivativeType == DT_ALPHA_left, for example, that means we have alpha-left-derivative, which is ( f(alpha) - f(alpha-1) ) / ( delta alpha ). So it returns dAlpha = -1, as a derivative direction vector.
 */
void GetDerivativeVector_2D(string derivativeType, int &dx, int &dy) {
	if (derivativeType == "x_left") {
		dx = -1;
	} else if (derivativeType == "x_right") {
		dx = 1;
	} else if (derivativeType == "y_left") {
		dy = -1;
	} else if (derivativeType == "y_right") {
		dy = 1;
	} else {
		printf("DERIVATIVE_TYPE: Unknown derivative approximation: %s\n",
				derivativeType.c_str());
	}
}

/**
 * Second derivative approximation, returns coefficients to be putted into the model matrix.
 *  \f$L_{\alpha \beta}(y) = (D_{\alpha \beta} \cdot y_{\bar{x}_\alpha})_{x_{\beta}}\f$
 *
 * Samarskiy, page 261
 *
 * Returns coefficients to be put into model matrix for an approximation of a second derivative.
 */
void SecondDerivativeApproximation_2D(CalculationMatrix &matr_A,
		int ix, int iy,
		string FirstDerivative, string SecondDerivative,
		Matrix2D<double> &x, ///< Coordinate x
		Matrix2D<double> &y, ///< Coordinate y
		Matrix2D<double> &D, ///< Diffusion coefficient
		Matrix2D<double> &G, ///< Jacobian
		double multiplicator) {

	int dx1 = 0, dy1 = 0;
	int dx2 = 0, dy2 = 0;
	// calculating derivative directions according derivative types
	// each second derivative has two directions, like
	// (alpha, alpha) - is a second derivative in alpha direction
	// (pc, pc) - second in pc direction
	// (pc, alpha) - mixed term
	// here we are getting that directions for derivative
	GetDerivativeVector_2D(FirstDerivative, dx1, dy1);
	GetDerivativeVector_2D(SecondDerivative, dx2, dy2);

	// matr_A line number
	int in = matr_A.index1d(ix, iy);

	double dh1, dh21, dh22;
	// The following trick work only for orthogonal grid
	// first and second dh in derivatove calculation

	// Grid steps.
	dh1  = (x[ix][iy]             + y[ix][iy])             - (x[ix + dx1][iy + dy1]             + y[ix + dx1][iy + dy1]);
	dh21 = (x[ix][iy]             + y[ix][iy])             - (x[ix + dx2][iy + dy2]             + y[ix + dx2][iy + dy2]);
	dh22 = (x[ix + dx1][iy + dy1] + y[ix + dx1][iy + dy1]) - (x[ix + dx1 + dx2][iy + dy1 + dy2] + y[ix + dx1 + dx2][iy + dy1 + dy2]);
	// !!!! dh1  = (x[ix][iy]       + y[ix][iy])       - (x[ix + dx1][iy]       + y[ix][iy + dy1]);
	// !!!! dh21 = (x[ix][iy]       + y[ix][iy])       - (x[ix + dx2][iy]       + y[ix][iy + dy2]);
	// !!!! dh22 = (x[ix + dx1][iy] + y[ix][iy + dy1]) - (x[ix + dx1 + dx2][iy] + y[ix][iy + dy1 + dy2]);


	int id;

	// The same for all four coefficients, what's stending befire derivatives pretty much
	double common_part = multiplicator / G[ix][iy] / dh1;

	// getting model matrix diagonal number according to derivative
	id = matr_A.index1d(ix, iy) - in;
	// add corresponding coefficient for corresponding matrix cell
	matr_A[id][in] += +common_part / dh21 * D[ix][iy] * G[ix][iy];

	id = matr_A.index1d(ix + dx2, iy + dy2) - in;
	matr_A[id][in] += -common_part / dh21 * D[ix][iy] * G[ix][iy];

	id = matr_A.index1d(ix + dx1, iy + dy1) - in;
	matr_A[id][in] += -common_part / dh22 * D[ix + dx1][iy + dy1] * G[ix + dx1][iy + dy1];

	id = matr_A.index1d(ix + dx1 + dx2, iy + dy1 + dy2) - in;
	matr_A[id][in] += +common_part / dh22 * D[ix + dx1][iy + dy1] * G[ix + dx1][iy + dy1];

}

void SecondDerivativeApproximation_2D_y(CalculationMatrix &matr_A,
		int ix, int iy,
		string FirstDerivative, string SecondDerivative,
		Matrix2D<double> &x, ///< Coordinate x
		Matrix2D<double> &y, ///< Coordinate y
		Matrix2D<double> &D, ///< Diffusion coefficient
		Matrix2D<double> &G, ///< Jacobian
		double multiplicator) {

	int dx1 = 0, dy1 = 0;
	int dx2 = 0, dy2 = 0;
	// calculating derivative directions according derivative types
	// each second derivative has two directions, like
	// (alpha, alpha) - is a second derivative in alpha direction
	// (pc, pc) - second in pc direction
	// (pc, alpha) - mixed term
	// here we are getting that directions for derivative
	GetDerivativeVector_2D(FirstDerivative, dx1, dy1);
	GetDerivativeVector_2D(SecondDerivative, dx2, dy2);

	// matr_A line number
	int in = matr_A.index1d(iy, ix);

	double dh1, dh21, dh22;
	// The following trick work only for orthogonal grid
	// first and second dh in derivatove calculation

	// Grid steps.
	dh1  = (x[ix][iy]             + y[ix][iy])             - (x[ix + dx1][iy + dy1]             + y[ix + dx1][iy + dy1]);
	dh21 = (x[ix][iy]             + y[ix][iy])             - (x[ix + dx2][iy + dy2]             + y[ix + dx2][iy + dy2]);
	dh22 = (x[ix + dx1][iy + dy1] + y[ix + dx1][iy + dy1]) - (x[ix + dx1 + dx2][iy + dy1 + dy2] + y[ix + dx1 + dx2][iy + dy1 + dy2]);

	int id;

	// The same for all four coefficients, what's stending befire derivatives pretty much
	double common_part = multiplicator / G[ix][iy] / dh1;

	// getting model matrix diagonal number according to derivative
	id = matr_A.index1d(iy, ix) - in;
	// add corresponding coefficient for corresponding matrix cell
	matr_A[id][in] += +common_part / dh21 * D[ix][iy] * G[ix][iy];

	id = matr_A.index1d(iy + dy2, ix + dx2) - in;
	matr_A[id][in] += -common_part / dh21 * D[ix][iy] * G[ix][iy];

	id = matr_A.index1d(iy + dy1, ix + dx1) - in;
	matr_A[id][in] += -common_part / dh22 * D[ix + dx1][iy + dy1] * G[ix + dx1][iy + dy1];

	id = matr_A.index1d(iy + dy1 + dy2, ix + dx1 + dx2) - in;
	matr_A[id][in] += +common_part / dh22 * D[ix + dx1][iy + dy1] * G[ix + dx1][iy + dy1];

}


/**
 * Second derivative approximation, returns coefficients to be putted into the model matrix.
 *  Coef1 * d/dx * Coef2 * df/dy
 *
 * Samarskiy, page 261
 *
 * Returns coefficients to be put into model matrix for an approximation of a second derivative.
 */
void AnySecondDerivativeApproximation_2D(CalculationMatrix &matr_A,
		int ix, int iy,
		string FirstDerivative, string SecondDerivative,
		Matrix2D<double> &x, ///< Coordinate x
		Matrix2D<double> &y, ///< Coordinate y
		Matrix2D<double> &Coef1, ///< Coefficient outside of the term
		Matrix2D<double> &Coef2, ///< Coefficient inside the first derivative
		double multiplicator) {

	int dx1 = 0, dy1 = 0;
	int dx2 = 0, dy2 = 0;
	// calculating derivative directions according derivative types
	// each second derivative has two directions, like
	// (alpha, alpha) - is a second derivative in alpha direction
	// (pc, pc) - second in pc direction
	// (pc, alpha) - mixed term
	// here we are getting that directions for derivative
	GetDerivativeVector_2D(FirstDerivative, dx1, dy1);
	GetDerivativeVector_2D(SecondDerivative, dx2, dy2);

	// matr_A line number
	int in = matr_A.index1d(ix, iy);

	double dh1, dh21, dh22;
	// The following trick work only for orthogonal grid
	// first and second dh in derivatove calculation

	// Grid steps.
	dh1  = (x[ix][iy]             + y[ix][iy])             - (x[ix + dx1][iy + dy1]             + y[ix + dx1][iy + dy1]);
	dh21 = (x[ix][iy]             + y[ix][iy])             - (x[ix + dx2][iy + dy2]             + y[ix + dx2][iy + dy2]);
	dh22 = (x[ix + dx1][iy + dy1] + y[ix + dx1][iy + dy1]) - (x[ix + dx1 + dx2][iy + dy1 + dy2] + y[ix + dx1 + dx2][iy + dy1 + dy2]);

	int id;

	// The same for all four coefficients, what's stending befire derivatives pretty much
	double common_part = multiplicator * Coef1[ix][iy] / dh1;

	// getting model matrix diagonal number according to derivative
	id = matr_A.index1d(ix, iy) - in;
	// add corresponding coefficient for corresponding matrix cell
	matr_A[id][in] += +common_part / dh21 * Coef2[ix][iy];

	id = matr_A.index1d(ix + dx2, iy + dy2) - in;
	matr_A[id][in] += -common_part / dh21 * Coef2[ix][iy];

	id = matr_A.index1d(ix + dx1, iy + dy1) - in;
	matr_A[id][in] += -common_part / dh22 * Coef2[ix + dx1][iy + dy1];

	id = matr_A.index1d(ix + dx1 + dx2, iy + dy1 + dy2) - in;
	matr_A[id][in] += +common_part / dh22 * Coef2[ix + dx1][iy + dy1];

}

void AnySecondDerivativeApproximation_2D_y(CalculationMatrix &matr_A,
		int ix, int iy,
		string FirstDerivative, string SecondDerivative,
		Matrix2D<double> &x, ///< Coordinate x
		Matrix2D<double> &y, ///< Coordinate y
		Matrix2D<double> &Coef1, ///< Coefficient outside of the term
		Matrix2D<double> &Coef2, ///< Coefficient inside the first derivative
		double multiplicator
		) {

	int dx1 = 0, dy1 = 0;
	int dx2 = 0, dy2 = 0;
	// calculating derivative directions according derivative types
	// each second derivative has two directions, like
	// (alpha, alpha) - is a second derivative in alpha direction
	// (pc, pc) - second in pc direction
	// (pc, alpha) - mixed term
	// here we are getting that directions for derivative
	GetDerivativeVector_2D(FirstDerivative, dx1, dy1);
	GetDerivativeVector_2D(SecondDerivative, dx2, dy2);

	// matr_A line number
	int in = matr_A.index1d(iy, ix);

	double dh1, dh21, dh22;
	// The following trick work only for orthogonal grid
	// first and second dh in derivatove calculation

	// Grid steps.
	dh1  = (x[ix][iy]             + y[ix][iy])             - (x[ix + dx1][iy + dy1]             + y[ix + dx1][iy + dy1]);
	dh21 = (x[ix][iy]             + y[ix][iy])             - (x[ix + dx2][iy + dy2]             + y[ix + dx2][iy + dy2]);
	dh22 = (x[ix + dx1][iy + dy1] + y[ix + dx1][iy + dy1]) - (x[ix + dx1 + dx2][iy + dy1 + dy2] + y[ix + dx1 + dx2][iy + dy1 + dy2]);

	int id;

	// The same for all four coefficients, what's standing before derivatives pretty much
	double common_part = multiplicator * Coef1[ix][iy] / dh1;

	// getting model matrix diagonal number according to derivative
	id = matr_A.index1d(iy, ix) - in;
	// add corresponding coefficient for corresponding matrix cell
	matr_A[id][in] += +common_part / dh21 * Coef2[ix][iy];

	id = matr_A.index1d(iy + dy2, ix + dx2) - in;
	matr_A[id][in] += -common_part / dh21 * Coef2[ix][iy];

	id = matr_A.index1d(iy + dy1, ix + dx1) - in;
	matr_A[id][in] += -common_part / dh22 * Coef2[ix + dx1][iy + dy1];

	id = matr_A.index1d(iy + dy1 + dy2, ix + dx1 + dx2) - in;
	matr_A[id][in] += +common_part / dh22 * Coef2[ix + dx1][iy + dy1];
}


/// Solve the AU=R system of equations, where A - tridiagonal matrix nxn with diagonals a[], b[], c[].
bool tridag(double a[], double b[], double c[], double r[], double u[], long n) {
	long j;
	double bet, *gam;
	gam = new double[n];

	if (b[0] == 0.0) {
		printf("TRIDAG_MATRIX_ERROR: tridag: error, b[0] = 0");
		//log1::fcout << "error 1 in tridag\n";
		//return -1;
	}

	u[0]=r[0]/(bet=b[0]);

	for (j=1;j<=n-1;j++) {
		gam[j]=c[j-1]/bet;
		bet=b[j]-a[j]*gam[j];

		if (bet == 0.0) {
			printf("TRIDAG_MATRIX_ERROR: tridag: error, bet = 0");
			//log1::fcout << "tridag(: Error 2 in tridag\n";
			//return -1;
		}

		u[j]=(r[j]-a[j]*u[j-1])/bet;
	}

	for (j=(n-2);j>=0;j--)
		u[j] -= gam[j+1]*u[j+1];

	delete gam;

	return true;
}
