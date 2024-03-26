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
 * \brief Solves model matrices for diffusion calculations and has functionality for derivative approximations of matrices
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
#include <cmath>
#include <stdlib.h>
#include <string>
#include <ctime>
#include <cmath>
#include <iostream>
#include <stdlib.h>

using namespace std;

/** Supportive sub-function to add boundary conditions to model matrix
 */
void AddBoundary(DiagMatrix &matr_A, BoundaryConditionType type, int in, int id1, double dh) 
{
	int id0 = 0;
	if (type == BoundaryConditionType::ConstantValue) { // for condition on value
		matr_A[id0][in] = 1;
		matr_A[id1][in] = 0;
	} else if (type == BoundaryConditionType::ConstantDerivative) { // for condition on derivative
		// !!! Works incorrectly for a derivative != 0
		matr_A[id0][in] = 1 / dh;
		matr_A[id1][in] = -1 / dh;
	} else {
		std::cout << "2D_DIFF_BOUNDARY: unknown boundary type: " << type << std::endl;
		exit(EXIT_FAILURE);
	}
}

/**
* Adds in the upper and lower boundaries for the 1D case when x is at first or last index
* \param matr_A,B,C - calculation matrices
* \param x_LBC, x_UBC, x_LBC_type, x_UBC_type  - predefined 1D boundaries/types
*/
bool AddBoundaries_1D(
	CalculationMatrix &matr_A, CalculationMatrix &matr_C,
	const Matrix1D<double> &x,
	int x_size,
	double x_LBC, double x_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	int ix) 
{
	int in, id1;
	double dh;

	in = ix;

	if (ix == 0 && x_size >= 3) {

		matr_C[0][in] = x_LBC;
		id1 = matr_A.index1d(ix + 1) - in;
		dh = x[ix + 1] - x[ix];
		//AddBoundary(matr_A, x_LBC_type, in, id, dh);
		int id0 = 0;
		if (x_LBC_type == BoundaryConditionType::ConstantValue) { // for condition on value
			matr_A[id0][in] = 1;
			matr_A[id1][in] = 0;
		} else if (x_LBC_type == BoundaryConditionType::ConstantDerivative) { // for condition on derivative
			// !!! Works incorrectly for a derivative != 0
			matr_A[id0][in] = 1 / dh;
			matr_A[id1][in] = -1 / dh;
		} else {
			std::cout << "2D_DIFF_BOUNDARY: unknown boundary type: " << x_LBC_type << std::endl;
			exit(EXIT_FAILURE);
		}

	} else if (ix == x_size - 1 && x_size >= 3) {

		matr_C[0][in] = x_UBC;
		id1 = matr_A.index1d(ix - 1) - in;
		dh = x[ix] - x[ix - 1];
		//AddBoundary(matr_A, x_UBC_type, in, id, dh);
		int id0 = 0;
		if (x_UBC_type == BoundaryConditionType::ConstantValue) { // for condition on value
			matr_A[id0][in] = 1;
			matr_A[id1][in] = 0;
		} else if (x_UBC_type == BoundaryConditionType::ConstantDerivative) { // for condition on derivative
			// !!! Works incorrectly for a derivative != 0
			matr_A[id0][in] = 1 / dh;
			matr_A[id1][in] = -1 / dh;
		} else {
			std::cout << "2D_DIFF_BOUNDARY: unknown boundary type: " << x_UBC_type << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	return true;
}


/**
* Adds in the upper and lower boundaries for the 2D case when x or y is at first or last index
* \param matr_A,B,C - calculation matrices
* \param LBC, UBC for x and y - predefined 2D boundaries/types
*/
bool AddBoundaries_2D(
		CalculationMatrix &matr_A, Matrix1D<double> &vec_C,
		const Matrix2D<double> &x, const Matrix2D<double> &y,
		int x_size, int y_size,
		const Matrix1D<double>& x_LBC, const Matrix1D<double>& x_UBC,
		const Matrix1D<double>& y_LBC, const Matrix1D<double>& y_UBC,
		BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
		BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
		int ix, int iy, int in) 
{
	int id1;
	double dh;

	// Boundary conditions
	if (ix == 0 && x_size >= 3) {

		vec_C[in] = x_LBC[iy];
		id1 = matr_A.index1d(ix + 1, iy) - in;
		dh = x[ix + 1][iy] - x[ix][iy];
		//AddBoundary(matr_A, x_LBC_type, in, id, dh);
		int id0 = 0;
		if (x_LBC_type == BoundaryConditionType::ConstantValue) { // for condition on value
			matr_A[id0][in] = 1;
			matr_A[id1][in] = 0;
		} else if (x_LBC_type == BoundaryConditionType::ConstantDerivative) { // for condition on derivative
			// !!! Works incorrectly for a derivative != 0
			matr_A[id0][in] = 1 / dh;
			matr_A[id1][in] = -1 / dh;
		} else {
			std::cout << "2D_DIFF_BOUNDARY: unknown boundary type: " << x_LBC_type << std::endl;
			exit(EXIT_FAILURE);
		}

	} else if (ix == x_size - 1 && x_size >= 3) {

		vec_C[in] = x_UBC[iy];
		id1 = matr_A.index1d(ix - 1, iy) - in;
		dh = x[ix][iy] - x[ix - 1][iy];
		//AddBoundary(matr_A, x_UBC_type, in, id, dh);
		int id0 = 0;
		if (x_UBC_type == BoundaryConditionType::ConstantValue) { // for condition on value
			matr_A[id0][in] = 1;
			matr_A[id1][in] = 0;
		} else if (x_UBC_type == BoundaryConditionType::ConstantDerivative) { // for condition on derivative
			// !!! Works incorrectly for a derivative != 0
			matr_A[id0][in] = 1 / dh;
			matr_A[id1][in] = -1 / dh;
		} else {
			std::cout << "2D_DIFF_BOUNDARY: unknown boundary type: " << x_UBC_type << std::endl;
			exit(EXIT_FAILURE);
		}

	} else if (iy == 0 && y_size >= 3) {

		vec_C[in] = y_LBC[ix];
		id1 = matr_A.index1d(ix, iy + 1) - in;
		dh = y[ix][iy + 1] - y[ix][iy];
		//AddBoundary(matr_A, y_LBC_type, in, id, dh);
		int id0 = 0;
		if (y_LBC_type == BoundaryConditionType::ConstantValue) { // for condition on value
			matr_A[id0][in] = 1;
			matr_A[id1][in] = 0;
		} else if (y_LBC_type == BoundaryConditionType::ConstantDerivative) { // for condition on derivative
			// !!! Works incorrectly for a derivative != 0
			matr_A[id0][in] = 1 / dh;
			matr_A[id1][in] = -1 / dh;
		} else {
			std::cout << "2D_DIFF_BOUNDARY: unknown boundary type: " << y_LBC_type << std::endl;
			exit(EXIT_FAILURE);
		}

	} else if (iy == y_size - 1 && y_size >= 3) {

		vec_C[in] = y_UBC[ix];
		id1 = matr_A.index1d(ix, iy - 1) - in;
		dh = y[ix][iy] - y[ix][iy - 1];
		//AddBoundary(matr_A, y_UBC_type, in, id, dh);
		int id0 = 0;
		if (y_UBC_type == BoundaryConditionType::ConstantValue) { // for condition on value
			matr_A[id0][in] = 1;
			matr_A[id1][in] = 0;
		} else if (y_UBC_type == BoundaryConditionType::ConstantDerivative) { // for condition on derivative
			// !!! Works incorrectly for a derivative != 0
			matr_A[id0][in] = 1 / dh;
			matr_A[id1][in] = -1 / dh;
		} else {
			std::cout << "2D_DIFF_BOUNDARY: unknown boundary type: " << y_UBC_type << std::endl;;
			exit(EXIT_FAILURE);
		}

	} else {
		// error - not a boundary
	}

	return true;
}

/**
 * Lapack solver for general banded matrix equations A * x = rhs
 * \param A Matrix A as DiagMatrix
 * \param rhs right-hand side vector, overwritten with solution
 */
void Lapack(const DiagMatrix &A, Matrix1D<double> &rhs) {

	long mat_rows = (long)A.at(0).size_q1;
	long rhs_rows = (long)rhs.size_q1;
	if(mat_rows != rhs_rows) {
		std::cerr << "Number of rows in matrix A must be equal to elements of rhs, but was " << mat_rows << " and " << rhs_rows << "!\n";
		exit(EXIT_FAILURE);
	}

	long kl = (long)(-A.cbegin()->first); // (absolute value of) lowest off-diagonal index
	long ku = (long)(A.crbegin()->first); // highest off-diagonal index
	long mat_cols = 2 * kl + ku + 1;
	long rhs_cols = 1;
	
	long* ipiv = new long[mat_rows]; // lapack array to store permutation indices
	long info = 1; // lapack error code

	double *flat_matrix = new double[mat_cols * mat_rows];
	for (int i = 0; i < mat_cols * mat_rows; i++) {
		flat_matrix[i] = 0;
	}

	// rearrange CalculationMatrix A to Lapack's general banded matrix format 
	for (const auto& [idx, diagonal] : A) {
		int diag_start = std::max(0, -idx);
		int diag_end = std::min(mat_rows, mat_rows - idx);
		int offset = ku + kl + idx * (mat_cols - 1);
		for (int in = diag_start; in < diag_end; in++) {
			flat_matrix[in * mat_cols + offset] = diagonal[in];
		}
	}

	// std::ofstream output("Lap_mat.dat");
	// for (int i = 0; i < mat_rows; i++) {
	// 	for (int j = 0; j < mat_cols; j++) {
	// 		output << flat_matrix[i * mat_cols + j] << '\t';
	// 	}
	// 	output << std::endl;
	// }
	// output.close();
	
	// Lapack's double-precision general banded (dgb) matrix solver
	dgbsv_(&mat_rows, &kl, &ku, &rhs_cols, flat_matrix, &mat_cols, ipiv, &rhs[0], &rhs_rows, &info);

	if (info != 0) {
		printf("Lapack inversion Error! INFO = %ld.\n", info);
		exit(EXIT_FAILURE);
	}

	delete[] ipiv;
	delete[] flat_matrix;
}

/**
* numerical derivative approximation of a matrix in 1D
*
*  \f$ Coef1 \frac{d}{dx} \f$
* which is approximatied using the method:
*
* \f$ \frac{m G(x) D(x)}{[f(x) - f'(x)]G(x)[(f(x) - f''(x))]} \f$
*
* Done for a couple locations near x[ix] each location has a slightly modified equation
*/
void SecondDerivativeApproximation_1D(CalculationMatrix &matr_A,
		int ix, ///< index for current location
		const string& FirstDerivative, ///< determines whether first derivative approx. is with the point to the right or left of ix
		const string& SecondDerivative, ///< determines whether second derivative approx. is with the point to the right or left of ix
		const Matrix1D<double> &x, ///< x matrix
		const Matrix1D<double> &D, ///< Diffusion coefficient
		const Matrix1D<double> &G, ///< Jacobian
		double multiplicator) ///< multiplier based on number of iterations done (m)
	{

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
void GetDerivativeVector_2D(const string& derivativeType, int &dx, int &dy) {
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
 * Second derivative approximation, returns coefficients to be put into the model matrix.
 *
 *  \f$ Coef1 \frac{d}{dx} * Coef2 \frac{df}{dy}  \f$
 *
 * The approximation is done in a similar fashion to SecondDerivativeApproximation_1D
 *
 * \f$ \frac{m G(x,y) D(x,y)}{[f(x,y) - f'(x , y)]G(x,y)[(f(x,y) - f''(x , y)]} \f$
 *
 * Returns coefficients to be put into model matrix for an approximation of a second derivative.
 */
void SecondDerivativeApproximation_2D(CalculationMatrix &matr_A,
		int ix, int iy,
		const string& FirstDerivative, const string& SecondDerivative,
		const Matrix2D<double> &x, ///< Coordinate x
		const Matrix2D<double> &y, ///< Coordinate y
		const Matrix2D<double> &D, ///< Diffusion coefficient
		const Matrix2D<double> &G, ///< Jacobian
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

	// The same for all four coefficients, what's standing before derivatives pretty much
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


/**
 * Second derivative approximation in the y direction first, returns coefficients to be put into the model matrix.
 *
*  \f$ Coef1 \frac{d}{dx} * Coef2 \frac{df}{dy}  \f$
 *
 * The approximation is done in a similar fashion to SecondDerivativeApproximation_1D
 *
 * \f$ \frac{m G(x,y) D(x,y)}{[f(x,y) - f'(x , y )]G(x,y)[(f(x,y) - f''(x , y)]} \f$
 *
 * Returns coefficients to be put into model matrix for an approximation of a second derivative.
 */
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
	// first and second dh in derivative calculation

	// Grid steps.
	dh1  = (x[ix][iy]             + y[ix][iy])             - (x[ix + dx1][iy + dy1]             + y[ix + dx1][iy + dy1]);
	dh21 = (x[ix][iy]             + y[ix][iy])             - (x[ix + dx2][iy + dy2]             + y[ix + dx2][iy + dy2]);
	dh22 = (x[ix + dx1][iy + dy1] + y[ix + dx1][iy + dy1]) - (x[ix + dx1 + dx2][iy + dy1 + dy2] + y[ix + dx1 + dx2][iy + dy1 + dy2]);

	int id;

	// The same for all four coefficients, what's standing before derivatives pretty much
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
 * Second derivative approximation, returns coefficients to be put into the model matrix.
 *  \f$ \frac {m Coef1 Coef2}{\frac{d}{dx} \frac{d}{dy}} \f$
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

	// The same for all four coefficients, what's standing before derivatives pretty much
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

/**
 * Second derivative approximation in the y direction first, returns coefficients to be put into the model matrix.
 *  \f$ \frac {m Coef1 Coef2}{\frac{d}{dx} \frac{d}{dy}} \f$
 */
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


/**
* Solver for a system of equations with 3-diagonal matrix.
*
* Au = r Where A = diag(a, b, c),
*
* \param a[]	- array, diagonal '-1' of the matrix
* \param b[]	- array, diagonal '0' of the matrix
* \param c[]	- array, diagonal '+1' of the matrix
* \param r[]	- array, r-vector - the Right Hand Side (RHS)
* \param u[]	- array, result
* \param n	- size of the matrix
*/
bool tridag(double a[], double b[], double c[], double r[], double u[], size_t n) {
	double bet, *gam;
	gam = new double[n];

	if (b[0] == 0.0) {
		printf("TRIDAG_MATRIX_ERROR: tridag: error, b[0] = 0");
		//log1::fcout << "error 1 in tridag\n";
		//return -1;
	}

	u[0]=r[0]/(bet=b[0]);

	for (size_t j = 1; j <= n - 1; j++) {
		gam[j] = c[j - 1] / bet;
		bet = b[j] - a[j] * gam[j];

		if (bet == 0.0) {
			printf("TRIDAG_MATRIX_ERROR: tridag: error, bet = 0");
			//log1::fcout << "tridag(: Error 2 in tridag\n";
			//return -1;
		}

		u[j] = (r[j] - a[j] * u[j - 1]) / bet;
	}

	for (size_t j = (n - 1); j >= 1; j--) {
		u[j - 1] -= gam[j] * u[j];
	}
	delete[] gam;

	return true;
}
