/**
 * \file MatrixSolver.h
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

#ifndef MatrixSolver_H
#define MatrixSolver_H

#include "Matrix.h"
#include "BoundaryConditionType.hpp"

void AddBoundary(DiagMatrix &Matrix_A, BoundaryConditionType type, int in, int id1, double dh);

bool AddBoundaries_1D(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	const Matrix1D<double> &x,
	int x_size,
	double x_LBC, double x_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	int ix);

bool AddBoundaries_2D(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	const Matrix2D<double> &x, const Matrix2D<double> &y,
	int x_size, int y_size,
	const Matrix1D<double>& x_LBC, const Matrix1D<double>& x_UBC,
	const Matrix1D<double>& y_LBC, const Matrix1D<double>& y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	int ix, int iy, int in
);

/// FUNCITON NOT IMPLEMENTED
bool MakeModelMatrix_2D(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	Matrix2D<double> &x, Matrix2D<double> &y,
	int x_size, int y_size,
	Matrix1D<double> &x_LBC, Matrix1D<double> &x_UBC,
	Matrix1D<double> &y_LBC, Matrix1D<double> &y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> &Dxx, Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
	Matrix2D<double> &G, Matrix2D<double> Sources, Matrix2D<double> Losses, double dt
);

void Lapack(DiagMatrix &A, Matrix1D<double> &B, Matrix1D<double> &X);

void SecondDerivativeApproximation_1D(
	CalculationMatrix &matr_A, int ix,
	const std::string& FirstDerivative,
	const std::string& SecondDerivative,
	const Matrix1D<double> &x, ///< Coordinate x
	const Matrix1D<double> &D, ///< Diffusion coefficient
	const Matrix1D<double> &G, ///< Jacobian
	double multiplicator
);

void GetDerivativeVector_2D(const std::string& derivativeType, int &dx, int &dy);

void SecondDerivativeApproximation_2D(
	CalculationMatrix &matr_A,
	int ix, int iy,
	const std::string& FirstDerivative, const std::string& SecondDerivative,
	const Matrix2D<double> &x, 				///< Coordinate x
	const Matrix2D<double> &y,  				///< Coordinate y
	const Matrix2D<double> &D, 				///< Diffusion coefficient
	const Matrix2D<double> &G, 				///< Jacobian
	double multiplicator
);

void SecondDerivativeApproximation_2D_y(
	CalculationMatrix &matr_A,
	int ix, int iy,
	std::string FirstDerivative, std::string SecondDerivative,
	Matrix2D<double> &x, 				///< Coordinate x
	Matrix2D<double> &y,  				///< Coordinate y
	Matrix2D<double> &D, 				///< Diffusion coefficient
	Matrix2D<double> &G, 				///< Jacobian
	double multiplicator
);

void AnySecondDerivativeApproximation_2D(
	CalculationMatrix &matr_A,
	int ix, int iy,
	std::string FirstDerivative, std::string SecondDerivative,
	Matrix2D<double> &x, ///< Coordinate x
	Matrix2D<double> &y, ///< Coordinate y
	Matrix2D<double> &Coef1, ///< Coefficient outside of the term
	Matrix2D<double> &Coef2, ///< Coefficient inside the first derivative
	double multiplicator
);

void AnySecondDerivativeApproximation_2D_y(
	CalculationMatrix &matr_A,
	int ix, int iy,
	std::string FirstDerivative, std::string SecondDerivative,
	Matrix2D<double> &x, ///< Coordinate x
	Matrix2D<double> &y, ///< Coordinate y
	Matrix2D<double> &Coef1, ///< Coefficient outside of the term
	Matrix2D<double> &Coef2, ///< Coefficient inside the first derivative
	double multiplicator
);


/// Solve the AU=R system of equations, where A - tridiagonal matrix nxn with diagonals a[], b[], c[].
bool tridag(double a[], double b[], double c[], double r[], double u[], long n);

/// Lapack function for matrix inversion declaration
extern "C" {
	//extern void dgesv_(int *,int *,double *,int *,int*,double *,int*,int*);
	//extern void sgesv_(int *,int *,float *,int *,int*, float *,int*,int*);
	extern void dgbsv_(long int *n, long int *kl, long int *ku, long int *nrhs, double *ab, long int *ldab, long int *ipiv, double *b, long int *ldb, long int *info);
};


bool MakeModelMatrix_2D_ADI1_x(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	Matrix2D<double> &x, Matrix2D<double> &y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> &Dxx, Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
	Matrix2D<double> &G, double dt
);

bool MakeModelMatrix_2D_ADI1_y(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	Matrix2D<double> &x, Matrix2D<double> &y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> &Dxx, Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
	Matrix2D<double> &G, double dt
);

bool MakeModelMatrix_2D_ADI2_x(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	Matrix2D<double> &x, Matrix2D<double> &y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> &Dxx, Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
	Matrix2D<double> &G, double dt
);

bool MakeModelMatrix_2D_ADI2_y(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	Matrix2D<double> &x, Matrix2D<double> &y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> &Dxx, Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
	Matrix2D<double> &G, double dt
);

bool MakeModelMatrix_2D_ADI3_y(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	Matrix2D<double> &x, Matrix2D<double> &y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> &Dxx, Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
	Matrix2D<double> &G, double dt
);

bool MakeModelMatrix_2D_ADI3_x(
	CalculationMatrix &matr_A, CalculationMatrix &matr_B, CalculationMatrix &matr_C,
	Matrix2D<double> &x, Matrix2D<double> &y,
	int x_size, int y_size,
	Matrix1D<double> x_LBC, Matrix1D<double> x_UBC,
	Matrix1D<double> y_LBC, Matrix1D<double> y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	Matrix2D<double> &Dxx, Matrix2D<double> &Dyy, Matrix2D<double> &Dxy, Matrix2D<double> &Dyx,
	Matrix2D<double> &G, double dt
);

void mkl_sparse_solve(double* values, int* columns, int* rowB, int* rowE, const double* rhs, double* solution, int m_size);
#endif

