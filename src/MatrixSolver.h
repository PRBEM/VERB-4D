/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

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
	CalculationMatrix &matr_A, CalculationMatrix &matr_C,
	const Matrix1D<double> &x,
	int x_size,
	double x_LBC, double x_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	int ix
);

bool AddBoundaries_2D(
	CalculationMatrix &matr_A, CalculationMatrix &matr_C,
	const Matrix2D<double> &x, const Matrix2D<double> &y,
	int x_size, int y_size,
	const Matrix1D<double>& x_LBC, const Matrix1D<double>& x_UBC,
	const Matrix1D<double>& y_LBC, const Matrix1D<double>& y_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
	int ix, int iy, int in
);

#ifdef LU_CACHING
void Lapack(
	const CalculationMatrix& A,
	const CalculationMatrix& B,
	const CalculationMatrix& C,
	Matrix2D<double>& psd,
	int num_substeps,
	double* lu_cache,
	long* index_cache,
	bool recompute_lu
);
#else
void Lapack(
	const CalculationMatrix& A,
	const CalculationMatrix& B,
	const CalculationMatrix& C,
	Matrix2D<double>& psd,
	int num_substeps
);
#endif


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
bool tridag(double a[], double b[], double c[], double r[], double u[], size_t n);

/// Lapack function for matrix inversion declaration
extern "C" {
	//extern void dgesv_(int *,int *,double *,int *,int*,double *,int*,int*);
	//extern void sgesv_(int *,int *,float *,int *,int*, float *,int*,int*);
	extern void dgbtrs_(char const* trans, long* m, long* kl, long* ku, long* nrhs, double* ab, long* ldab, long* ipiv, double* b, long* ldb, long* info);
	extern void dgbtrf_(long* m, long* n, long* kl, long* ku, double* ab , long* ldab, long* ipiv, long* info );
}

#endif

