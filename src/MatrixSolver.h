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

/**
 * @brief Adds boundary condition constraints to a diagonal matrix
 * 
 * @details This function modifies the coefficient matrix to enforce boundary conditions
 * by setting appropriate matrix elements. It handles both Dirichlet (constant value)
 * and Neumann (constant derivative) boundary conditions.
 * 
 * **Mathematical Implementation:**
 * 
 * For **ConstantValue** (Dirichlet): \f$ f_i = f_{BC} \f$
 * - Matrix row becomes: [0, ..., 0, 1, 0, ..., 0] with 1 at position i
 * 
 * For **ConstantDerivative** (Neumann): \f$ \frac{\partial f}{\partial x}\big|_i = g_{BC} \f$
 * - Finite difference: \f$ \frac{f_{i+1} - f_i}{\Delta x} = g_{BC} \f$
 * - Matrix row becomes: [0, ..., 1/Δx, -1/Δx, 0, ...]
 * 
 * @param Matrix_A [in,out] Diagonal matrix to modify (coefficient matrix)
 * @param type [in] Type of boundary condition to apply
 * @param in [in] Matrix row index for the boundary point
 * @param id1 [in] Relative index offset to neighboring point
 * @param dh [in] Grid spacing for derivative calculations
 * 
 * @warning ConstantDerivative implementation only works correctly for zero derivative
 */
void AddBoundary(DiagMatrix &Matrix_A, BoundaryConditionType type, int in, int id1, double dh);

/**
 * @brief Adds boundary conditions for 1D problems
 * 
 * @details This function applies boundary conditions at the endpoints of a 1D domain
 * by modifying the coefficient matrices. It handles both lower and upper boundaries
 * and supports all boundary condition types.
 * 
 * **Boundary Detection:**
 * - Lower boundary: ix == 0
 * - Upper boundary: ix == x_size - 1
 * - Interior points: no boundary conditions applied
 * 
 * **Matrix Modifications:**
 * - **matr_A**: Coefficient matrix for implicit terms (modified for BC constraints)
 * - **matr_C**: Constant term vector (set to boundary values)
 * 
 * **Mathematical Implementation:**
 * 
 * For each boundary point, the corresponding matrix row is replaced with
 * the appropriate boundary condition constraint, and the RHS vector is
 * set to the boundary value.
 * 
 * @param matr_A [in,out] Coefficient matrix for implicit terms
 * @param matr_C [in,out] Constant term matrix (RHS contributions)
 * @param x [in] Spatial coordinate array (size x_size)
 * @param x_size [in] Number of grid points in x-direction
 * @param x_LBC [in] Lower boundary condition value
 * @param x_UBC [in] Upper boundary condition value
 * @param x_LBC_type [in] Type of lower boundary condition
 * @param x_UBC_type [in] Type of upper boundary condition
 * @param ix [in] Current grid point index to check for boundary
 * 
 * @return true on successful application of boundary conditions
 */
bool AddBoundaries_1D(
	CalculationMatrix &matr_A, CalculationMatrix &matr_C,
	const Matrix1D<double> &x,
	int x_size,
	double x_LBC, double x_UBC,
	BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
	int ix
);

/**
 * @brief Adds boundary conditions for 2D problems
 * 
 * @details This function applies boundary conditions at the boundaries of a 2D domain
 * by modifying the coefficient matrices. It handles boundaries in both x and y directions
 * and supports all boundary condition types.
 * 
 * **Boundary Detection:**
 * - X-boundaries: ix == 0 (lower) or ix == x_size-1 (upper)
 * - Y-boundaries: iy == 0 (lower) or iy == y_size-1 (upper)
 * - Corner points: handled by the first applicable boundary condition
 * - Interior points: no boundary conditions applied
 * 
 * **Matrix Modifications:**
 * Similar to 1D case, but with 2D indexing and boundary value arrays.
 * 
 * **Boundary Value Arrays:**
 * - x_LBC, x_UBC: Arrays of size y_size (values along y for x-boundaries)
 * - y_LBC, y_UBC: Arrays of size x_size (values along x for y-boundaries)
 * 
 * **Mathematical Implementation:**
 * 
 * For 2D problems, boundary conditions are applied along entire boundary edges.
 * The matrix row corresponding to each boundary point is replaced with the
 * appropriate constraint equation.
 * 
 * @param matr_A [in,out] Coefficient matrix for implicit terms
 * @param matr_C [in,out] Constant term matrix (RHS contributions)
 * @param x [in] X-coordinate grid (x_size × y_size)
 * @param y [in] Y-coordinate grid (x_size × y_size)
 * @param x_size [in] Number of grid points in x-direction
 * @param y_size [in] Number of grid points in y-direction
 * @param x_LBC [in] Lower x boundary values (size y_size)
 * @param x_UBC [in] Upper x boundary values (size y_size)
 * @param y_LBC [in] Lower y boundary values (size x_size)
 * @param y_UBC [in] Upper y boundary values (size x_size)
 * @param x_LBC_type [in] Type of lower x boundary condition
 * @param x_UBC_type [in] Type of upper x boundary condition
 * @param y_LBC_type [in] Type of lower y boundary condition
 * @param y_UBC_type [in] Type of upper y boundary condition
 * @param ix [in] Current x grid index
 * @param iy [in] Current y grid index
 * @param in [in] Linear matrix index for point (ix, iy)
 * 
 * @return true on successful application of boundary conditions
 */
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

