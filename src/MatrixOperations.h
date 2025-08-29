/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * @file MatrixOperations.h
 * @brief Matrix operations and linear algebra functions
 * 
 * This file provides comprehensive matrix operations including concatenation,
 * arithmetic operations, comparison operators, logical operations, and
 * linear algebra routines using LAPACK for high performance.
 */

#pragma once

#include "Matrix.h"
#include <bitset>
#include <vector>

/**
 * @brief Concatenate two 1D matrices along the first dimension
 * @param a First matrix to concatenate
 * @param b Second matrix to concatenate
 * @return Concatenated matrix with elements from a followed by elements from b
 */
Matrix1D<double> cat(const Matrix1D<double>& a, const Matrix1D<double>& b);

/**
 * @brief Concatenate two 2D matrices along the first dimension (rows)
 * @param a First matrix to concatenate
 * @param b Second matrix to concatenate (must have same number of columns)
 * @return Concatenated matrix with rows from a followed by rows from b
 */
Matrix2D<double> cat(const Matrix2D<double>& a, const Matrix2D<double>& b);

/**
 * @brief Concatenate two 3D matrices along the first dimension
 * @param a First matrix to concatenate
 * @param b Second matrix to concatenate (must have same size in dimensions 2 and 3)
 * @return Concatenated matrix with slices from a followed by slices from b
 */
Matrix3D<double> cat(const Matrix3D<double>& a, const Matrix3D<double>& b);
/**
 * @brief Compute matrix transpose
 * @param A Input matrix to transpose
 * @return Transposed matrix where A[i][j] becomes result[j][i]
 */
Matrix2D<double> transpose(const Matrix2D<double>& A);

/**
 * @brief Compute element-wise base-10 logarithm of 1D matrix
 * @param A Input matrix
 * @return Matrix with log10 of each element (minimum value clamped to -31)
 */
Matrix1D<double> log10(const Matrix1D<double>& A);

/**
 * @brief Compute element-wise base-10 logarithm of 2D matrix
 * @param A Input matrix
 * @return Matrix with log10 of each element (minimum value clamped to -31)
 */
Matrix2D<double> log10(const Matrix2D<double>& A);

/**
 * @brief Compute element-wise power operation: base^A[i]
 * @param base Base value for power operation
 * @param A Matrix containing exponents
 * @return Matrix where result[i] = base^A[i]
 */
Matrix1D<double> pow(double base, const Matrix1D<double>& A);

/**
 * @brief Create diagonal matrix from vector
 * @param d Vector containing diagonal elements
 * @return Square matrix with d[i] on diagonal and zeros elsewhere
 */
Matrix2D<double> diag(const Matrix1D<double>& d); 
/**
 * @brief Element-wise equality comparison: A[i] == number
 * @param A Input matrix
 * @param number Scalar value to compare against
 * @return Boolean matrix with comparison results
 */
Matrix1D<bool> operator==(const Matrix1D<double>& A, double number);

/**
 * @brief Element-wise greater than comparison: A[i] > number
 * @param A Input matrix
 * @param number Scalar value to compare against
 * @return Boolean matrix with comparison results
 */
Matrix1D<bool> operator>(const Matrix1D<double>& A, double number);

/**
 * @brief Element-wise greater than or equal comparison: A[i] >= number
 * @param A Input matrix
 * @param number Scalar value to compare against
 * @return Boolean matrix with comparison results
 */
Matrix1D<bool> operator>=(const Matrix1D<double>& A, double number);

/**
 * @brief Element-wise less than comparison: A[i] < number
 * @param A Input matrix
 * @param number Scalar value to compare against
 * @return Boolean matrix with comparison results
 */
Matrix1D<bool> operator<(const Matrix1D<double>& A, double number);

/**
 * @brief Element-wise less than or equal comparison: A[i] <= number
 * @param A Input matrix
 * @param number Scalar value to compare against
 * @return Boolean matrix with comparison results
 */
Matrix1D<bool> operator<=(const Matrix1D<double>& A, double number);

/**
 * @brief Element-wise equality comparison: number == A[i]
 * @param number Scalar value to compare against
 * @param A Input matrix
 * @return Boolean matrix with comparison results
 */
Matrix1D<bool> operator==(double number, const Matrix1D<double>& A);

/**
 * @brief Element-wise greater than comparison: number > A[i]
 * @param number Scalar value to compare against
 * @param A Input matrix
 * @return Boolean matrix with comparison results
 */
Matrix1D<bool> operator>(double number, const Matrix1D<double>& A);

/**
 * @brief Element-wise greater than or equal comparison: number >= A[i]
 * @param number Scalar value to compare against
 * @param A Input matrix
 * @return Boolean matrix with comparison results
 */
Matrix1D<bool> operator>=(double number, const Matrix1D<double>& A);

/**
 * @brief Element-wise less than comparison: number < A[i]
 * @param number Scalar value to compare against
 * @param A Input matrix
 * @return Boolean matrix with comparison results
 */
Matrix1D<bool> operator<(double number, const Matrix1D<double>& A);

/**
 * @brief Element-wise less than or equal comparison: number <= A[i]
 * @param number Scalar value to compare against
 * @param A Input matrix
 * @return Boolean matrix with comparison results
 */
Matrix1D<bool> operator<=(double number, const Matrix1D<double>& A);         
/**
 * @brief Element-wise logical AND operation: a[i] && b[i]
 * @param a First boolean matrix
 * @param b Second boolean matrix (must have same size as a)
 * @return Boolean matrix with logical AND results
 */
Matrix1D<bool> operator&&(const Matrix1D<bool>& a, const Matrix1D<bool>& b);

/**
 * @brief Element-wise logical OR operation: a[i] || b[i]
 * @param a First boolean matrix
 * @param b Second boolean matrix (must have same size as a)
 * @return Boolean matrix with logical OR results
 */
Matrix1D<bool> operator||(const Matrix1D<bool>& a, const Matrix1D<bool>& b); 
/**
 * @brief Element-wise equality comparison for 2D matrix: A[i][j] == number
 * @param A Input matrix
 * @param number Scalar value to compare against
 * @return Boolean matrix with comparison results
 */
Matrix2D<bool> operator==(const Matrix2D<double>& A, double number);

/**
 * @brief Element-wise greater than comparison for 2D matrix: A[i][j] > number
 * @param A Input matrix
 * @param number Scalar value to compare against
 * @return Boolean matrix with comparison results
 */
Matrix2D<bool> operator>(const Matrix2D<double>& A, double number);

/**
 * @brief Element-wise greater than or equal comparison for 2D matrix: A[i][j] >= number
 * @param A Input matrix
 * @param number Scalar value to compare against
 * @return Boolean matrix with comparison results
 */
Matrix2D<bool> operator>=(const Matrix2D<double>& A, double number);

/**
 * @brief Element-wise less than comparison for 2D matrix: A[i][j] < number
 * @param A Input matrix
 * @param number Scalar value to compare against
 * @return Boolean matrix with comparison results
 */
Matrix2D<bool> operator<(const Matrix2D<double>& A, double number);

/**
 * @brief Element-wise less than or equal comparison for 2D matrix: A[i][j] <= number
 * @param A Input matrix
 * @param number Scalar value to compare against
 * @return Boolean matrix with comparison results
 */
Matrix2D<bool> operator<=(const Matrix2D<double>& A, double number);

/**
 * @brief Element-wise equality comparison for 2D matrix: number == A[i][j]
 * @param number Scalar value to compare against
 * @param A Input matrix
 * @return Boolean matrix with comparison results
 */
Matrix2D<bool> operator==(double number, const Matrix2D<double>& A);

/**
 * @brief Element-wise greater than comparison for 2D matrix: number > A[i][j]
 * @param number Scalar value to compare against
 * @param A Input matrix
 * @return Boolean matrix with comparison results
 */
Matrix2D<bool> operator>(double number, const Matrix2D<double>& A);

/**
 * @brief Element-wise greater than or equal comparison for 2D matrix: number >= A[i][j]
 * @param number Scalar value to compare against
 * @param A Input matrix
 * @return Boolean matrix with comparison results
 */
Matrix2D<bool> operator>=(double number, const Matrix2D<double>& A);

/**
 * @brief Element-wise less than comparison for 2D matrix: number < A[i][j]
 * @param number Scalar value to compare against
 * @param A Input matrix
 * @return Boolean matrix with comparison results
 */
Matrix2D<bool> operator<(double number, const Matrix2D<double>& A);

/**
 * @brief Element-wise less than or equal comparison for 2D matrix: number <= A[i][j]
 * @param number Scalar value to compare against
 * @param A Input matrix
 * @return Boolean matrix with comparison results
 */
Matrix2D<bool> operator<=(double number, const Matrix2D<double>& A);

/**
 * @brief Element-wise greater than or equal comparison between matrices: A[i][j] >= B[i][j]
 * @param A First input matrix
 * @param B Second input matrix (must have same dimensions as A)
 * @return Boolean matrix with comparison results
 */
Matrix2D<bool> operator>=(const Matrix2D<double>& A, const Matrix2D<double>& B);

/**
 * @brief Element-wise logical AND operation for 2D matrices: A[i][j] && B[i][j]
 * @param A First boolean matrix
 * @param B Second boolean matrix (must have same dimensions as A)
 * @return Boolean matrix with logical AND results
 */
Matrix2D<bool> operator&&(const Matrix2D<bool>& A, const Matrix2D<bool>& B); 
/**
 * @brief Element-wise subtraction for 1D matrices: a[i] - b[i]
 * @param a First matrix (minuend)
 * @param b Second matrix (subtrahend, must have same size as a)
 * @return Matrix with element-wise differences
 */
Matrix1D<double> operator-(const Matrix1D<double>& a, const Matrix1D<double>& b); 
/**
 * @brief Compute the matrix vector product of matrix A and vector b using Lapack
 * 
 * @param A left factor of product
 * @param b right factor of product
 * @return Product of A and B 
 */
Matrix1D<double> operator*(const Matrix2D<double>& A, const Matrix1D<double>& b);
/**
 * @brief Compute the matrix product of A and B using Lapack
 * 
 * @param A left factor of product
 * @param B right factor of product
 * @return Product of A and B 
 */ 
Matrix2D<double> operator*(const Matrix2D<double>& A, const Matrix2D<double>& B); 
/**
 * @brief Matrix exponentiation: compute A^power using efficient binary exponentiation
 * @param A Square matrix to raise to power
 * @param power Integer exponent (must be positive)
 * @return Matrix A raised to the specified power
 */
Matrix2D<double> operator^(Matrix2D<double> A, const int power);

/**
 * @brief Element-wise addition for 2D matrices: a[i][j] + b[i][j]
 * @param a First matrix
 * @param b Second matrix (must have same dimensions as a)
 * @return Matrix with element-wise sums
 */
Matrix2D<double> operator+(const Matrix2D<double>& a, const Matrix2D<double>& b);

/**
 * @brief Element-wise subtraction for 2D matrices: a[i][j] - b[i][j]
 * @param a First matrix (minuend)
 * @param b Second matrix (subtrahend, must have same dimensions as a)
 * @return Matrix with element-wise differences
 */
Matrix2D<double> operator-(const Matrix2D<double>& a, const Matrix2D<double>& b);

/**
 * @brief Compute matrix inverse using LAPACK
 * @param A Square matrix to invert (modified in-place)
 * @return Inverse of matrix A
 */
Matrix2D<double> inv(Matrix2D<double> A);

/**
 * @brief Solve transposed linear system: A^T * X^T = B^T (equivalent to X * A = B)
 * @param A Coefficient matrix
 * @param B Right-hand side matrix (modified in-place to contain solution)
 * @return Reference to B containing the solution matrix X
 */
Matrix2D<double>& trans_solve(const Matrix2D<double>& A, Matrix2D<double>& B);

/**
 * @brief Convert 2D matrix to 1D matrix by flattening row-wise
 * @param A Input 2D matrix
 * @return 1D matrix with elements in row-major order
 */
Matrix1D<double> toMatrix1D(const Matrix2D<double>& A);

/**
 * @brief Convert 1D matrix to 2D matrix with specified dimensions
 * @param A Input 1D matrix
 * @param size_q1 Number of rows in output matrix
 * @param size_q2 Number of columns in output matrix
 * @return 2D matrix with elements arranged in row-major order
 */
Matrix2D<double> toMatrix2D(const Matrix1D<double>& A, size_t size_q1, size_t size_q2);

/**
 * @brief Compute matrix product A * B^T using LAPACK
 * @param A First matrix
 * @param B Second matrix (transposed during multiplication)
 * @return Product A * B^T
 */
Matrix2D<double> abtrans(const Matrix2D<double>& A, const Matrix2D<double>& B);