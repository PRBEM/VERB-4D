#pragma once

#include "Matrix.h"

Matrix1D<double> cat(const Matrix1D<double>& a, const Matrix1D<double>& b);
Matrix2D<double> cat(const Matrix2D<double>& a, const Matrix2D<double>& b);
Matrix3D<double> cat(const Matrix3D<double>& a, const Matrix3D<double>& b);
Matrix2D<double> transpose(const Matrix2D<double>& A);
Matrix1D<double> log10(const Matrix1D<double>& A);
Matrix2D<double> log10(const Matrix2D<double>& A);
Matrix1D<double> pow(double base, const Matrix1D<double>& A);
Matrix2D<double> diag(const Matrix1D<double>& d); 
Matrix1D<bool> operator==(const Matrix1D<double>& A, double number);       
Matrix1D<bool> operator>(const Matrix1D<double>& A, double number);         
Matrix1D<bool> operator>=(const Matrix1D<double>& A, double number);         
Matrix1D<bool> operator<(const Matrix1D<double>& A, double number);         
Matrix1D<bool> operator<=(const Matrix1D<double>& A, double number);         
Matrix1D<bool> operator==(double number, const Matrix1D<double>& A);         
Matrix1D<bool> operator>(double number, const Matrix1D<double>& A);         
Matrix1D<bool> operator>=(double number, const Matrix1D<double>& A);         
Matrix1D<bool> operator<(double number, const Matrix1D<double>& A);         
Matrix1D<bool> operator<=(double number, const Matrix1D<double>& A);         
Matrix1D<bool> operator&&(const Matrix1D<bool>& a, const Matrix1D<bool>& b); 
Matrix1D<bool> operator||(const Matrix1D<bool>& a, const Matrix1D<bool>& b); 
Matrix2D<bool> operator==(const Matrix2D<double>& A, double number);         
Matrix2D<bool> operator>(const Matrix2D<double>& A, double number);         
Matrix2D<bool> operator>=(const Matrix2D<double>& A, double number);         
Matrix2D<bool> operator<(const Matrix2D<double>& A, double number);         
Matrix2D<bool> operator<=(const Matrix2D<double>& A, double number);         
Matrix2D<bool> operator==(double number, const Matrix2D<double>& A);         
Matrix2D<bool> operator>(double number, const Matrix2D<double>& A);         
Matrix2D<bool> operator>=(double number, const Matrix2D<double>& A);         
Matrix2D<bool> operator<(double number, const Matrix2D<double>& A);         
Matrix2D<bool> operator<=(double number, const Matrix2D<double>& A);         
Matrix2D<bool> operator&&(const Matrix2D<bool>& A, const Matrix2D<bool>& B); 
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
Matrix2D<double> operator+(const Matrix2D<double>& a, const Matrix2D<double>& b); 
Matrix2D<double> operator-(const Matrix2D<double>& a, const Matrix2D<double>& b); 
Matrix2D<double> inv(Matrix2D<double> A); 
Matrix2D<double>& trans_solve(const Matrix2D<double>& A, Matrix2D<double>& B); 
Matrix1D<double> toMatrix1D(const Matrix2D<double>& A); 
Matrix2D<double> toMatrix2D(const Matrix1D<double>& A, size_t size_q1, size_t size_q2);
Matrix2D<double> abtrans(const Matrix2D<double>& A, const Matrix2D<double>& B);