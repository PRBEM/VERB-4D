/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/** Matrix 1D, 2D, 3D and 4D and operations with them
 *
 * \file Matrix.h
 *
 * File has 1D-class, 2D-class 3D-class and 4D-class of matrixes and various functions to work with them.
 *
 * \author Developed under supervision of the PI Yuri Shprits
 *
 * \brief Matrix 1D, 2D, 3D and 4D and operations with them
 *
 */
#pragma once

#include <assert.h>
#include <string>
#include <fstream>
#include <memory.h>
#include <cmath>
#include <map>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <assert.h>
#include <ctime>
#include <algorithm>
#include "Parameters.h"

// Matlab library which will have to be linked at compile time
// Usually found in matlabroot/extern/include where matlabroot is the result of typing the matlabroot command into matlab
// Used for reading/writing functions for .mat files
#if (MATLAB_CAPABLE)
#include <mat.h>
#endif

// Stupid hook to make it works in Visual Studio with Maltab 2013-2015 libs
//#define printf printf

// Another way is to not include mex.h into the matrix.h
//#include <mex.h>

/**
 * Matrix 1D class
 *
 * Matrixes and operations.
 *
 * @brief A matrix of 1 dimensions that is immutable
 */
template <typename T>
class Matrix1D {
private:
	size_t num_elements = 0;
public:
	T *plane_array = nullptr;
	T *matrix_array = nullptr;						///< Array to keep the values
	bool initialized = false;						///< Flag, equal true if initialized
	size_t size_q1 = 0;								///< size x
	std::string name;							///< name of the Matrix

	// constructors and destructors
	Matrix1D() = default;
	Matrix1D( size_t size_q1);
	Matrix1D( const Matrix1D<T> &M );
    Matrix1D<T>(size_t size_q1, const T* data);     // pybind is using this constructor
	~Matrix1D();
	
	virtual void AllocateMemory( size_t size_q1 );

	// Operators
	T& operator[](size_t i);													///< Return the i-th value of matrix
	T& operator[](size_t i) const;												///< const function to return the i-th value of matrix
	T& operator()(int q1);													///< Return the x-th value of matrix
	T& Value (int q1) { return operator()(q1); }								///< Return the (x,y) value of matrix
	Matrix1D<T>& MatrixArray () { return *this; }							///< Return pointer to the instance of the class.
	T* MatrixArrayPointer () { return matrix_array; }							///< Return pointer to the instance of the class.

	// unary
	const Matrix1D& operator+() const { return *this; }						 ///< Return the matrix
	const Matrix1D operator-() const { return ((*this)*(-1)); } 				///< Return the matrix with all values multiplied by -1

	// The following operators modify the matrix they applied to
	Matrix1D& operator= (const Matrix1D<T> &M);
	Matrix1D& operator= (const T val);

	// \todo Some of the matrix operators still need to be implimented
	// I didn't have time yet to write these functions - these are matrix opearations
	//Matrix1D& operator*= (const Matrix1D<T> &M); 					// reserved for something good
	//Matrix1D& operator/= (const Matrix1D<T> &M); 					// reserved for something good
	template <class S>
	Matrix1D& operator+= (const Matrix1D<S> &M);
	Matrix1D& operator-= (const Matrix1D<T> &M);
	Matrix1D& operator*= (const T Val);
	Matrix1D& operator/= (const T Val);
	//Matrix1D& operator+= (const T Val); 								///< Add the Val to each matrix element, stores result in the matrix it's applied to
	//Matrix1D& operator-= (const T Val); 								///< Substract the Val from each matrix element, stores result in the matrix it's applied to

	Matrix1D& times_equal (const Matrix1D<T> &M); 					///< Arraywise multiplication (A.*B), stores result in the matrix it's applied to
	//Matrix1D& divide_equal (const Matrix1D<T> &M); 					///< Arraywise division (A./B), stores result in the matrix it's applied to

	// The following operators save the result to a new matrix
	//Matrix1D operator* (const Matrix1D<T> &M) const; 				// reserved for something good
	//Matrix1D operator/ (const Matrix1D<T> &M) const;  				// reserved for something good
	//Matrix1D operator+ (const Matrix1D<T> &M) const;
	//Matrix1D operator- (const Matrix1D<T> &M) const;
	Matrix1D operator* (const T Val) const;
	Matrix1D operator/ (const T Val) const;
	//Matrix1D operator+ (const T Val) const; 							///< Add the Val to each matrix element, stores result in a new matrix
	//Matrix1D operator- (const T Val) const; 							///< Substract the Val from each matrix element, stores result in a new matrix

	Matrix1D times (const Matrix1D<T> &M) const; 					///< Arraywise multiplication (A.*B), stores result in a new matrix
	Matrix1D divide (const Matrix1D<T> &M) const; 					///< Arraywise division (A./B), stores result in a new matrix
	Matrix1D sqrt () const; 	                    				///< Arraywise square root sqrt(A), stores result in a new matrix
	Matrix1D<bool> is_finite () const; 							    ///< Bool matrix reflecting if a 'real' element is present (not NAN) 

	T dot (const Matrix1D<T> &M) const; 					///< Dot product
	T norm () const; 					///< Norm

	//T max();
	//T maxabs();

	// writeToFileting
	virtual void writeToFile(const std::string& filename) const;
	virtual void writeToFile(const std::string& filename, const Matrix1D<T>& grid_q1) const;
	virtual void readFromFile(const std::string& filename);
	virtual void readFromFile(const std::string& filename, const Matrix1D<T>& grid_q1);

	//ADDED
	virtual void readFromMatlabFile(const std::string& filename, int columnNumber = 1);
	virtual void readFromMatlabFile(const std::string& filename, const Matrix1D<T>& grid_q1);

	T min() const;
	T max() const;
	T maxabs() const;
	Matrix1D<T> abs() const;
	Matrix1D<T>& max_of(T val);
};

/**
 * Matrix 2D class
 *
 * Matrixes and operations.
 *
 * @brief A matrix of 2 dimensions that is immutable
 */
template <typename T> class Matrix2D {
private:
	/// Matrix array (array of links to other arrays). Final links pointed to the memory addresses of the plane array. Matrix[x][y] can be used.
	/// Also, all rows saved in the memory one after another as a big array. So Matrix[x+x_size*y] can be also used.
	T **matrix_array = nullptr;
	size_t num_elements = 0;
public:
	// const static int N_of_dimentions2 = 2;									///< Not used anywhere

	T *plane_array = nullptr;
	bool initialized = false;														///< Flag, equal true if initialized
	size_t size_q1 = 0;															///< size x
	size_t size_q2 = 0;															 ///< size_y
	std::string name;															///< name of the Matrix

	// Constructors and destructors
	Matrix2D() = default;
	Matrix2D( const Matrix2D<T> &M );
	Matrix2D( size_t size_q1, size_t size_q2 );
    Matrix2D(const ssize_t* sizes, const T* data); 	// pybind is using this constructor
	~Matrix2D();

	virtual void AllocateMemory(size_t size_q1, size_t size_q2);

	// Operators

	T* operator[](size_t i) { return matrix_array[i]; }						///< Return the i-th pointer to 1d-array. Next [j] can be applied, so we have regular [i][j].
	T* operator[](size_t i) const { return matrix_array[i]; }				///< const function to return the i-th pointer to 1d-array. Next [j] can be applied, so we have regular [i][j].
	T& operator()(int q1, int q2) { return matrix_array[0][q1*size_q2 + q2]; }	///< Return the (x,y)-th value of matrix
	T& Value (int q1, int q2) { return operator()(q1, q2); }				///< Return the (x,y) value of matrix
	Matrix2D<T>& MatrixArray () { return *this; }					///< Return pointer to the instance of the class.

	// unary
	const Matrix2D& operator+() const { return *this; } ///< unary : return the matrix
	const Matrix2D operator-() const { return ((*this)*(-1)); } ///< unary : return the matrix with all values multiplied by -1

	// The following operators modify the matrix they applied to
	Matrix2D& operator= (const Matrix2D<T> &M);
	Matrix2D& operator= (const T val);

	//Matrix2D& operator*= (const Matrix2D<T> &M); 					// reserved for something good
	//Matrix2D& operator/= (const Matrix2D<T> &M); 					// reserved for something good
	template <class S>
	Matrix2D& operator+= (const Matrix2D<S> &M);
	Matrix2D& operator-= (const Matrix2D<T> &M);
	Matrix2D& operator*= (const T Val);
	//Matrix2D& operator/= (const T Val);
	//Matrix2D& operator+= (const T Val); 								///< Add the Val to each matrix element, stores result in the matrix it's applied to
	//Matrix2D& operator-= (const T Val); 								///< Substract the Val from each matrix element, stores result in the matrix it's applied to

	//Matrix2D& times_equal (const Matrix2D<T> &M); 					///< Arraywise multiplication (A.*B), stores result in the matrix it's applied to
	//Matrix2D& divide_equal (const Matrix2D<T> &M); 					///< Arraywise division (A./B), stores result in the matrix it's applied to

	// The following operators save the result to a new matrix
	//Matrix2D operator* (const Matrix2D<T> &M) const; 				// Matrix multiplication
	//Matrix2D operator/ (const Matrix2D<T> &M) const;  				// reserved for something good
	//Matrix2D operator+ (const Matrix2D<T> &M) const;
	//Matrix2D operator- (const Matrix2D<T> &M) const;
	Matrix2D operator* (const T Val) const;
	Matrix2D operator/ (const T Val) const;
	//Matrix2D operator+ (const T Val) const; 							///< Add the Val to each matrix element, stores result in a new matrix
	//Matrix2D operator- (const T Val) const; 							///< Substract the Val from each matrix element, stores result in a new matrix

	Matrix2D times (const Matrix2D<T> &M) const; 					///< Arraywise multiplication (A.*B), stores result in a new matrix
	Matrix2D divide (const Matrix2D<T> &M) const; 					///< Arraywise division (A./B), stores result in a new matrix
	Matrix2D<bool> is_finite () const; 							    ///< Bool matrix reflecting if a 'real' element is present (not NAN) 

	// Return corresponding index of 1d array
	size_t index1d(size_t q1, size_t q2) const;

	T max() const;
	T maxabs() const;
	T min() const;
	Matrix2D<T> abs() const;
	Matrix2D<T> exp(double multiplicator=1) const;
	// It returns maximum between values from class psd2DSlice and argument (VC::zero_f in that case)
	Matrix2D<T>& max_of(T val);

	// writeToFileting
	virtual void writeToFile(const std::string& filename) const;
	virtual void writeToFile(const std::string& filename, const Matrix2D<T>& grid_q1, const Matrix2D<T>& grid_q2) const;
	virtual void readFromFile(const std::string& filename, int column = 1);
	virtual void readFromFile(const std::string& filename, const Matrix2D<T>& grid_q1, const Matrix2D<T>& grid_q2);

	//ADDED
	virtual void readFromMatlabFile(const std::string& filename, int column = 1);
	virtual void readFromMatlabFile(const std::string& filename, const Matrix2D<T>& grid_q1, const Matrix2D<T>& grid_q2);

	void writeToBinaryFile(const std::string& filename) const;

	// slices - get 1D slice from 2D array
	Matrix1D<T> xSlice(size_t p_q1) const;
	void xSlice(Matrix1D<T>& out, size_t p_q1) const;
	
	Matrix1D<T> ySlice(size_t p_q2) const;
	void ySlice(Matrix1D<T>& out, size_t p_q2) const;

};

/**
 * Matrix 3D class
 *
 * Matrixes and operations.
 *
 * @brief A matrix of 3 dimensions that is immutable
 */
template <typename T>
class Matrix3D {
private:
	/// Plane array of values. All rows saved in the memory one after anouther as a big array.
	/// Matrix array (array of links to other arrays). Final links pointed to the memory addresses of the plane array. Matrix[x][y][z] can be used.
	T ***matrix_array = nullptr;
	size_t num_elements = 0;
public:
	T *plane_array = nullptr;
	bool initialized = false;														///< Flag, equal true if initialized
	size_t size_q1 = 0;															///< size x
	size_t size_q2 = 0;															///< size y
	size_t size_q3 = 0;															///< size z
	std::string name;															///< name of the Matrix

	// constructors and destructors
	/// Default constructor. Do nothing.
	Matrix3D() = default;
	Matrix3D( const Matrix3D<T> &M );
	Matrix3D( size_t size_q1, size_t size_q2, size_t size_q3 );
	~Matrix3D();

	virtual void AllocateMemory(size_t size_q1, size_t size_q2, size_t size_q3);

	// Operators
	T** operator[] (size_t i); 											///< Return the i-th pointer to 2d-array. Next [j][k] can be applied, so we have regular [i][j][k].
	T** operator[] (size_t i) const { return matrix_array[i]; }			///< const function to return the i-th pointer to 2d-array. Next [j][k] can be applied, so we have regular [i][j][k].
	T& operator() (int q1, int q2, int q3); 							///< Return the (x,y,z) value of matrix
	T& Value (int q1, int q2, int q3) { return operator()(q1, q2, q3); }	///< Return the (x,y,z) value of matrix
	Matrix3D<T>& MatrixArray () { return *this; }					///< Return pointer to the instance of the class.

	// The following operators modify the matrix they applied to
	Matrix3D& operator= (const Matrix3D<T> &M);
//	Matrix3D& operator= (const Matrix2D<T> &M); 						// confusing function
	Matrix3D& operator= (const T Val);

	// unary
	const Matrix3D& operator+() const { return *this; } ///< unary : return the matrix
	const Matrix3D operator-() const { return ((*this)*(-1)); } ///< unary : return the matrix with all values multiplied by -1


	//Matrix3D& operator*= (const Matrix3D<T> &M); 					// reserved for something good
	//Matrix3D& operator/= (const Matrix3D<T> &M); 					// reserved for something good
	template <class S>
	Matrix3D& operator+= (const Matrix3D<S> &M);
	Matrix3D& operator-= (const Matrix3D<T> &M);
	Matrix3D& operator*= (const T Val);
	Matrix3D& operator/= (const T Val);
	Matrix3D& operator+= (const T Val); 								///< Add the Val to each matrix element, stores result in the matrix it's applied to
	Matrix3D& operator-= (const T Val); 								///< Substract the Val from each matrix element, stores result in the matrix it's applied to

	Matrix3D& times_equal (const Matrix3D<T> &M); 					///< Arraywise multiplication (A.*B), stores result in the matrix it's applied to
	Matrix3D& divide_equal (const Matrix3D<T> &M); 					///< Arraywise division (A./B), stores result in the matrix it's applied to

	// The following operators save the result to a new matrix
	//Matrix3D operator* (const Matrix3D<T> &M) const; 				// reserved for something good
	//Matrix3D operator/ (const Matrix3D<T> &M) const;  				// reserved for something good
	Matrix3D operator+ (const Matrix3D<T> &M) const;
	Matrix3D operator- (const Matrix3D<T> &M) const;
	Matrix3D operator* (const T Val) const;
	Matrix3D operator/ (const T Val) const;
	//Matrix3D operator+ (const T Val) const; 							///< Add the Val to each matrix element, stores result in a new matrix
	//Matrix3D operator- (const T Val) const; 							///< Substract the Val from each matrix element, stores result in a new matrix

	Matrix3D times (const Matrix3D<T> &M) const; 					///< Arraywise multiplication (A.*B), stores result in a new matrix
	Matrix3D divide (const Matrix3D<T> &M) const; 					///< Arraywise division (A./B), stores result in a new matrix
	Matrix3D sqrt () const; 	                    				///< Arraywise square root sqrt(A), stores result in a new matrix
	Matrix3D<bool> is_finite () const; 							    ///< Bool matrix reflecting if a 'real' element is present (not NAN) 

	// Saving (loading) of a matrix into (from) file
	virtual void writeToFile(const std::string& filename, const std::string& info = "") const; 										///< Save matrix to a file
	virtual void writeToFile(const std::string& filename, const Matrix3D<T>& grid_q1, const Matrix3D<T>& grid_q2, const Matrix3D<T>& grid_q3) const; ///< Save matrix to a file, including grid
	virtual void readFromFile(const std::string& filename, int column = 1);  									///< Load matrix from a file
	virtual void readFromFile(const std::string& filename, const Matrix3D<T>& grid_q1, const Matrix3D<T>& grid_q2, const Matrix3D<T>& grid_q3); ///< Load matrix to a file

	// ADDED
	virtual void readFromMatlabFile(const std::string& filename, int column = 1);
	virtual void readFromMatlabFile(const std::string& filename, const Matrix3D<T>& grid_q1, const Matrix3D<T>& grid_q2, const Matrix3D<T>& grid_q3);

    virtual void writeToBinaryFile(const std::string& filename) const;
    virtual void readFromBinaryFile(const std::string& filename);

    virtual void writeToAnyFile(const std::string& filename, const IOMethod& io_method, const std::string& info) const;
    virtual void readFromAnyFile(const std::string& filename, const IOMethod& io_method);
    virtual void readFromAnyFile(const std::string& filename, const IOMethod& io_method, const Matrix3D<T>& grid_q1, const Matrix3D<T>& grid_q2, const Matrix3D<T>& grid_q3);

	// Some other stuff
	std::string change_ind;														///< Variables useful for tracking of changes (time of change can be stored here)

	size_t index1d(size_t q1, size_t q2, size_t q3);								///< Returns index of the element (x,y,z) in 1d array

	T min() const;
	T max() const;
	T maxabs() const;
	Matrix3D<T> abs() const;
	Matrix3D<T> exp(double multiplicator=1) const;
	Matrix3D<T>& max_of(T val);

	// slices - get 2D slice from 3D array
	Matrix2D<T> xSlice(size_t p_q1) const;
	void xSlice(Matrix2D<T>& out, size_t p_q1) const;

	Matrix2D<T> ySlice(size_t p_q2) const;
	void ySlice(Matrix2D<T>& out, size_t p_q2) const;

	Matrix2D<T> zSlice(size_t p_q3) const;
	void zSlice(Matrix2D<T>& out, size_t p_q3) const;

	// slices - get 1D slice from 3D array
	Matrix1D<T> xySlice(size_t p_i1, size_t p_i2) const;
	void xySlice(Matrix1D<T>& out, size_t p_i1, size_t p_i2) const;

	Matrix1D<T> yzSlice(size_t p_i2, size_t p_i3) const;
	void yzSlice(Matrix1D<T>& out, size_t p_i2, size_t p_i3) const;

	Matrix1D<T> xzSlice(size_t p_i1, size_t p_i3) const;
	void xzSlice(Matrix1D<T>& out, size_t p_i1, size_t p_i3) const;


};


/**
 * Matrix 4D class
 *
 * Matrixes and operations.
 *
 * @brief A matrix of 4 dimensions that is immutable
 */
template <typename T> class Matrix4D {
private:
	/// Plane array of values. All rows saved in the memory one after anouther as a big array.
	/// Matrix array (array of links to other arrays). Final links pointed to the memory addresses of the plane array. Matrix[w][x][y][z] can be used.
	T ****matrix_array = nullptr;
	size_t num_elements = 0;
public:
	T *plane_array = nullptr;
	bool initialized = false;														///< Flag, equal true if initialized
	size_t size_w = 0;																///< size w
	size_t size_x = 0;																///< size x
	size_t size_y = 0; 															///< size y
	size_t size_z = 0;																///< size z
	std::string name;															///< name of the Matrix

	// constructors and destructors
	/// Default constructor. Do nothing.
	Matrix4D() = default;
	Matrix4D( const Matrix4D<T> &M );
	Matrix4D( size_t size_w, size_t size_x, size_t size_y, size_t size_z );
	~Matrix4D();

	virtual void AllocateMemory(size_t size_w, size_t size_x, size_t size_y, size_t size_z);

	// Operators
	T*** operator[] (size_t i); 											///< Return the i-th pointer to 3d-array. Next [j][k][l] can be applied, so we have regular [i][j][k][l].
	T*** operator[] (size_t i) const { return matrix_array[i]; }            ///< const function to return the i-th pointer to 3d-array. Next [j][k][l] can be applied, so we have regular [i][j][k][l].
	T& operator() (int w, int x, int y, int z); 							///< Return the (w,x,y,z) value of matrix
	T& Value (int w, int x, int y, int z) { return operator()(w, x, y, z); }	///< Return the (w,x,y,z) value of matrix
	Matrix4D<T>& MatrixArray () { return *this; }					///< Return pointer to the instance of the class.

	// The following operators modify the matrix they applied to
	Matrix4D& operator= (const Matrix4D<T> &M);
	Matrix4D& operator= (const T Val);

	// unary
	const Matrix4D& operator+() const { return *this;}  ///< unary : return the matrix
	const Matrix4D operator-() const { return ((*this)*(-1)); } ///< unary : return the matrix with all values multiplied by -1

	//Matrix4D& operator*= (const Matrix4D<T> &M); 					// reserved for something good
	//Matrix4D& operator/= (const Matrix4D<T> &M); 					// reserved for something good
	template <class S>
	Matrix4D& operator+= (const Matrix4D<S> &M);
	Matrix4D& operator-= (const Matrix4D<T> &M);
	Matrix4D& operator*= (const T Val);
	Matrix4D& operator/= (const T Val);
	Matrix4D& operator+= (const T Val); 								///< Add Val to each matrix element, stores result in the matrix it's applied to
	Matrix4D& operator-= (const T Val); 								///< Substract Val from each element in this matrix, stores result in the matrix it's applied to

	Matrix4D& times_equal (const Matrix4D<T> &M); 					///< Arraywise multiplication (A.*B), stores result in the matrix it's applied to
	Matrix4D& divide_equal (const Matrix4D<T> &M); 					///< Arraywise division (A./B), stores result in the matrix it's applied to

	// The following operators save the result to a new matrix
	//Matrix4D operator* (const Matrix4D<T> &M) const; 				// reserved for something good
	//Matrix4D operator/ (const Matrix4D<T> &M) const;  				// reserved for something good
	Matrix4D operator+ (const Matrix4D<T> &M) const;
	Matrix4D operator- (const Matrix4D<T> &M) const;
	Matrix4D operator* (const T Val) const;
	Matrix4D operator/ (const T Val) const;
	//Matrix4D operator+ (const T Val) const; 							///< Add the Val to each matrix element, stores result in a new matrix
	//Matrix4D operator- (const T Val) const; 							///< Substract the Val from each matrix element, stores result in a new matrix

	Matrix4D times (const Matrix4D<T> &M) const; 					///< Arraywise multiplication (A.*B), stores result in a new matrix
	Matrix4D divide (const Matrix4D<T> &M) const; 					///< Arraywise division (A./B), stores result in a new matrix
	Matrix4D<bool> is_finite () const; 							    ///< Bool matrix reflecting if a 'real' element is present (not NAN) 

	// Saving (loading) of a matrix into (from) file
	virtual void writeToFile(const std::string& filename, const std::string& info = "") const;								///< Save matrix to a file
	virtual void writeToFile(const std::string& filename, const Matrix4D<T>& grid_w, const Matrix4D<T>& grid_x, const Matrix4D<T>& grid_y, const Matrix4D<T>& grid_z) const; ///< Save matrix to a file, including grid
	virtual void readFromFile(const std::string& filename, int column = 1);  									///< Load matrix from a file
	virtual void readFromFile(const std::string& filename, const Matrix4D<T>& grid_w, const Matrix4D<T>& grid_x, const Matrix4D<T>& grid_y, const Matrix4D<T>& grid_z); ///< Load matrix to a file

	// ADDED
#if (MATLAB_CAPABLE)
	virtual mxArray* createStructMatrix(const std::string& info = "") const; ///< Package matrix into matlab variable
#endif

	virtual void writeToMatlabFile(const std::string& filename, const std::string& info = "") const;									///< Save matrix to a file
	virtual void writeToMatlabFile(const std::string& filename, const Matrix4D<T>& grid_w, const Matrix4D<T>& grid_x, const Matrix4D<T>& grid_y, const Matrix4D<T>& grid_z) const;
	virtual void readFromMatlabFile(const std::string& filename, int column = 1);
	virtual void readFromMatlabFile(const std::string& filename, const Matrix4D<T>& grid_w, const Matrix4D<T>& grid_x, const Matrix4D<T>& grid_y, const Matrix4D<T>& grid_z);

    virtual void writeToBinaryFile(const std::string& filename) const;
    virtual void readFromBinaryFile(const std::string& filename);

    virtual void writeToAnyFile(const std::string& filename, const IOMethod& io_method, const std::string& info) const;
    virtual void readFromAnyFile(const std::string& filename, const IOMethod& io_method);
    virtual void readFromAnyFile(const std::string& filename, const IOMethod& io_method, const Matrix4D<T>& grid_w, const Matrix4D<T>& grid_x, const Matrix4D<T>& grid_y, const Matrix4D<T>& grid_z);

	virtual std::string getExtByIoMethod(const IOMethod& io_method) const;
	virtual void writeToLstFile(const std::string& filename, const IOMethod& io_method, const std::string& info, const std::string& output_folder) const;

	// Some other stuff
	std::string change_ind;														///< Variables useful for tracking of changes (time of change can be stored here)

	size_t index1d(size_t w, size_t x, size_t y);								///< Returns index of the element (x,y,z) in 1d array

	T min() const;
	T max() const;
	T maxabs() const;
	Matrix4D<T> abs() const;
	Matrix4D<T>& max_of(T val);
	// slices - get 3D slice from 4D array
	Matrix3D<T> wSlice(size_t p_w) const;
	void wSlice(Matrix3D<T>& out, size_t p_w) const;

	Matrix3D<T> xSlice(size_t p_x) const;
	void xSlice(Matrix3D<T>& out, size_t p_x) const;

	Matrix3D<T> ySlice(size_t p_y) const;
	void ySlice(Matrix3D<T>& out, size_t p_y) const;

	Matrix3D<T> zSlice(size_t p_z) const;
	void zSlice(Matrix3D<T>& out, size_t p_z) const;

	// slices - get 2D slice from 4D array
	Matrix2D<T> wxSlice(size_t p_w, size_t p_x) const;
	void wxSlice(Matrix2D<T>& out, size_t p_w, size_t p_y) const;

	Matrix2D<T> wySlice(size_t p_w, size_t p_y) const;
	void wySlice(Matrix2D<T>& out, size_t p_w, size_t p_y) const;

	Matrix2D<T> wzSlice(size_t p_w, size_t p_z) const;
	void wzSlice(Matrix2D<T>& out, size_t p_w, size_t p_y) const;
	
	Matrix2D<T> xySlice(size_t p_w, size_t p_z) const;
	void xySlice(Matrix2D<T>& out, size_t p_x, size_t p_y) const;

	Matrix2D<T> xzSlice(size_t p_x, size_t p_z) const;
	void xzSlice(Matrix2D<T>& out, size_t p_x, size_t p_y) const;
	
	Matrix2D<T> yzSlice(size_t p_x, size_t p_z) const;
	void yzSlice(Matrix2D<T>& out, size_t p_y, size_t p_z) const;

	// slices - get 1D slice from 4D array
	Matrix1D<T> wxySlice(size_t p_w, size_t p_x, size_t p_y) const;
	void wxySlice(Matrix1D<T>& out, size_t p_w, size_t p_x, size_t p_z) const;
	
	Matrix1D<T> wxzSlice(size_t p_w, size_t p_x, size_t p_y) const;
	void wxzSlice(Matrix1D<T>& out, size_t p_w, size_t p_x, size_t p_z) const;

	Matrix1D<T> wyzSlice(size_t p_w, size_t p_y, size_t p_z) const;
	void wyzSlice(Matrix1D<T>& out, size_t p_w, size_t p_x, size_t p_z) const;

	Matrix1D<T> xyzSlice(size_t p_w, size_t p_y, size_t p_z) const;
	void xyzSlice(Matrix1D<T>& out, size_t p_x, size_t p_y, size_t p_z) const;
};


/** @class DiagMatrix
 * This method of storage for matrices is convenient for diagonal (spread) matrices.
 * Stored as map (diagonal number, 1d diagonal array)
 * The USED diagonals of the matrix are stored in 1d arrays.
 * CalculationMatrix inherits from this
 * @brief Diagonals of matrix stored as map (diagonal number, 1d diagonal array)
 *
 * @typedef DiagMatrix
 * This method of storage for matrices is convenient for diagonal (spread) matrices.
 * Stored as map (diagonal number, 1d diagonal array)
 * The USED diagonals of the matrix are stored in 1d arrays.
 * This typedef is used in CalculationMatrix
 * @brief Diagonals of matrix stored as map (diagonal number, 1d diagonal array)
 */
typedef std::map<int, Matrix1D<double>> DiagMatrix;

/** Model matrix (or related matrices)
 * It is based on DiagMatrix and have methods for conversion from 3D or 2D PSD (and related) arrays into 1d array of unknown elements
 *
 * @brief A matrix used primarily for diagonalized calculations
 */
class CalculationMatrix : public DiagMatrix {
public:

	bool initialized; ///< Variables used for denoting whether it is initialized or not
	int x_size; ///< size in x direction of matrix
	int y_size; ///< size in y direction of matrix
	int total_size; ///< total size of matrix
	// flag, if needs to be recalculated
	clock_t change_ind;						///< Variables useful for changes tracking (store here time when changed)

	// Constructors
	CalculationMatrix();
	// !!! CalculationMatrix(int x_size, int y_size = 0, int z_size = 0, int n_of_diags = 1);
	CalculationMatrix(int x_size, int y_size = 1, int z_size = 1, int n_of_diags = 1);

	// Initialization
	// !!! void Initialize(int x_size, int y_size = 0, int z_size = 0, int n_of_diags = 1);
	void Initialize(int x_size, int y_size = 1, int z_size = 1, int n_of_diags = 1);

	// Returns 1d index for multiple dimension array (2D or 3D)
	int index1d(int x, int y = 0, int z = 0) const;

	// Save to a file
	void writeToFile(const std::string& filename) const;

	// Operators
	/// FUNCTION NOT IMPLEMENTED
	Matrix1D<double> operator* (Matrix1D<double> &V) const;

};

/// FUNCTION NOT IMPLEMENTED
int index1d2(int x, int y = 0, int z = 0);
