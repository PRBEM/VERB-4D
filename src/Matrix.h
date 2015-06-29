/** Matrix 1D, 2D and 3D and operations with them
 *
 * \file Matrix.h
 *
 * File has 1D-class, 2D-class and 3D-class of matrixes and variouse functions to work with them.
 *
 * \author Developed under supervision of the PI Yuri Shprits
 *
 */
#ifndef matrix_array_MATRIX_H
#define matrix_array_MATRIX_H

#include <assert.h>
#include <string>
#include <string.h>
#include <fstream>
#include <memory.h>
#include <math.h>
#include <map>
#include <stdlib.h>
#include <iostream>
#include <sstream>

using namespace std;

/**
 * Matrix 1D class
 *
 * Matrixes and operations.
 */
template <typename T>
class Matrix1D {
public:

	T *matrix_array;						///< Array to keep the values

	bool initialized;						///< Flag, equal true if initialized
	int size_q1;								///< size x
	string name;							///< name of the Matrix

	// constructors and destructors
	Matrix1D() { initialized = false; };
	Matrix1D( int size_q1 , string name = "f");
	Matrix1D( const Matrix1D<T> &M );
	~Matrix1D();

	virtual void AllocateMemory( int size_q1 );

	// Operators
	inline T& operator[](int i);													///< Return the i-th value of matrix
	inline T& operator[](int i) const;
	inline T& operator()(int q1);													///< Return the x-th value of matrix
	inline T& Value (int q1) { return operator()(q1); }								///< Return the (x,y) value of matrix
	inline Matrix1D<T>& MatrixArray () { return *this; }							///< Return pointer to the instance of the class.
	inline T* MatrixArrayPointer () { return matrix_array; }							///< Return pointer to the instance of the class.

	// unary
	inline const Matrix1D& operator+() const { return *this; }
	inline const Matrix1D operator-() const { return ((*this)*(-1)); }

	// The following operators modify the matrix they applied to
	inline Matrix1D& operator= (const Matrix1D<T> &M);
	inline Matrix1D& operator= (const T val);

	// \todo Some of the matrix operators still need to be implimented
	// I didn't have time yet to write these functions - these are matrix opearations
	//inline Matrix1D& operator*= (const Matrix1D<T> &M); 					// reserved for something good
	//inline Matrix1D& operator/= (const Matrix1D<T> &M); 					// reserved for something good
	inline Matrix1D& operator+= (const Matrix1D<T> &M);
	inline Matrix1D& operator-= (const Matrix1D<T> &M);
	inline Matrix1D& operator*= (const T Val);
	inline Matrix1D& operator/= (const T Val);
	inline Matrix1D& operator+= (const T Val); 								///< Add the Val to each matrix element, stores result in the matrix it's applied to
	inline Matrix1D& operator-= (const T Val); 								///< Substract the Val from each matrix element, stores result in the matrix it's applied to

	//inline Matrix1D& times_equal (const Matrix1D<T> &M); 					///< Arraywise multiplication (A.*B), stores result in the matrix it's applied to
	//inline Matrix1D& divide_equal (const Matrix1D<T> &M); 					///< Arraywise division (A./B), stores result in the matrix it's applied to

	// The following operators save the result to a new matrix
	//inline Matrix1D operator* (const Matrix1D<T> &M) const; 				// reserved for something good
	//inline Matrix1D operator/ (const Matrix1D<T> &M) const;  				// reserved for something good
	//inline Matrix1D operator+ (const Matrix1D<T> &M) const;
	//inline Matrix1D operator- (const Matrix1D<T> &M) const;
	inline Matrix1D operator* (const T Val) const;
	inline Matrix1D operator/ (const T Val) const;
	//inline Matrix1D operator+ (const T Val) const; 							///< Add the Val to each matrix element, stores result in a new matrix
	//inline Matrix1D operator- (const T Val) const; 							///< Substract the Val from each matrix element, stores result in a new matrix

	inline Matrix1D times (const Matrix1D<T> &M) const; 					///< Arraywise multiplication (A.*B), stores result in a new matrix
	inline Matrix1D divide (const Matrix1D<T> &M) const; 					///< Arraywise division (A./B), stores result in a new matrix

	inline T dot (const Matrix1D<T> &M) const; 					///< Dot product
	inline T norm () const; 					///< Norm

	//T max();
	//T maxabs();

	// writeToFileting
	virtual void writeToFile(string filename);
	virtual void writeToFile(string filename, Matrix1D<T> &grid_q1);
	virtual void readFromFile(string filename);
	virtual void readFromFile(string filename, const Matrix1D<T> grid_q1);

	T min();
	T max();
	T maxabs();
	Matrix1D<T> abs();

};

/**
 * Matrix 2D class
 *
 * Matrixes and operations.
 */
template <typename T> class Matrix2D {
private:
	/// Matrix array (array of links to other arrays). Final likns pointed to the memory addresses of the plane array. Matrix[x][y] can be used.
	/// Also, all rows saved in the memory one after anouther as a big array. So Matrix[x+x_size*y] can be also used.
	T **matrix_array;
public:
	const static int N_of_dimentions2 = 2;

	bool initialized;														///< Flag, equal true if initialized
	int size_q1, size_q2;														///< size x, size_y
	string name;															///< name of the Matrix

	// Constructors and destructors
	Matrix2D() { initialized = false; };
	Matrix2D( const Matrix2D<T> &M );
	Matrix2D( int size_q1, int size_q2 );
	~Matrix2D();

	virtual void AllocateMemory(int size_q1, int size_q2);

	// Operators

	inline T* operator[](int i) { return matrix_array[i]; }						///< Return the i-th pointer to 1d-array. Next [j] can be applied, so we have regular [i][j].
	inline T* operator[](int i) const { return matrix_array[i]; }
	inline T& operator()(int q1, int q2) { return matrix_array[0][q1*size_q2 + q2]; }	///< Return the (x,y)-th value of matrix
	inline T& Value (int q1, int q2) { return operator()(q1, q2); }				///< Return the (x,y) value of matrix
	inline Matrix2D<T>& MatrixArray () { return *this; }					///< Return pointer to the instance of the class.

	// unary
	inline const Matrix2D& operator+() const { return *this; }
	inline const Matrix2D operator-() const { return ((*this)*(-1)); }

	// The following operators modify the matrix they applied to
	inline Matrix2D& operator= (const Matrix2D<T> &M);
	inline Matrix2D& operator= (const T val);

	//inline Matrix2D& operator*= (const Matrix2D<T> &M); 					// reserved for something good
	//inline Matrix2D& operator/= (const Matrix2D<T> &M); 					// reserved for something good
	inline Matrix2D& operator+= (const Matrix2D<T> &M);
	inline Matrix2D& operator-= (const Matrix2D<T> &M);
	inline Matrix2D& operator*= (const T Val);
	//inline Matrix2D& operator/= (const T Val);
	//inline Matrix2D& operator+= (const T Val); 								///< Add the Val to each matrix element, stores result in the matrix it's applied to
	//inline Matrix2D& operator-= (const T Val); 								///< Substract the Val from each matrix element, stores result in the matrix it's applied to

	//inline Matrix2D& times_equal (const Matrix2D<T> &M); 					///< Arraywise multiplication (A.*B), stores result in the matrix it's applied to
	//inline Matrix2D& divide_equal (const Matrix2D<T> &M); 					///< Arraywise division (A./B), stores result in the matrix it's applied to

	// The following operators save the result to a new matrix
	//inline Matrix2D operator* (const Matrix2D<T> &M) const; 				// Matrix multiplication
	//inline Matrix2D operator/ (const Matrix2D<T> &M) const;  				// reserved for something good
	//inline Matrix2D operator+ (const Matrix2D<T> &M) const;
	//inline Matrix2D operator- (const Matrix2D<T> &M) const;
	inline Matrix2D operator* (const T Val) const;
	inline Matrix2D operator/ (const T Val) const;
	//inline Matrix2D operator+ (const T Val) const; 							///< Add the Val to each matrix element, stores result in a new matrix
	//inline Matrix2D operator- (const T Val) const; 							///< Substract the Val from each matrix element, stores result in a new matrix

	inline Matrix2D times (const Matrix2D<T> &M) const; 					///< Arraywise multiplication (A.*B), stores result in a new matrix
	inline Matrix2D divide (const Matrix2D<T> &M) const; 					///< Arraywise division (A./B), stores result in a new matrix

	// Return corresponding index of 1d array
	inline int index1d(int q1, int q2) const;

	T max();
	T maxabs();
	T min();
	Matrix2D<T> abs();
	// It returns maximum between values from class psd2DSlice and argument (VC::zero_f in that case)
	Matrix2D max_of(T val);

	// writeToFileting
	virtual void writeToFile(string filename);
	virtual void writeToFile(string filename, Matrix2D<T> &grid_q1, Matrix2D<T> &grid_q2);
	virtual void readFromFile(string filename, int column = 1);
	virtual void readFromFile(string filename, const Matrix2D<T> grid_q1, const Matrix2D<T> grid_q2);

	// slices - get 1D slice from 2D array
	Matrix1D<T> xSlice(int p_q1) const;
	Matrix1D<T> ySlice(int p_q2) const;

};

/**
 * Matrix 3D class
 *
 * Matrixes and operations.
 */
template <typename T>
class Matrix3D {
private:
	/// Plane array of values. All rows saved in the memory one after anouther as a big array.
	T *plane_array;
	/// Matrix array (array of links to other arrays). Final likns pointed to the memory addresses of the plane array. Matrix[x][y] can be used.
	T ***matrix_array;
public:
	bool initialized;														///< Flag, equal true if initialized
	int size_q1, size_q2, size_q3;												///< size x, size y, size z
	string name;															///< name of the Matrix

	// constructors and destructors
	/// Default constructor. Do nothing.
	Matrix3D() { initialized = false; };
	Matrix3D( const Matrix3D<T> &M );
	Matrix3D( int size_q1, int size_q2, int size_q3 );
	~Matrix3D();

	virtual void AllocateMemory(int size_q1, int size_q2, int size_q3);

	// Operators
	inline T** operator[] (int i); 											///< Return the i-th pointer to 2d-array. Next [j][k] can be applied, so we have regular [i][j][k].
	inline T** operator[] (int i) const { return matrix_array[i]; }
	inline T& operator() (int q1, int q2, int q3); 							///< Return the (x,y,z) value of matrix
	inline T& Value (int q1, int q2, int q3) { return operator()(q1, q2, q3); }	///< Return the (x,y,z) value of matrix
	inline Matrix3D<T>& MatrixArray () { return *this; }					///< Return pointer to the instance of the class.

	// The following operators modify the matrix they applied to
	inline Matrix3D& operator= (const Matrix3D<T> &M);
//	inline Matrix3D& operator= (const Matrix2D<T> &M); 						// confusing function
	inline Matrix3D& operator= (const T Val);

	// unary
	inline const Matrix3D& operator+() const { return *this; }
	inline const Matrix3D operator-() const { return ((*this)*(-1)); }

	//inline Matrix3D& operator*= (const Matrix3D<T> &M); 					// reserved for something good
	//inline Matrix3D& operator/= (const Matrix3D<T> &M); 					// reserved for something good
	inline Matrix3D& operator+= (const Matrix3D<T> &M);
	inline Matrix3D& operator-= (const Matrix3D<T> &M);
	inline Matrix3D& operator*= (const T Val);
	inline Matrix3D& operator/= (const T Val);
	inline Matrix3D& operator+= (const T Val); 								///< Add the Val to each matrix element, stores result in the matrix it's applied to
	inline Matrix3D& operator-= (const T Val); 								///< Substract the Val from each matrix element, stores result in the matrix it's applied to

	inline Matrix3D& times_equal (const Matrix3D<T> &M); 					///< Arraywise multiplication (A.*B), stores result in the matrix it's applied to
	inline Matrix3D& divide_equal (const Matrix3D<T> &M); 					///< Arraywise division (A./B), stores result in the matrix it's applied to

	// The following operators save the result to a new matrix
	//inline Matrix3D operator* (const Matrix3D<T> &M) const; 				// reserved for something good
	//inline Matrix3D operator/ (const Matrix3D<T> &M) const;  				// reserved for something good
	inline Matrix3D operator+ (const Matrix3D<T> &M) const;
	inline Matrix3D operator- (const Matrix3D<T> &M) const;
	inline Matrix3D operator* (const T Val) const;
	inline Matrix3D operator/ (const T Val) const;
	//inline Matrix3D operator+ (const T Val) const; 							///< Add the Val to each matrix element, stores result in a new matrix
	//inline Matrix3D operator- (const T Val) const; 							///< Substract the Val from each matrix element, stores result in a new matrix

	inline Matrix3D times (const Matrix3D<T> &M) const; 					///< Arraywise multiplication (A.*B), stores result in a new matrix
	inline Matrix3D divide (const Matrix3D<T> &M) const; 					///< Arraywise division (A./B), stores result in a new matrix

	// Saving (loading) of a matrix into (from) file
	virtual void writeToFile(string filename, string info = ""); 										///< Save matrix to a file
	virtual void writeToFile(string filename, Matrix3D<T> &grid_q1, Matrix3D<T> &grid_q2, Matrix3D<T> &grid_q3); ///< Save matrix to a file, including grid
	virtual void readFromFile(string filename, int column = 1);  									///< Load matrix from a file
	virtual void readFromFile(string filename, const Matrix3D<T> grid_q1, const Matrix3D<T> grid_q2, const Matrix3D<T> grid_q3); ///< Load matrix to a file

	// Some other stuff
	string change_ind;														///< Variables useful for tracking of changes (time of change can be stored here)

	inline int index1d(int q1, int q2, int q3);								///< Returns index of the element (x,y,z) in 1d array

	T min();
	T max();
	T maxabs();
	Matrix3D<T> abs();

	// slices - get 2D slice from 3D array
	Matrix2D<T> xSlice(int p_q1) const;
	Matrix2D<T> ySlice(int p_q2) const;
	Matrix2D<T> zSlice(int p_q3) const;

	// slices - get 1D slice from 3D array
	Matrix1D<T> xySlice(int p_i1, int p_i2) const;
	Matrix1D<T> yzSlice(int p_i2, int p_i3) const;
	Matrix1D<T> xzSlice(int p_i1, int p_i3) const;


};


/**
 * Matrix 3D class
 *
 * Matrixes and operations.
 */
template <typename T> class Matrix4D {
private:
	/// Plane array of values. All rows saved in the memory one after anouther as a big array.
	T *plane_array;
	/// Matrix array (array of links to other arrays). Final likns pointed to the memory addresses of the plane array. Matrix[x][y] can be used.
	T ****matrix_array;
public:
	bool initialized;														///< Flag, equal true if initialized
	int size_w, size_x, size_y, size_z;												///< size x, size y, size z
	string name;															///< name of the Matrix

	// constructors and destructors
	/// Default constructor. Do nothing.
	Matrix4D() { initialized = false; };
	Matrix4D( const Matrix4D<T> &M );
	Matrix4D( int size_w, int size_x, int size_y, int size_z );
	~Matrix4D();

	virtual void AllocateMemory(int size_w, int size_x, int size_y, int size_z);

	// Operators
	inline T*** operator[] (int i); 											///< Return the i-th pointer to 2d-array. Next [j][k] can be applied, so we have regular [i][j][k].
	inline T*** operator[] (int i) const { return matrix_array[i]; }
	inline T& operator() (int w, int x, int y, int z); 							///< Return the (x,y,z) value of matrix
	inline T& Value (int w, int x, int y, int z) { return operator()(w, x, y, z); }	///< Return the (x,y,z) value of matrix
	inline Matrix4D<T>& MatrixArray () { return *this; }					///< Return pointer to the instance of the class.

	// The following operators modify the matrix they applied to
	inline Matrix4D& operator= (const Matrix4D<T> &M);
	inline Matrix4D& operator= (const T Val);

	// unary
	inline const Matrix4D& operator+() const { return *this; }
	inline const Matrix4D operator-() const { return ((*this)*(-1)); }

	//inline Matrix4D& operator*= (const Matrix4D<T> &M); 					// reserved for something good
	//inline Matrix4D& operator/= (const Matrix4D<T> &M); 					// reserved for something good
	inline Matrix4D& operator+= (const Matrix4D<T> &M);
	inline Matrix4D& operator-= (const Matrix4D<T> &M);
	inline Matrix4D& operator*= (const T Val);
	inline Matrix4D& operator/= (const T Val);
	inline Matrix4D& operator+= (const T Val); 								///< Add the Val to each matrix element, stores result in the matrix it's applied to
	inline Matrix4D& operator-= (const T Val); 								///< Substract the Val from each matrix element, stores result in the matrix it's applied to

	inline Matrix4D& times_equal (const Matrix4D<T> &M); 					///< Arraywise multiplication (A.*B), stores result in the matrix it's applied to
	inline Matrix4D& divide_equal (const Matrix4D<T> &M); 					///< Arraywise division (A./B), stores result in the matrix it's applied to

	// The following operators save the result to a new matrix
	//inline Matrix4D operator* (const Matrix4D<T> &M) const; 				// reserved for something good
	//inline Matrix4D operator/ (const Matrix4D<T> &M) const;  				// reserved for something good
	inline Matrix4D operator+ (const Matrix4D<T> &M) const;
	inline Matrix4D operator- (const Matrix4D<T> &M) const;
	inline Matrix4D operator* (const T Val) const;
	inline Matrix4D operator/ (const T Val) const;
	//inline Matrix4D operator+ (const T Val) const; 							///< Add the Val to each matrix element, stores result in a new matrix
	//inline Matrix4D operator- (const T Val) const; 							///< Substract the Val from each matrix element, stores result in a new matrix

	inline Matrix4D times (const Matrix4D<T> &M) const; 					///< Arraywise multiplication (A.*B), stores result in a new matrix
	inline Matrix4D divide (const Matrix4D<T> &M) const; 					///< Arraywise division (A./B), stores result in a new matrix

	// Saving (loading) of a matrix into (from) file
	virtual void writeToFile(string filename, string info = ""); 										///< Save matrix to a file
	virtual void writeToFile(string filename, Matrix4D<T> &grid_w, Matrix4D<T> &grid_x, Matrix4D<T> &grid_y, Matrix4D<T> &grid_z); ///< Save matrix to a file, including grid
	virtual void readFromFile(string filename, int column = 1);  									///< Load matrix from a file
	virtual void readFromFile(string filename, const Matrix4D<T> grid_w, const Matrix4D<T> grid_x, const Matrix4D<T> grid_y, const Matrix4D<T> grid_z); ///< Load matrix to a file

	// Some other stuff
	string change_ind;														///< Variables useful for tracking of changes (time of change can be stored here)

	inline int index1d(int w, int x, int y, int z);								///< Returns index of the element (x,y,z) in 1d array

	T min();
	T max();
	T maxabs();
	Matrix4D<T> abs();

	// slices - get 3D slice from 4D array
	Matrix3D<T> wSlice(int p_w) const;
	Matrix3D<T> xSlice(int p_x) const;
	Matrix3D<T> ySlice(int p_y) const;
	Matrix3D<T> zSlice(int p_z) const;

	// slices - get 2D slice from 4D array
	Matrix2D<T> wxSlice(int p_w, int p_x) const;
	Matrix2D<T> wySlice(int p_w, int p_y) const;
	Matrix2D<T> wzSlice(int p_w, int p_z) const;
	Matrix2D<T> xySlice(int p_x, int p_y) const;
	Matrix2D<T> xzSlice(int p_x, int p_z) const;
	Matrix2D<T> yzSlice(int p_y, int p_z) const;

	// slices - get 1D slice from 4D array
	Matrix1D<T> wxySlice(int p_w, int p_x, int p_y) const;
	Matrix1D<T> wxzSlice(int p_w, int p_x, int p_z) const;
	Matrix1D<T> wyzSlice(int p_w, int p_y, int p_z) const;
	Matrix1D<T> xyzSlice(int p_x, int p_y, int p_z) const;
};


/** Diagonal matrix.
 * This method of storage for matrices is convenient for diagonal (spread) matrices.
 * Stored as map (diagonal number, 1d diagonal array)
 * The USED diagonals of the matrix are stored in 1d arrays.
 */
typedef map <int , Matrix1D<double> > DiagMatrix;

/** Model matrix (or related matrices)
 * It is based on Diagonal matrix and have methods for conversion from 3D or 2D PSD (and related) arrays into 1d array of unknown elements
 */
class CalculationMatrix : public DiagMatrix {
public:

	bool initialized;
	int x_size, y_size, total_size;
	// flag, if needs to be recalculated
	string change_ind;						///< Variables useful for changes tracking (store here time when changed)

	// Constructors
	CalculationMatrix() { this->initialized = false; }
	// !!! CalculationMatrix(int x_size, int y_size = 0, int z_size = 0, int n_of_diags = 1);
	CalculationMatrix(int x_size, int y_size = 1, int z_size = 1, int n_of_diags = 1);

	// Initialization
	// !!! void Initialize(int x_size, int y_size = 0, int z_size = 0, int n_of_diags = 1);
	void Initialize(int x_size, int y_size = 1, int z_size = 1, int n_of_diags = 1);

	// Returns 1d index for multiple dimension array (2D or 3D)
	int index1d(int x, int y = 0, int z = 0);

	// Save to a file
	void writeToFile(string filename);

	// Operators
	Matrix1D<double> operator* (Matrix1D<double> &V) const;

};

int index1d2(int x, int y = 0, int z = 0);

#endif
