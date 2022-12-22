/** Matrix 1D, 2D, 3D and 4D and operations with them
 *
 * \file Matrix.cpp
 *
 * File has 1D-class, 2D-class 3D-class and 4D-class of matrixes and various functions to work with them.
 *
 * \author Developed under supervision of the PI Yuri Shprits
 *
 * \brief Matrix 1D, 2D, 3D and 4D and operations with them
 *
 */
#ifndef matrix_array_MATRIX_CPP
#define matrix_array_MATRIX_CPP


#if defined(_WIN32) || defined(_WIN64)
	 // Matrix2D<T>::writeToBinaryFile uses fopen function. Which apparently is old and deprecated in CRT (C runtime library)
	 // This directive disables the warning in Visual Studio. (see https://docs.microsoft.com/en-us/cpp/c-runtime-library/security-features-in-the-crt?view=msvc-170)
	#define _CRT_SECURE_NO_WARNINGS
#endif

#include "Matrix.h"
#include "Logger.h"
#include <inttypes.h>
#include <limits>
#include <algorithm>

using namespace std;

#if defined(_WIN32) || defined(_WIN64)
	#define strncasecmp _strnicmp
	#define strcasecmp _stricmp
#endif



const double err = 1e-6;

// #define DEBUG_MODE

// Memory related functions

/// Allocating memory for 1D matrix
template<class T>inline T* matrix(size_t Rows)
{
	T *m=new T[Rows];
	//	assert(m!=NULL);
	if (m == NULL) {
		printf("MEMORY_ERRROR: Memory can't be initialized: %zu size", Rows*sizeof(T));
		exit(EXIT_FAILURE);
	}
	return m;
}

/// Initilizing memory for 2D matrix
template<class T>inline T** matrix(size_t Rows, size_t Columns)
{
	// allocating memory for array of pinters
	T **m=new T*[Rows];
	//	assert(m!=NULL);
	if (m == NULL) {
		printf("MEMORY_ERRROR: Memory can't be initialized: %zu size", Rows * sizeof(T));
		exit(EXIT_FAILURE);
	}
	// allocating memory for data array
	m[0] = new T[Rows * Columns];
	//	assert(m[0]!=NULL);
	if (m[0] == NULL) {
		printf("MEMORY_ERRROR: Memory can't be initialized: %zu size", Rows * Columns * sizeof(T));
		exit(EXIT_FAILURE);
	}
	// assign pointers to data ranges
	for(size_t i=1; i<Rows; i++) m[i] = m[i-1] + Columns;
	return m;
}

/// Initializing memory for 3D matrix
template<class T>inline T*** matrix(size_t size_x, size_t size_y, size_t size_z)
{
	// allocating memory for array of pointers to pointers
	T ***m=new T**[size_x];
	//	assert(m!=NULL);
	if (m == NULL) {
		printf("MEMORY_ERRROR: Memory can't be initialized: %zu size", size_x * sizeof(T));
		exit(EXIT_FAILURE);
	}
	for (size_t x = 0; x < size_x; x++) {
		// for each pointer allocating memory for array of pointers
		m[x] = new T*[size_y];
		//		assert(m[x]!=NULL);
		if (m[x] == NULL) {
			printf("MEMORY_ERRROR: Memory can't be initialized: %zu size", size_y * sizeof(T));
			exit(EXIT_FAILURE);
		}
	}
	//allocating memory for data array
	m[0][0] = new T[size_x * size_y * size_z];
	//	assert(m[0][0]!=NULL);
	if (m[0][0] == NULL) {
		printf("MEMORY_ERRROR: Memory can't be initialized: %zu size", size_x * size_y * size_z * sizeof(T));
		exit(EXIT_FAILURE);
	}
	for (size_t x = 0; x < size_x; x++) {
		for (size_t y = 0; y < size_y; y++) {
			// assign pointers to data ranges
			m[x][y] = m[0][0] + (x*size_y + y)*size_z;
		}
	}

	return m;
}

/// Initializing memory for 4D matrix
template<class T>inline T**** matrix(size_t size_w, size_t size_x, size_t size_y, size_t size_z)
{
	// allocating memory for array of pointers to pointers
	T ****m=new T***[size_w];
	//	assert(m!=NULL);

	if (m == NULL) {
		printf("MEMORY_ERRROR: Memory can't be initialized: %zu size", size_w * sizeof(T));
		exit(EXIT_FAILURE);
	}

	for (size_t w = 0; w < size_w; w++) {
		// for each pointer allocating memory for array of pointers
		m[w] = new T**[size_x];
		//		assert(m[x]!=NULL);
		if (m[w] == NULL) {
			printf("MEMORY_ERRROR: Memory can't be initialized: %zu size", size_x * sizeof(T));
			exit(EXIT_FAILURE);
		}

		for (size_t x = 0; x < size_x; x++) {
			// for each pointer allocating memory for array of pointers
			m[w][x] = new T*[size_y];
			//		assert(m[x]!=NULL);
			if (m[w][x] == NULL) {
				printf("MEMORY_ERRROR: Memory can't be initialized: %zu size", size_y * sizeof(T));
				exit(EXIT_FAILURE);
			}
		}
	}

	//allocating memory for data array
	m[0][0][0] = new T[size_w * size_x * size_y * size_z];
	//	assert(m[0][0]!=NULL);
	if (m[0][0][0] == NULL) {
		printf("MEMORY_ERRROR: Memory can't be initialized: %zu size", size_w * size_x * size_y * size_z * sizeof(T));
		exit(EXIT_FAILURE);
	}
	for (size_t w = 0; w < size_w; w++) {
		for (size_t x = 0; x < size_x; x++) {
			for (size_t y = 0; y < size_y; y++) {
				// assign pointers to data ranges
				m[w][x][y] = m[0][0][0] + ((w*size_x + x)*size_y + y)*size_z;
			}
		}
	}

	return m;
}


/// Freeing memory for 1D matrix
template<class T>inline void free_matrix(T* m) {
	delete[] m;
}

/// Freeing memory for 2D matrix
template<class T>inline void free_matrix(T** m) {
	delete[](m[0]);
	delete[](m);
}

/// Freeing memory for 3D matrix
template<class T>inline void free_matrix(T*** m, size_t size_x, size_t size_y) {
	delete[](m[0][0]);
	for (size_t x = 0; x < size_x; x++) {
		delete[](m[x]);
	}
	delete[](m);
}


/// Freeing memory for 4D matrix
template<class T>inline void free_matrix(T**** m, size_t size_w, size_t size_x, size_t size_y) {
	delete[](m[0][0][0]);
	for (size_t w = 0; w < size_w; w++) {
		for (size_t x = 0; x < size_x; x++) {
			delete[](m[w][x]);
		}
		delete[](m[w]);
	}
	delete[](m);
}


//////////////////////////////////////////////////
//
//	Matrix 1D
//
//////////////////////////////////////////////////

/**
* Constructor.
*
* Runs allocating memory function and store matrix name.
*
* \param size_q1 - size of the matrix
* \param name - name of the matrix
*
*/
template<class T>
Matrix1D<T>::Matrix1D( size_t size_t , const std::string& name) {
	initialized = false;
	// !!! this->name = name;
	// allocating memory
	AllocateMemory(size_t);
}

/**
* Constructor.
*
* Make new matrix equal to Matrix M.
*
* \param &M - matrix M
*
*/
template<class T>
Matrix1D<T>::Matrix1D( const Matrix1D<T> &M ) {
	initialized = false;
	this->operator = (M);
}

/**
* Destructor. Destruct the class.
*/
template<class T>
Matrix1D<T>::~Matrix1D() {
	if (initialized) free_matrix<T>(matrix_array);
}

/**
* Allocating memory
*
* \param size_q1 - size x
*/
template<class T>
void Matrix1D<T>::AllocateMemory(size_t size_q1) {
	this->size_q1 = size_q1;
	num_elements = size_q1;
	// using inline template for memory allocation
	matrix_array = matrix<T>(size_q1);
	plane_array = matrix_array;
	initialized = true;
}


/**
* Operator [i], returns value of element i.
* If  DEBUG_MODE defined, check if matrix has been initialized.
*
* \param i - number of element to return
*/
template<class T>
inline T& Matrix1D<T>::operator[](size_t i1) {
#ifdef DEBUG_MODE
	if (!initialized) {
		printf("MATRIX_ERROR: Using unitialized matrix");
		exit(EXIT_FAILURE);
	}
#endif
	// return i-value
	return matrix_array[i1];
}

/**
* Operator [i], returns value of element i, version returns 'const' value, can not be later modified.
* If  DEBUG_MODE defined, check if matrix has been initialized.
*
* \param i - number of element to return
*/
template<class T>
inline T& Matrix1D<T>::operator[](size_t i1) const {
#ifdef DEBUG_MODE
	if (!initialized) {
		printf("MATRIX_ERROR: Using unitialized matrix");
		exit(EXIT_FAILURE);
	}
#endif
	// return i-value
	return matrix_array[i1];
}

/**
* Operator (x), returns value of element x.
* If  DEBUG_MODE defined, check if matrix has been initialized.
* No dofference between [] and () operators for 1d-matrix class.
*
* \param x - number of element to return
*/
template<class T>
inline T& Matrix1D<T>::operator()(int i1) {
#ifdef DEBUG_MODE
	if (!initialized) {
		printf("MATRIX_ERROR: Using unitialized matrix");
		exit(EXIT_FAILURE);
	}
#endif
	// return x-value
	return matrix_array[i1];
}

/**
* Make matrix equal to matrix M
* Return the same instance of class Matrix.
*
* \param &M - matrix M
*/
template<class T>
inline Matrix1D<T>& Matrix1D<T>::operator= (const Matrix1D<T> &M) {
	// Check, if both matrix aren't the same one
	if (this!=&M) {
		// check, if RHS matrix was initialized
		if (M.initialized) {
			// free LHS matrix if it was initialized
			if (initialized && (size_q1 != M.size_q1)) {
				free_matrix<T>(matrix_array);
				initialized = false;
			}
			if (!initialized) {
				this->size_q1 = M.size_q1;
				// allocating memory for LHS matrix (and creating the correct pointers M[0] - M[N] to the matrix rows)
				this->AllocateMemory(this->size_q1);
			}

			// XXX: This might cause problems with openmp!
			//this->name = M.name;
			// this->name = M.name.c_str();

			// fast values copying as a memory range
			memcpy( matrix_array, M.matrix_array, this->size_q1 * sizeof( T ) );
		} else {
			this->initialized = false;
		}
	}
	return *this;
}

/**
* Matrix summation, result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix1D<T>& Matrix1D<T>::operator+= (const Matrix1D<T> &M) {
	for (size_t i1 = 0; i1 < size_q1; i1++)
		matrix_array[i1] += M.matrix_array[i1];
	return *this;
}

/**
* Matrix subtraction, result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix1D<T>& Matrix1D<T>::operator-= (const Matrix1D<T> &M) {
	for (size_t i1 = 0; i1 < size_q1; i1++)
		matrix_array[i1] -= M.matrix_array[i1];
	return *this;
}

/**
* Multiplication to a value. Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix1D<T>& Matrix1D<T>::operator*= (const T Val) {
	for (size_t i1 = 0; i1 < size_q1; i1++)
		matrix_array[i1] *= Val;
	return *this;
}

/**
* Division by a value. Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix1D<T>& Matrix1D<T>::operator/= (const T Val) {
	for (size_t i1 = 0; i1 < size_q1; i1++)
		matrix_array[i1] /= Val;
	return *this;
}


/**
* Make matrix equal to value Val.
* Return the same instance of class Matrix.
*
* \param Val - value val
*/
template<class T>
inline Matrix1D<T>& Matrix1D<T>::operator= (const T Val) {
	for (size_t i1 = 0; i1 < size_q1; i1++)
		matrix_array[i1] = Val;
	return *this;
}


/**
* Multiply a matrix to a value Val.
* Return new instance of class Matrix.
*
* \param Val - value Val
*/
template<class T>
inline Matrix1D<T> Matrix1D<T>::operator* (const T Val) const {
	Matrix1D<T> Tmp(*this);
	for (size_t i1 = 0; i1 < size_q1; i1++)
		Tmp[i1] = matrix_array[i1] * Val;
	return Tmp;
}

/**
* Divide a matrix to a value Val.
* Return new instance of class Matrix.#define DEBUG_MODE
*
* \param Val - value Val
*/
template<class T>
inline Matrix1D<T> Matrix1D<T>::operator/ (const T Val) const {
	Matrix1D<T> Tmp(*this);
	for (size_t i1 = 0; i1 < size_q1; i1++)
		Tmp[i1] = matrix_array[i1] / Val;
	return Tmp;
}

/**
* Multiply all values of matrix to values of matrix M.
*
* \param &M - matrix M.
*/
template<class T>
inline Matrix1D<T> Matrix1D<T>::times (const Matrix1D<T> &M) const {
	Matrix1D<T> Tmp(size_q1);
	for (size_t i1 = 0; i1 < size_q1; i1++)
		Tmp[i1] = matrix_array[i1] * M.matrix_array[i1];
	return Tmp;
}

/**
* Divide all values of matrix to values of matrix M.
*
* \param &M - matrix M.
*/
template<class T>
inline Matrix1D<T> Matrix1D<T>::divide (const Matrix1D<T> &M) const {
	Matrix1D<T> Tmp(size_q1);
	for (size_t i1 = 0; i1 < size_q1; i1++)
		Tmp[i1] = matrix_array[i1] / M.matrix_array[i1];
	return Tmp;
}

/**
* Norm of vector
*/
template<class T>
inline T Matrix1D<T>::norm() const {
	T res = 0;
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		res += matrix_array[i1] * matrix_array[i1];
	}
	return sqrt(res);
}

/*
* Dot product
*/
template<class T>
inline T Matrix1D<T>::dot( const Matrix1D<T> &W ) const {
	if (this->size_q1 != W.size_q1) {
		printf("DOT_PRODUCT: Size is different");
		exit(EXIT_FAILURE);
	}
	T res = 0;
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		res += matrix_array[i1] * W[i1];
	}
	return res;
}


/**
* Write matrix data to file.
*/
template<class T>
void Matrix1D<T>::writeToFile(const std::string& filename) const {
	ofstream output(filename.c_str());
	//if (output==NULL && (filename.find("Debug") == std::string::npos)) {
	if (!output.is_open() && (filename.find("Debug") == std::string::npos)) {
		printf("FILE: Unable to output file: %s", filename.c_str());
		exit(EXIT_FAILURE);
	}
	output << "VARIABLES = \"" << ((this->name!="")?this->name:"function") << "\" "<< endl;
	output << "ZONE T=\"" << filename << "\", I=" << size_q1 << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		output << matrix_array[i1] << endl;
	}
	output.close();
}

/**
* Write matrix data to file with grid.
*/
template<class T>
void Matrix1D<T>::writeToFile(const std::string& filename, const Matrix1D<T> &grid_q1) const {
	ofstream output(filename.c_str());
	//if (output==NULL && (filename.find("Debug") == std::string::npos)) {
	if (!output.is_open() && (filename.find("Debug") == std::string::npos)) {
		printf("FILE: Unable to output file: %s", filename.c_str());
		exit(EXIT_FAILURE);
	}
	output << "VARIABLES = \"" << ((grid_q1.name!="")?grid_q1.name:"x") << "\", \""  << ((this->name!="")?this->name:"function") << "\" "<< endl;
	output << "ZONE T=\"" << filename << "\", I=" << size_q1 << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		output << "\t" << grid_q1[i1] << "\t" << "\t" << matrix_array[i1] << endl;
	}
	output.close();
}


/**
* Read matrix data from file.
*/
template<class T>
void Matrix1D<T>::readFromFile(const std::string& filename) {
	std::string inBuf;
	if (!initialized) {
		printf("MATRIX_ERROR: Using unitialized matrix");
		exit(EXIT_FAILURE);
	} else {
		this->name = filename;
		ifstream input(filename.c_str());
		//if (input != NULL && !input.eof()) {
		if (input.is_open() && !input.eof()) {
			// Skipping first two lines.
			input >> inBuf;
			// !!!			while (strcmp(strupr((char *) inBuf.c_str()), "ZONE") != 0 && !input.eof() ) {
			while (strcasecmp(inBuf.c_str(), "ZONE") != 0 && !input.eof() ) {
				input >> inBuf;
			}
			// read to the end of the line with 'zone'
			input.ignore(9999, '\n');
			for (size_t i1 = 0; i1 < size_q1; i1++) {
				input >> matrix_array[i1];
			}
		} else {
			printf("MATRIX_LOAD_ERROR: Error reading file %s.\n", filename.c_str());
			exit(EXIT_FAILURE);
		}
		input.close();
	}
}



/**
* Function for reading from matlab file in 1-dimension
* Will check the variables, order them in (P, R/L, V, K, Val) format and then set matrix_array to be the variable with the corresponding column number
* This is the same as the readFromFile() function although only compatible with .mat files instead of .plt or other text files
*/
template<class T>
void Matrix1D<T>::readFromMatlabFile(const std::string& file , int columnNumber)
{
#if (MATLAB_CAPABLE)

	Logger::message << "Reading " << file << ": " << endl;

	MATFile *mfPtr; /* MAT-file pointer */
    mxArray *aPtr;  /* mxArray pointer */
    double *realPtr = nullptr; /* pointer to data */
	double *PtrW = nullptr; /* pointer to data */
	double *PtrX = nullptr; /* pointer to data */
    double *PtrY = nullptr; /* pointer to data */
    double *PtrZ = nullptr; /* pointer to data */
	double *PtrFinal = nullptr; /* pointer to data */
	double *PtrReturn = nullptr; /* pointer to data */
	double *PtrL = nullptr; /* pointer to data */
	std::string arr; /*name of variable*/
	std::string field = "arr"; // name of field
	mwSize nElements;       /* number of elements in array */
    mwIndex eIdx;           /* element index */
    const mxArray *fPtr;    /* field pointer */
	const mxArray *fPtrW;    /* field pointer */
	const mxArray *fPtrX;    /* field pointer */
	const mxArray *fPtrY;    /* field pointer */
	const mxArray *fPtrZ;    /* field pointer */
    size_t x; 			/* for index*/
	const char* name;		/* for getting variable names */
	char nameTemp;		/* for getting variable names */
	bool wReached = false;
	bool xReached = false;
	bool yReached = false;
	bool zReached = false;
	bool defaultReached = false;

	mfPtr = matOpen(file.c_str(), "r");
   	if (mfPtr == NULL) {
       	printf("Error opening file %s\n", file.c_str());
	}

	for (int i= 0; i < 3; i++)
	{
		aPtr = matGetNextVariableInfo(mfPtr, &name);
		std::string temp = name;
		if (temp.length() == 1)
		{
			nameTemp = temp[0];
		}
		else
		{
			if (temp.substr(1, 1) == "_")
			{
				nameTemp = temp[0];
			}
			else if ((temp.substr(0, 2) == "pc"))
			{
				nameTemp = 'V';
			}
			else if ((temp.substr(0, 2) == "al"))
			{
				nameTemp = 'K';
			}
			else
			{
				nameTemp = '!';
			}
		}
		switch (nameTemp)
		{
			case 'P':
				wReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrW = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrW != NULL) && (mxGetClassID(fPtrW) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrW)))
       		 				{
								PtrW = mxGetPr(fPtrW);
								//printf("%e \n", PtrW[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'R':
				xReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrX = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrX != NULL) && (mxGetClassID(fPtrX) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrX)))
       		 				{
								PtrX = mxGetPr(fPtrX);
								//printf("%e \n", PtrX[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'L':
				xReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrX = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrX != NULL) && (mxGetClassID(fPtrX) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrX)))
       		 				{
								PtrX = mxGetPr(fPtrX);
								PtrL = mxGetPr(fPtrX);
								//printf("%e \n", PtrX[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'V':
				yReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrY = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrY != NULL) && (mxGetClassID(fPtrY) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrY)))
       		 				{
								PtrY = mxGetPr(fPtrY);
								//printf("%e \n", PtrY[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'K':
				zReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrZ = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrZ != NULL) && (mxGetClassID(fPtrZ) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrZ)))
       		 				{
								PtrZ = mxGetPr(fPtrZ);
								//printf("%e \n", PtrZ[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			default:
				defaultReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtr = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtr != NULL) && (mxGetClassID(fPtr) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtr)))
       		 				{
								PtrFinal = mxGetPr(fPtr);
								//printf("%e \n", PtrFinal[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
		}
	}


	// To get the column number of the return
	if (wReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrW;
		}
	}
	if (xReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrX;
		}
	}
	if (yReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrY;
		}
	}
	if (zReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrZ;
		}
	}
	if (defaultReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrFinal;
		}
	}
	if (columnNumber > 0)
	{
		printf("column number too high \n");
		EXIT_FAILURE;
	}


	//sets the matrix array to be equal to an array of doubles
	for (x = 0; x < size_q1; x++) {
		matrix_array[x] = PtrReturn[x];
	}

    if (matClose(mfPtr) != 0) {
        printf("Error closing file %s\n", file.c_str());
    }

	mxDestroyArray(aPtr);

#endif

}





/**
* Read matrix data from file with grid,
* Checks if the matrix data in the file is the same as the grids that were sent in with error < 1e-8,
* if not within error range will signal error and exit
*
* Overloaded readFromFile function
* \param filename - file to read grids from
* \param grid_q1 - checks grid data against the file data
*/
template<class T>
void Matrix1D<T>::readFromFile(const std::string& filename, const Matrix1D<T>& grid_q1) {
	std::string inBuf;
	double loaded_q1;
	if (!initialized) {
		printf("MATRIX_ERROR: Using unitialized matrix");
		exit(EXIT_FAILURE);
	} else {
		this->name = filename;
		ifstream input(filename.c_str());
		//if (input != NULL && !input.eof()) {
		if (input.is_open() && !input.eof()) {
			// Skipping first two lines. Should search for 'ZONE' and read from following line better - done.
			input >> inBuf;
			// !!!			while (strcmp(strupr((char *) inBuf.c_str()), "ZONE") != 0 && !input.eof() ) {
			while (strcasecmp(inBuf.c_str(), "ZONE") != 0 && !input.eof() ) {
				input >> inBuf;
			}
			// read to the end of the line with 'zone'
			input.ignore(9999, '\n');
			for (size_t i1 = 0; i1 < size_q1; i1++) {
				input >> loaded_q1;
				// check if grid is the same
				if (fabs(log10(loaded_q1) - log10(grid_q1[i1])) > err) {
					printf("MATRIX_LOAD_GRID_ERR: Loading %s: grid mismatch.\nLoaded: %e\nGrid: %e\n", filename.c_str(), loaded_q1, grid_q1[i1]);
					exit(EXIT_FAILURE);
				} else {
					input >> matrix_array[i1];
				}
			}
		} else {
			printf("MATRIX_LOAD_ERROR: Error reading file %s.\n", filename.c_str());
			exit(EXIT_FAILURE);
		}
		input.close();
	}
}








/**
* Function for reading from matlab file in 1-dimension
* Will check the variables in the order they are saved in matlab, thus (P R V K Var) should be the standard
* The variables will be checked against the input grid parameters in order to make sure the right variables/values are being loaded
* This is the same as the readFromFile() function although only compatible with .mat files instead of .plt or other text files
*/
template<class T>
void Matrix1D<T>::readFromMatlabFile(const std::string& file , const Matrix1D<T>& grid_x)
{

#if (MATLAB_CAPABLE)

	Logger::message << "Reading " << file << ": " << endl;

	MATFile *mfPtr; /* MAT-file pointer */
    mxArray *aPtr;  /* mxArray pointer */
	double *PtrW = nullptr; /* pointer to data */
	double *PtrX = nullptr; /* pointer to data */
	double *PtrY = nullptr; /* pointer to data */
	double *PtrZ = nullptr; /* pointer to data */
	double *PtrFinal = nullptr; /* pointer to data */
	double *PtrL = nullptr; /* pointer to data */
	std::string arr; /*name of variable*/
	std::string field = "arr"; // name of field
	mwSize nElements;       /* number of elements in array */
    mwIndex eIdx;           /* element index */
    const mxArray *fPtr;    /* field pointer */
	const mxArray *fPtrW;    /* field pointer */
	const mxArray *fPtrX;    /* field pointer */
	const mxArray *fPtrY;    /* field pointer */
	const mxArray *fPtrZ;    /* field pointer */
	size_t x;
	const char* name;		/* for getting variable names */
	char nameTemp;		/* for getting variable names */
	bool wReached = false;
	bool xReached = false;
	bool yReached = false;
	bool zReached = false;
	bool defaultReached = false;

	mfPtr = matOpen(file.c_str(), "r");
   	if (mfPtr == NULL) {
       	printf("Error opening file %s\n", file.c_str());
	}

	for (int i= 0; i < 3; i++)
	{

		aPtr = matGetNextVariableInfo(mfPtr, &name);
		std::string temp = name;
		if (temp.length() == 1)
		{
			nameTemp = temp[0];
		}
		else
		{
			if (temp.substr(1, 1) == "_")
			{
				nameTemp = temp[0];
			}
			else if ((temp.substr(0, 2) == "pc"))
			{
				nameTemp = 'V';
			}
			else if ((temp.substr(0, 2) == "al"))
			{
				nameTemp = 'K';
			}
			else
			{
				nameTemp = '!';
			}
		}
		switch (nameTemp)
		{
			case 'P':
				wReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrW = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrW != NULL) && (mxGetClassID(fPtrW) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrW)))
       		 				{
								PtrW = mxGetPr(fPtrW);
								//printf("%e \n", PtrW[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'R':
				xReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrX = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrX != NULL) && (mxGetClassID(fPtrX) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrX)))
       		 				{
								PtrX = mxGetPr(fPtrX);
								//printf("%e \n", PtrX[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'L':
				xReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrX = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrX != NULL) && (mxGetClassID(fPtrX) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrX)))
       		 				{
								PtrX = mxGetPr(fPtrX);
								PtrL = mxGetPr(fPtrX);
								//printf("%e \n", PtrX[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'V':
				yReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrY = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrY != NULL) && (mxGetClassID(fPtrY) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrY)))
       		 				{
								PtrY = mxGetPr(fPtrY);
								//printf("%e \n", PtrY[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'K':
				zReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrZ = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrZ != NULL) && (mxGetClassID(fPtrZ) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrZ)))
       		 				{
								PtrZ = mxGetPr(fPtrZ);
								//printf("%e \n", PtrZ[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			default:
				defaultReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtr = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtr != NULL) && (mxGetClassID(fPtr) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtr)))
       		 				{
								PtrFinal = mxGetPr(fPtr);
								//printf("%e \n", PtrFinal[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
		}
	}

	if (wReached)
	{
		PtrX = PtrW;
	}
	if (yReached)
	{
		PtrX = PtrY;
	}
	if (zReached)
	{
		PtrX = PtrZ;
	}





	//sets the matrix array to be equal to an array of doubles
	for (x = 0; x < size_q1; x++) {
		if (fabs(log10(PtrX[x]) - log10(grid_x[x])) > err ) {
		 			printf("MATRIX_LOAD_GRID_ERR: Loading %s: grid mismatch [%zu].\nLoaded: %e\nGrid: %e\n", file.c_str(), x, PtrX[x], grid_x[x]);
					//printf("grid error \n");
					exit(EXIT_FAILURE);
		}
		matrix_array[x] = PtrFinal[x];
	}

    if (matClose(mfPtr) != 0) {
        printf("Error closing file %s\n", file.c_str());
    }

	mxDestroyArray(aPtr);

#endif
}






// // ADDED
// /* Analyze field FNAME in struct array SPTR. */
// static void
// analyzestructarray(const mxArray *sPtr, const char *fName)
// {
//     mwSize nElements;       /* number of elements in array */
//     mwIndex eIdx;           /* element index */
//     mwIndex fPTRIdx;           /* fptr index */
//     const mxArray *fPtr;    /* field pointer */
//     double *realPtr;        /* pointer to data */
//     mwSize nElementsInRealData; /* number of elements in array */


// 	// Goes through all of the structs
// 	// I think all of them only have one struct consisting of a couple fields (one or more of which is the double array which we want)
//     nElements = (mwSize)mxGetNumberOfElements(sPtr);
//     for (eIdx = 0; eIdx < nElements; eIdx++) {
//         fPtr = mxGetField(sPtr, eIdx, fName);
//         nElementsInRealData = (mwSize)mxGetNumberOfElements(fPtr);
//         printf("number of elements in %s: %.2d\n", fName, nElementsInRealData);
//         if ((fPtr != NULL)
//             && (mxGetClassID(fPtr) == mxDOUBLE_CLASS)
//             && (!mxIsComplex(fPtr)))
//         {
//             realPtr = mxGetPr(fPtr);
//             for (fPTRIdx = 0; fPTRIdx < nElementsInRealData; fPTRIdx++)
//             {
//                 printf("%.2f \n", realPtr[fPTRIdx]);
//             }
//         }
//     }
// }

// /* Find struct array ARR in MAT-file FILE.
//  * Pass field name FIELD to analyzestructarray function. */
// /*int findstructure(
//         const char *file,
//         const char *arr,
//         const char *field) {



// }*/















/**
* Return minimum value of the 1d matrix.
* Default value set to 1e99
*/
template<class T>
T Matrix1D<T>::min() const {
	T tmp = 1e99;
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		tmp = (tmp<matrix_array[i1])?tmp:matrix_array[i1];
	}
	return tmp;
}

/**
* Return maximum value of the 1d matrix.
* Default value seet to 0
*/
template<class T>
T Matrix1D<T>::max() const {
	T tmp = 0;
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		tmp = (tmp>matrix_array[i1])?tmp:matrix_array[i1];
	}
	return tmp;
}

/**
* Return absolute maximum value of the 1d matrix.
* Default value set to 0
*/
template<class T>
T Matrix1D<T>::maxabs() const {
	T tmp = 0;
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		tmp = (tmp>fabs((double)matrix_array[i1]))?tmp:fabs((double)matrix_array[i1]);
	}
	return tmp;
}

/**
* Return absolute value of the 1d matrix.
* Changes every element to a positive value with the same magnitude
*/
template<class T>
Matrix1D<T> Matrix1D<T>::abs() const {
	Matrix1D<T> tmp(this->size_q1);
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		tmp[i1] = (matrix_array[i1]>0)?matrix_array[i1]:-matrix_array[i1];
	}
	return tmp;
}

/**
* Replace all matrix malues smaller than val with val.
*
* \param val - new minimal value.
*/
template<class T>
Matrix1D<T>& Matrix1D<T>::max_of(T val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] = std::max(plane_array[i], val);
	}
	return *this;
}

//////////////////////////////////////////////////
//
//	Matrix 2D
//
//////////////////////////////////////////////////

/**
* Constructor.
* Allocate memory.
*
* \param size_q1 - x size
* \param size_q2 - y size
*/
template<class T>
Matrix2D<T>::Matrix2D(size_t size_q1, size_t size_q2) {
	initialized = false;
	// allocating memory
	AllocateMemory(size_q1, size_q2);
}

/**
* Constructor.
* Create new matrix from the Matrix M.
*
* \param &M - Matrix M
*/
template<class T>
Matrix2D<T>::Matrix2D( const Matrix2D<T> &M ) {
	initialized = false;
	this->operator = (M);
}

/**
* Destructor.
*/
template<class T>
Matrix2D<T>::~Matrix2D() {
	if (initialized) free_matrix<T>(matrix_array);
}

/**
* Allocate memory.
*
* \param size_q1 - x size
* \param size_q2 - y size
*/
template<class T>
void Matrix2D<T>::AllocateMemory(size_t size_q1, size_t size_q2) {
	this->size_q1 = size_q1;
	this->size_q2 = size_q2;
	num_elements = size_q1 * size_q2;
	// using matrix inline template to allocate memory
	matrix_array = matrix<T>(size_q1, size_q2);
	plane_array = matrix_array[0];
	initialized = true;
}

/**
* Make matrix equal to Matrix M.
*
* \param &M - Matrix M.
*/
template<class T>
Matrix2D<T>& Matrix2D<T>::operator= (const Matrix2D<T> &M) {
	// check, if both matrix are the same one
	if (this!=&M) {
		// check, if RHS matrix in initialized
		if (M.initialized) {
			// free LHS matrix if it was initialized
			if (initialized && (size_q1 != M.size_q1 || size_q2 != M.size_q2)) {
				free_matrix<T>(matrix_array);
				initialized = false;
			}
			if (!initialized) {
				this->size_q1 = M.size_q1;
				this->size_q2 = M.size_q2;
				// allocating memory for LHS matrix (and creating the correct pointers M[0] - M[N] to the matrix rows)
				this->AllocateMemory(this->size_q1, this->size_q2);
			}
			// XXX !!! this->name = M.name;

			// fast values copying as a memory range
			memcpy(matrix_array[0], M.matrix_array[0], this->size_q1 * this->size_q2 * sizeof(T) );
		} else {
			this->initialized = false;
		}
	}
	return *this;
}

/**
* Make matrix equal to value Val.
*
* \param val - value Val.
*/
template<class T>
inline Matrix2D<T>& Matrix2D<T>::operator= (const T val) {
	if (initialized) {
		for(size_t i = 0; i < num_elements; i++)
		{
			plane_array[i] = val;
		}
	} else {
		printf("MATRIX_ERROR: Using un-itialized matrix");
		exit(EXIT_FAILURE);
	}

	return *this;
}


/**
* Matrix summation, result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix2D<T>& Matrix2D<T>::operator+= (const Matrix2D<T> &M) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] += M.plane_array[i];
	}
	return *this;
}

/**
* Matrix subtraction, result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix2D<T>& Matrix2D<T>::operator-= (const Matrix2D<T> &M) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] -= M.plane_array[i];
	}
	return *this;
}


/**
* Multiplication to a value. Result is stored into the left hand side matrix
*/
template<class T>
inline Matrix2D<T>& Matrix2D<T>::operator*= (const T Val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] *= Val;
	}
	return *this;
}


/**
* Multiply matrix to Val.
*
* \param Val - value Val.
*/
template<class T>
inline Matrix2D<T> Matrix2D<T>::operator* (const T Val) const {
	Matrix2D<T> Tmp(size_q1, size_q2);
	for (size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] * Val;
	}
	return Tmp;
}

/**
* Divide matrix to Val.
*
* \param Val - value Val.
*/
template<class T>
inline Matrix2D<T> Matrix2D<T>::operator/(const T Val) const {
	Matrix2D<T> Tmp(size_q1, size_q2);
	for (size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] / Val;
	}
	return Tmp;
}

/**
* Divide all values of matrix to values of matrix M.
*
* \param &M - matrix M.
*/
template<class T>
inline Matrix2D<T> Matrix2D<T>::divide (const Matrix2D<T> &M) const {
	Matrix2D<T> Tmp(size_q1, size_q2);
	for (size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] / M.plane_array[i];
	}
	return Tmp;
}

/**
* Multiply all values of matrix to values of matrix M.
*
* \param &M - matrix M.
*/
template<class T>
inline Matrix2D<T> Matrix2D<T>::times (const Matrix2D<T> &M) const {
	Matrix2D<T> Tmp(size_q1, size_q2);
	for (size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] * M.plane_array[i];
	}
	return Tmp;
}


/**
* Return maximum value of the 2d matrix.
* Default value seet to 0
*/
template<class T>
T Matrix2D<T>::max() const {
	T tmp = -std::numeric_limits<T>::infinity();
	for(size_t i = 0; i < num_elements; i++)
	{
		tmp = std::max(tmp, plane_array[i]);
	}
	return tmp;
}

/**
* Return absolute maximum value of the 2d matrix.
* Default value set to 0
*/
template<class T>
T Matrix2D<T>::maxabs() const {
	T tmp = 0;
	for(size_t i = 0; i < num_elements; i++)
	{
		tmp = std::max(tmp, std::abs(plane_array[i]));
	}
	return tmp;
}

/**
* Return exp(matrix*multiplicator) of the 2d matrix.
*/
template<class T>
Matrix2D<T> Matrix2D<T>::exp(double multiplicator) const {
	Matrix2D<T> tmp(size_q1, size_q2);
	for(size_t i = 0; i < num_elements; i++)
	{
		tmp.plane_array[i] = std::exp(plane_array[i]*multiplicator);
	}
	return tmp;
}

/**
* Return minimum value of the 2d matrix.
* Default value set to infinity
*/
template<class T>
T Matrix2D<T>::min() const {
	T tmp = std::numeric_limits<T>::infinity();
	for(size_t i = 0; i < num_elements; i++)
	{
		tmp = std::min(tmp, plane_array[i]);
	}
	return tmp;
}

/**
* Return absolute value of the 2d matrix.
* Changes every element to a positive value with the same magnitude
*/
template<class T>
Matrix2D<T> Matrix2D<T>::abs() const {
	Matrix2D<T> tmp(size_q1, size_q2);
	for(size_t i = 0; i < num_elements; i++)
	{
		tmp.plane_array[i] = std::abs(plane_array[i]);
	}
	return tmp;
}


/**
* Replace all matrix malues smaller than val with val.
*
* \param val - new minimal value.
*/
template<class T>
Matrix2D<T>& Matrix2D<T>::max_of(T val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] = std::max(plane_array[i], val);
	}
	return *this;
}

/**
* Returns corresponding index of 2d matrix if represented as a 1d array
* \param x,y - index of element in every dimension for the 2d matrix
*/
template<class T>
inline size_t Matrix2D<T>::index1d(size_t x, size_t y) const {
	return x*size_q2 + y;
}



/**
* Writes the matrix to a file.
* File has two header lines.
*
* \param filename - file name
*/
template<class T>
void Matrix2D<T>::writeToFile(const std::string& filename) const {
	ofstream output(filename.c_str());
	//if (output==NULL && (filename.find("Debug") == std::string::npos)) {
	if (!output.is_open() && (filename.find("Debug") == std::string::npos)) {
		printf("FILE: Unable to output file: %s", filename.c_str());
		exit(EXIT_FAILURE);
	}
	output << "VARIABLES = \""<< ((this->name!="")?this->name:"f") << "\" "<< endl;
	output << "ZONE T=\"" << filename << "\", I=" << size_q2 << ", J= " << size_q1 << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		for (size_t i2 = 0; i2 < size_q2; i2++) {
			output << matrix_array[i1][i2] << endl;
		}
	}
	output.close();
}


/**
* Write the matrix to a file using other two matrixes as a grid.
* Simply that means - write all three matrixes to a file.
*/
template<class T>
void Matrix2D<T>::writeToFile(const std::string& filename, const Matrix2D<T> &grid_x, const Matrix2D<T> &grid_y) const {
	ofstream output(filename.c_str());
	output << "VARIABLES = \"" << ((grid_x.name!="")?grid_x.name:"x") << "\", \"" << ((grid_y.name!="")?grid_y.name:"y") << "\", \"" << ((this->name!="")?this->name:"f") << "\" "<< endl;
	output << "ZONE T=\"" << filename << "\", I=" << size_q2 << ", J=" << size_q1 << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		for (size_t i2 = 0; i2 < size_q2; i2++) {
			output << "\t" << grid_x[i1][i2] << "\t" << grid_y[i1][i2] << "\t" << matrix_array[i1][i2] << endl;
		}
	}
	output.close();
}

template<class T>
void Matrix2D<T>::writeToBinaryFile(const std::string& filename) const {
	std::ofstream outputFile(filename, std::ios::binary | std::ios::out);
    if (outputFile.is_open()) {
		int32_t size_array[2] =  { (int32_t)size_q1, (int32_t)size_q2 };
		outputFile.write((char*)size_array, 2 * sizeof(int32_t));
		if (!outputFile.good()){
			printf("Writing error");
			exit(EXIT_FAILURE);
		}

		outputFile.write((char*)plane_array, num_elements * sizeof(T));
		if (!outputFile.good()){
			printf("Writing error");
			exit(EXIT_FAILURE);
		}
    } else {
        printf("MATRIX_WRITE_ERROR: Error writing file %s.\n", filename.c_str());
        exit(EXIT_FAILURE);
    }
    outputFile.close();
}

/**
* Read matrix data from file with grid, by column
*
* Overloaded readFromFile function
* \param filename - file to read grids from
* \param read_column - read up to this column from file
*/
template<class T>
void Matrix2D<T>::readFromFile(const std::string& filename, int read_column) {
	std::string inBuf;
	if (!initialized) {
		printf("MATRIX_ERROR: Using unitialized matrix");
		exit(EXIT_FAILURE);
	} else {
		this->name = filename;
		ifstream input(filename.c_str());
		//if (input != NULL && !input.eof()) {
		if (input.is_open() && !input.eof()) {
			// Skipping first two lines.
			input >> inBuf;
			// !!!			while (strcmp(strupr((char *) inBuf.c_str()), "ZONE") != 0 && !input.eof() ) {
			while (strcasecmp(inBuf.c_str(), "ZONE") != 0 && !input.eof() ) {
				input >> inBuf;
			}
			// read to the end of the line with 'zone'
			input.ignore(9999, '\n');
			for (size_t i1 = 0; i1 < size_q1; i1++) {
				for (size_t i2 = 0; i2 < size_q2; i2++) {
					for (int column = 1; column < read_column; column++) input >> inBuf; // skip first columns
					input >> matrix_array[i1][i2];
				}
			}
		} else {
			printf("MATRIX_LOAD_ERROR: Error reading file %s.\n", filename.c_str());
			exit(EXIT_FAILURE);
		}
		input.close();
	}
}



/**
* Function for reading from matlab file in 2-dimensions
* Will check the variables, order them in (P, R/L, V, K, Val) format and then set matrix_array to be the variable with the corresponding column number
* This is the same as the readFromFile() function although only compatible with .mat files instead of .plt or other text files
*/
template<class T>
void Matrix2D<T>::readFromMatlabFile(const std::string& file ,  int columnNumber)
{

#if (MATLAB_CAPABLE)

	Logger::message << "Reading " << file << ": " << endl;

	MATFile *mfPtr; /* MAT-file pointer */
    mxArray *aPtr;  /* mxArray pointer */
	double *PtrW = nullptr; /* pointer to data */
	double *PtrX = nullptr; /* pointer to data */
    double *PtrY = nullptr; /* pointer to data */
    double *PtrZ = nullptr; /* pointer to data */
	double *PtrFinal = nullptr; /* pointer to data */
	double *PtrReturn = nullptr; /* pointer to data */
	double *PtrL = nullptr; /* pointer to data */
	std::string arr; /*name of variable*/
	std::string field = "arr"; // name of field
	mwSize nElements;       /* number of elements in array */
    mwIndex eIdx;           /* element index */
    const mxArray *fPtr;    /* field pointer */
	const mxArray *fPtrW;    /* field pointer */
	const mxArray *fPtrX;    /* field pointer */
	const mxArray *fPtrY;    /* field pointer */
	const mxArray *fPtrZ;    /* field pointer */
	size_t x, y;
	const char* name;		/* for getting variable names */
	char nameTemp;		/* for getting variable names */
	
	bool wReached = false;
	bool xReached = false;
	bool yReached = false;
	bool zReached = false;
	bool defaultReached = false;

	mfPtr = matOpen(file.c_str(), "r");
   	if (mfPtr == NULL) {
       	printf("Error opening file %s\n", file.c_str());
	}

	for (int i= 0; i < 3; i++)
	{
		aPtr = matGetNextVariableInfo(mfPtr, &name);
		std::string temp = name;
		if (temp.length() == 1)
		{
			nameTemp = temp[0];
		}
		else
		{
			if (temp.substr(1, 1) == "_")
			{
				nameTemp = temp[0];
			}
			else if ((temp.substr(0, 2) == "pc"))
			{
				nameTemp = 'V';
			}
			else if ((temp.substr(0, 2) == "al"))
			{
				nameTemp = 'K';
			}
			else
			{
				nameTemp = '!';
			}
		}
		switch (nameTemp)
		{
			case 'P':
				wReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrW = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrW != NULL) && (mxGetClassID(fPtrW) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrW)))
       		 				{
								PtrW = mxGetPr(fPtrW);
								//printf("%e \n", PtrW[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'R':
				xReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrX = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrX != NULL) && (mxGetClassID(fPtrX) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrX)))
       		 				{
								PtrX = mxGetPr(fPtrX);
								//printf("%e \n", PtrX[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'L':
				xReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrX = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrX != NULL) && (mxGetClassID(fPtrX) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrX)))
       		 				{
								PtrX = mxGetPr(fPtrX);
								PtrL = mxGetPr(fPtrX);
								//printf("%e \n", PtrX[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'V':
				yReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrY = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrY != NULL) && (mxGetClassID(fPtrY) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrY)))
       		 				{
								PtrY = mxGetPr(fPtrY);
								//printf("%e \n", PtrY[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'K':
				zReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrZ = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrZ != NULL) && (mxGetClassID(fPtrZ) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrZ)))
       		 				{
								PtrZ = mxGetPr(fPtrZ);
								//printf("%e \n", PtrZ[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			default:
				defaultReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtr = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtr != NULL) && (mxGetClassID(fPtr) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtr)))
       		 				{
								PtrFinal = mxGetPr(fPtr);
								//printf("%e \n", PtrFinal[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
		}
	}

	// To get the column number of the return
	if (wReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrW;
		}
	}
	if (xReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrX;
		}
	}
	if (yReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrY;
		}
	}
	if (zReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrZ;
		}
	}
	if (defaultReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrFinal;
		}
	}
	if (columnNumber > 0)
	{
		printf("column number too high \n");
		EXIT_FAILURE;
	}

	//sets the matrix array to be equal to an array of doubles
	for (x = 0; x < size_q1; x++) {
		for (y = 0; y < size_q2; y++) {
			matrix_array[x][y] = PtrReturn[ y*(size_q1) +  x];
		}
	}

    if (matClose(mfPtr) != 0) {
        printf("Error closing file %s\n", file.c_str());
    }

	mxDestroyArray(aPtr);

#endif

}







/**
* Read matrix data from file with grid,
* Checks if the matrix data in the file is the same as the grids that were sent in with error < 1e-8,
* if not within error range will signal error and exit
*
* Overloaded readFromFile function
* \param filename - file to read grids from
* \param grids x,y - checks grids data against the file data
*/
template<class T>
void Matrix2D<T>::readFromFile(const std::string& filename, const Matrix2D<T>& grid_x, const Matrix2D<T>& grid_y) {
	std::string inBuf;
	double loaded_x, loaded_y;
	if (!initialized) {
		printf("MATRIX_ERROR: Using un-itialized matrix");
		exit(EXIT_FAILURE);
	} else {
		this->name = filename;
		ifstream input(filename.c_str());
		//if (input != NULL && !input.eof()) {
		if (input.is_open() && !input.eof()) {
			// Skipping first two lines.
			input >> inBuf;
			//	!!!		while (strcmp(strupr((char *) inBuf.c_str()), "ZONE") != 0 && !input.eof() ) {
			while (strcasecmp(inBuf.c_str(), "ZONE") != 0 && !input.eof() ) {
				input >> inBuf;
			}
			// read to the end of the line with 'zone'
			input.ignore(9999, '\n');
			for (size_t i1 = 0; i1 < size_q1; i1++) {
				for (size_t i2 = 0; i2 < size_q2; i2++) {
					input >> loaded_x >> loaded_y;
					// check if grid is the same
					if (fabs(log10(loaded_x) - log10(grid_x[i1][i2])) > err || fabs(log10(loaded_y) - log10(grid_y[i1][i2])) > err) {
						printf("MATRIX_LOAD_GRID_ERR: Loading %s: grid mismatch.\nLoaded: %e, %e\nGrid: %e, %e\n", filename.c_str(), loaded_x, loaded_y, grid_x[i1][i2], grid_y[i1][i2]);
						exit(EXIT_FAILURE);
					} else {
						input >> matrix_array[i1][i2];
						// skip till to the end of the line
						input.ignore(9999, '\n');
					}
				}
			}
		} else {
			printf("MATRIX_LOAD_ERROR: Error reading file %s.\n", filename.c_str());
			exit(EXIT_FAILURE);
		}
		input.close();
	}
}




/**
* Function for reading from matlab file in 2-dimensions
* Will check the variables in the order they are saved in matlab, thus (P R V K Var) should be the standard
* The variables will be checked against the input grid parameters in order to make sure the right variables/values are being loaded
* This is the same as the readFromFile() function although only compatible with .mat files instead of .plt or other text files
*/
template<class T>
void Matrix2D<T>::readFromMatlabFile(const std::string& file , const Matrix2D<T>& grid_x, const Matrix2D<T>& grid_y)
{

#if (MATLAB_CAPABLE)

	Logger::message << "Reading " << file << ": " << endl;

	MATFile *mfPtr; /* MAT-file pointer */
    mxArray *aPtr;  /* mxArray pointer */
	double *PtrW = nullptr; /* pointer to data */
	double *PtrX = nullptr; /* pointer to data */
    double *PtrY = nullptr; /* pointer to data */
    double *PtrZ = nullptr; /* pointer to data */
	double *PtrFinal = nullptr; /* pointer to data */
	double *PtrL = nullptr; /* pointer to data */
	std::string arr; /*name of variable*/
	std::string field = "arr"; // name of field
	mwSize nElements;       /* number of elements in array */
    mwIndex eIdx;           /* element index */
    const mxArray *fPtr;    /* field pointer */
	const mxArray *fPtrW;    /* field pointer */
	const mxArray *fPtrX;    /* field pointer */
	const mxArray *fPtrY;    /* field pointer */
	const mxArray *fPtrZ;    /* field pointer */
	size_t x, y;
	const char* name;		/* for getting variable names */
	char nameTemp;		/* for getting variable names */
	
	bool wReached = false;
	bool xReached = false;
	bool yReached = false;
	bool zReached = false;
	bool defaultReached = false;

	mfPtr = matOpen(file.c_str(), "r");
   	if (mfPtr == NULL) {
       	printf("Error opening file %s\n", file.c_str());
	}

	for (int i= 0; i < 3; i++)
	{

		aPtr = matGetNextVariableInfo(mfPtr, &name);
		std::string temp = name;
		if (temp.length() == 1)
		{
			nameTemp = temp[0];
		}
		else
		{
			if (temp.substr(1, 1) == "_")
			{
				nameTemp = temp[0];
			}
			else if ((temp.substr(0, 2) == "pc"))
			{
				nameTemp = 'V';
			}
			else if ((temp.substr(0, 2) == "al"))
			{
				nameTemp = 'K';
			}
			else
			{
				nameTemp = '!';
			}
		}
		switch (nameTemp)
		{
			case 'P':
				wReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrW = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrW != NULL) && (mxGetClassID(fPtrW) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrW)))
       		 				{
								PtrW = mxGetPr(fPtrW);
								//printf("%e \n", PtrW[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'R':
				xReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrX = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrX != NULL) && (mxGetClassID(fPtrX) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrX)))
       		 				{
								PtrX = mxGetPr(fPtrX);
								//printf("%e \n", PtrX[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'L':
				xReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrX = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrX != NULL) && (mxGetClassID(fPtrX) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrX)))
       		 				{
								PtrX = mxGetPr(fPtrX);
								PtrL = mxGetPr(fPtrX);
								//printf("%e \n", PtrX[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'V':
				yReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrY = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrY != NULL) && (mxGetClassID(fPtrY) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrY)))
       		 				{
								PtrY = mxGetPr(fPtrY);
								//printf("%e \n", PtrY[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'K':
				zReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrZ = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrZ != NULL) && (mxGetClassID(fPtrZ) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrZ)))
       		 				{
								PtrZ = mxGetPr(fPtrZ);
								//printf("%e \n", PtrZ[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			default:
				defaultReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtr = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtr != NULL) && (mxGetClassID(fPtr) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtr)))
       		 				{
								PtrFinal = mxGetPr(fPtr);
								//printf("%e \n", PtrFinal[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
		}
	}

	if (!wReached)
	{
		if (!xReached)
		{
			PtrX = PtrY;
			PtrY = PtrZ;
		}
		if (!yReached)
		{
			PtrY = PtrZ;
		}
	}
	else if (!xReached)
	{
		if (!yReached)
		{
			PtrX = PtrW;
			PtrY = PtrZ;
		}
		if (!zReached)
		{
			PtrX = PtrW;
		}
	}
	else
	{
		PtrY = PtrX;
		PtrX = PtrW;
	}




	//sets the matrix array to be equal to an array of doubles
	for (x = 0; x < size_q1; x++) {
		for (y = 0; y < size_q2; y++) {
			if (fabs(log10(PtrX[y*(size_q1) +  x]) - log10(grid_x[x][y])) > err || fabs(log10(PtrY[ y*(size_q1) +  x]) - log10(grid_y[x][y])) > err) {
			 			printf("MATRIX_LOAD_GRID_ERR: Loading %s: grid mismatch [%zu, %zu].\nLoaded: %e, %e\nGrid: %e, %e\n", file.c_str(), x, y, PtrX[y*(size_q1) +  x],PtrY[ y*(size_q1) +  x], grid_x[x][y], grid_y[x][y]);
						//printf("grid error \n");
						exit(EXIT_FAILURE);
			}
			matrix_array[x][y] = PtrFinal[ y*(size_q1) +  x];
		}
	}

    if (matClose(mfPtr) != 0) {
        printf("Error closing file %s\n", file.c_str());
    }

	mxDestroyArray(aPtr);

#endif
}





/**
* Make x-slice of 2d matrix = 1d matrix.
*/
template<class T>
Matrix1D<T> Matrix2D<T>::xSlice(size_t p_x) const {
	if (p_x > this->size_q1) {
		printf("MATRIX_ERROR: xSlice outside or array boundary: %zu > %zu", p_x, this->size_q1);
		exit(EXIT_FAILURE);
	}

	Matrix1D<T> tmp(this->size_q2);
	// !!! tmp.name = this->name + "_slice";
	for (size_t i2 = 0; i2 < size_q2; i2++) {
		tmp[i2] = matrix_array[p_x][i2];
	}
	return tmp;
}

/**
* Take an x-slice of 2d matrix and write to out.
*/
template<class T>
void Matrix2D<T>::xSlice(Matrix1D<T>& out, size_t p_x) const {
	if (p_x > this->size_q1) {
		printf("MATRIX_ERROR: xSlice outside or array boundary: %zu > %zu", p_x, this->size_q1);
		exit(EXIT_FAILURE);
	}

	for (size_t i2 = 0; i2 < size_q2; i2++) {
		out[i2] = matrix_array[p_x][i2];
	}
}

/**
* Make y-slice of 2d matrix = 1d matrix.
*/
template<class T>
Matrix1D<T> Matrix2D<T>::ySlice(size_t p_y) const {
	if (p_y > this->size_q2) {
		printf("MATRIX_ERROR: ySlice outside or array boundary: %zu > %zu", p_y, this->size_q2);
		exit(EXIT_FAILURE);
	}

	Matrix1D<T> tmp(this->size_q1);
	// !!! tmp.name = this->name + "_slice";
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		tmp[i1] = matrix_array[i1][p_y];
	}
	return tmp;
}

/**
* Take a y-slice of 2d matrix and write to out.
*/
template<class T>
void Matrix2D<T>::ySlice(Matrix1D<T>& out, size_t p_y) const {
	if (p_y > this->size_q2) {
		printf("MATRIX_ERROR: ySlice outside or array boundary: %zu > %zu", p_y, this->size_q2);
		exit(EXIT_FAILURE);
	}

	for (size_t i1 = 0; i1 < size_q1; i1++) {
		out[i1] = matrix_array[i1][p_y];
	}
}
//////////////////////////////////////////////////
//
//	Matrix 3D
//
//////////////////////////////////////////////////

/**
* Constructor.
* Allocate memory.
*/
template<class T>
Matrix3D<T>::Matrix3D(size_t size_q1, size_t size_q2, size_t size_q3) {
	initialized = false;
	// allocating memory
	AllocateMemory(size_q1, size_q2, size_q3);
}

/**
* Constructor.
* Create matrix equal to Matrix M.
*
* \param &M - Matrix M.
*/
template<class T>
Matrix3D<T>::Matrix3D( const Matrix3D<T> &M ) {
	initialized = false;
	this->operator = (M);
}

/**
* Destructor
*/
template<class T>
Matrix3D<T>::~Matrix3D() {
	if (initialized) free_matrix<T>(matrix_array, size_q1, size_q2);
}

/**
* Allocating memory and filling it with zero-values.
*/
template<class T>
void Matrix3D<T>::AllocateMemory(size_t size_q1, size_t size_q2, size_t size_q3) {
	this->size_q1 = size_q1;
	this->size_q2 = size_q2;
	this->size_q3 = size_q3;
	num_elements = size_q1 * size_q2 * size_q3;
	matrix_array = matrix<T>(size_q1, size_q2, size_q3);
	plane_array = matrix_array[0][0];
	initialized = true;
#ifdef DEBUG_MODE
	// should not initialize matrix with zeros, it can slow the code greatly in some cases
	for (size_t i = 0; i < size_q1; i++)
		for (size_t j = 0; j < size_q2; j++)
			for (size_t k = 0; k < size_q3; k++)
				matrix_array[i][j][k] = 0;
#endif
}

/**
* Operator [i1], returns pointer to 2D array. Next [i2][k] can be applied to return value.
* If  DEBUG_MODE defined, check if matrix has been initialized.
*
* \param i - number of element to return
*/
template<class T>
inline T** Matrix3D<T>::operator[] (size_t i1) {
#ifdef DEBUG_MODE
	if (!initialized) {
		printf("MATRIX_ERROR: Using unitialized matrix");
		exit(EXIT_FAILURE);
	}
#endif
	// return i-th pointer
	return matrix_array[i1];
}

/**
* Operator (x, y, z), returns value of element [i1][i2][i3].
* If  DEBUG_MODE defined, check if matrix has been initialized.
*
*/
template<class T>
inline T& Matrix3D<T>::operator() (int x, int y, int z) {
#ifdef DEBUG_MODE
	if (!initialized) {
		printf("MATRIX_ERROR: Using not initialized matrix");
		exit(EXIT_FAILURE);
	}
	if ((x < 0 || x > size_q1-1) || (y < 0 || y > size_q2-1) || (z < 0 || z > size_q3-1)) {
		printf("MATRIX_ERROR: Index is out of bound");
		exit(EXIT_FAILURE);
	}
#endif
	// return (x,y,z) value
	return plane_array[(x*size_q2 + y)*size_q3 + z];
}

/**
* Makes matrix equal to Matrix M.
*
* \param &M - Matrix M.
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::operator= (const Matrix3D<T> &M) {
	// check if not LHS and RHS matrix are the same
	if (this!=&M) {
		// check if RHS matrix is initialized
		if (M.initialized) {
			// free LHS matrix, if it is initialized
			if (initialized && (size_q1 != M.size_q1 || size_q2 != M.size_q2 || size_q3 != M.size_q3)) {
				free_matrix<T>(matrix_array, this->size_q1, this->size_q2);
				initialized = false;
			}
			if (!initialized) {
				this->size_q1 = M.size_q1;
				this->size_q2 = M.size_q2;
				this->size_q3 = M.size_q3;
				// allocating memory for LHS matrix
				this->AllocateMemory(this->size_q1, this->size_q2, this->size_q3);
			}
			// !!! XXX this->name = M.name;
			//			for (int x = 0; x < this->size_q1; x++)
			//				for (int y = 0; y < this->size_q2; y++)
			//					memcpy( matrix_array[i1][i2], M.matrix_array[i1][i2], this->size_q3 * sizeof( T ) );
			// fast values copying as memory range
			memcpy( matrix_array[0][0], M.matrix_array[0][0], this->size_q1 * this->size_q2 * this->size_q3 * sizeof( T ) );
		} else {
			this->initialized = false;
		}
	}
	return *this;
}

/**
* Makes 3D matrix from 2D matrix.
* The 3rd dimension makes equal to 1.
*/
/* template<class T>
Matrix3D<T>& Matrix3D<T>::operator= (const Matrix2D<T> &M) {
double val;
if (M.initialized) {
if (initialized) free_matrix<T>(matrix_array, this->size_q1, this->size_q2);
this->size_q1 = 1;
this->size_q2 = M.size_q1;
this->size_q3 = M.size_q2;
//this->name = M.name;
this->AllocateMemory(this->size_q1, this->size_q2, this->size_q3);
for (int y = 0; y < this->size_q2; y++)
//memcpy( matrix_array[0][i2], M.matrix_array[i2], this->size_q3 * sizeof( T ) );
for (int z = 0; z < this->size_q3; z++) {
this->matrix_array[0][i2][i3] = M[i2][i3];
}
} else {
this->initialized = false;
}
return *this;
}*/

/**
* Makes Matrix equal to value Val.
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::operator= (const T Val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] = Val;
	}
	return *this;
}


/**
* Matrix summation, result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::operator+= (const Matrix3D<T> &M) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] += M.plane_array[i];
	}
	return *this;
}

/**
* Matrix subtraction, result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::operator-= (const Matrix3D<T> &M) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] -= M.plane_array[i];
	}
	return *this;
}

/**
* Multiplication to a value. Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::operator*= (const T Val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] *= Val;
	}
	return *this;
}

/**
* Division by a value. Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::operator/= (const T Val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] /= Val;
	}
	return *this;
}

/**
* Summation with a value. Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::operator+= (const T Val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] += Val;
	}
	return *this;
}

/**
* Subtraction of a value. Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::operator-= (const T Val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] -= Val;
	}
	return *this;
}

/**
* Multiplication between each element of the matrices (not a matrix multiplication). Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::times_equal (const Matrix3D<T> &M) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] *= M.plane_array[i];
	}
	return *this;
}

/**
* Division of each element of one matrices to the element of another. Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::divide_equal (const Matrix3D<T> &M) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] /= M.plane_array[i];
	}
	return *this;
}

/**
* Add each element of the matrix to corresponds element of matrix M.
*/
template<class T>
inline Matrix3D<T> Matrix3D<T>::operator+ (const Matrix3D<T> &M) const {
	Matrix3D<T> Tmp(size_q1, size_q2, size_q3);
	for(size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] + M.plane_array[i];
	}
	return Tmp;
}

/**
* Substract each element of the matrix to corresponds element of matrix M.
*/
template<class T>
inline Matrix3D<T> Matrix3D<T>::operator- (const Matrix3D<T> &M) const {
	Matrix3D<T> Tmp(size_q1, size_q2, size_q3);
	for(size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] - M.plane_array[i];
	}
	return Tmp;
}


/**
* Multiply each element of the matrix to Val, save result to a new matrix.
*/
template<class T>
inline Matrix3D<T> Matrix3D<T>::operator* (const T Val) const {
	Matrix3D<T> Tmp(size_q1, size_q2, size_q3);
	for(size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] * Val;
	}
	return Tmp;
}


/**
* Divide each element of the matrix to Val, save result to a new matrix.
*/
template<class T>
inline Matrix3D<T> Matrix3D<T>::operator/ (const T Val) const {
	Matrix3D<T> Tmp(size_q1, size_q2, size_q3);
	for(size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] / Val;
	}
	return Tmp;
}

/**
* Multiply each element of the matrix to corresponds element of matrix M.
*/
template<class T>
inline Matrix3D<T> Matrix3D<T>::times (const Matrix3D<T> &M) const {
	Matrix3D<T> Tmp(size_q1, size_q2, size_q3);
	for(size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] * M.plane_array[i];
	}
	return Tmp;
}

/**
* Divide each element of the matrix to corresponds element of matrix M.
*/
template<class T>
inline Matrix3D<T> Matrix3D<T>::divide (const Matrix3D<T> &M) const {
	Matrix3D<T> Tmp(size_q1, size_q2, size_q3);
	for(size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] / M.plane_array[i];
	}
	return Tmp;
}


/**
* Write matrix to file.
* File has two header lines.
*/
template<class T>
void Matrix3D<T>::writeToFile(const std::string& filename, const std::string& info) const {
	ofstream output(filename.c_str());
	output << "VARIABLES = \""<< ((this->name!="")?this->name:"f") <<"\" "<< endl;
	output << "ZONE T=\"" << ((info=="")?filename:info) << "\", I=" << size_q3 << ", J=" << size_q2 << ", K=" << size_q1 << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		for (size_t i2 = 0; i2 < size_q2; i2++) {
			for (size_t i3 = 0; i3 < size_q3; i3++) {
				output << matrix_array[i1][i2][i3] << endl;
			}
		}
	}
	output.close();
}

/**
* Write matrix to file, using 3 other matrixes as a grid (simply - write all 4 matrixes to the file).
* File has two header lines.
*/
template<class T>
void Matrix3D<T>::writeToFile(const std::string& filename, const Matrix3D<T>& grid_x, const Matrix3D<T>& grid_y, const Matrix3D<T>& grid_z) const {
	ofstream output(filename.c_str());
	//if (output==NULL && (filename.find("Debug") == std::string::npos)) {
	  if (!output.is_open() && (filename.find("Debug") == std::string::npos)) {
		printf("FILE: Unable to output file: %s", filename.c_str());
		exit(EXIT_FAILURE);
	}
	output << "VARIABLES = \"" << ((grid_x.name!="")?grid_x.name:"x") << "\", \"" << ((grid_y.name!="")?grid_y.name:"y") << "\", \"" << ((grid_z.name!="")?grid_z.name:"z") << "\", \"" << ((this->name!="")?this->name:"f") << "\" "<< endl;
	output << "ZONE T=\"" << filename << "\", I=" << size_q3 << ", J=" << size_q2 << ", K=" << size_q1 << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		for (size_t i2 = 0; i2 < size_q2; i2++) {
			for (size_t i3 = 0; i3 < size_q3; i3++) {
				output << "\t" << grid_x[i1][i2][i3] << "\t" << grid_y[i1][i2][i3] << "\t" << grid_z[i1][i2][i3] << "\t" << matrix_array[i1][i2][i3] << endl;
			}
		}
	}
	output.close();
}

/**
* Read matrix data from file with grid, by column
*
* Overloaded readFromFile function
* \param filename - file to read grids from
* \param read_column - read up to this column from file
*/
template<class T>
void Matrix3D<T>::readFromFile(const std::string& filename, int read_column) {
	std::string inBuf;
	int column;
	if (!initialized) {
		printf("MATRIX_ERROR: Using un-itialized matrix");
		exit(EXIT_FAILURE);
	} else {
		this->name = filename;
		ifstream input(filename.c_str());
		//if (input != NULL && !input.eof()) {
		if (input.is_open() && !input.eof()) {
			// Skipping first two lines.
			input >> inBuf;
			// !!!			while (strcmp(strupr((char *) inBuf.c_str()), "ZONE") != 0 && !input.eof() ) {
			while (strcasecmp(inBuf.c_str(), "ZONE") != 0 && !input.eof() )
				input >> inBuf;
			// read to the end of the line with 'zone'
			input.ignore(9999, '\n'); // read till the end of the line
			//getline(input, inBuf);  // read till the end of the line
			for (size_t i1 = 0; i1 < size_q1; i1++) {
				for (size_t i2 = 0; i2 < size_q2; i2++) {
					for (size_t i3 = 0; i3 < size_q3; i3++) {
						for (column = 1; column < read_column; column++) input >> inBuf; // skip first columns
						input >> matrix_array[i1][i2][i3];
						//getline(input, inBuf);  // read till the end of the line
						input.ignore(9999, '\n'); // read till the end of the line
					}
				}
			}
		} else {
			printf("MATRIX_LOAD_ERROR: Error reading file %s.\n", filename.c_str());
			exit(EXIT_FAILURE);
		}
		input >> inBuf;
		if (!input.eof()) {
			printf("Check grid size - %s.\n", filename.c_str());
			exit(EXIT_FAILURE);
		}
		input.close();
	}
}



/**
* Function for reading from matlab file in 3-dimensions
* Will check the variables, order them in (P, R/L, V, K, Val) format and then set matrix_array to be the variable with the corresponding column number
* This is the same as the readFromFile() function although only compatible with .mat files instead of .plt or other text files
*/
template<class T>
void Matrix3D<T>::readFromMatlabFile(const std::string& file , int columnNumber)
{

#if (MATLAB_CAPABLE)

	Logger::message << "Reading " << file << ": " << endl;

	MATFile *mfPtr; /* MAT-file pointer */
    mxArray *aPtr;  /* mxArray pointer */
    double *realPtr = nullptr; /* pointer to data */
	double *PtrW = nullptr; /* pointer to data */
	double *PtrX = nullptr; /* pointer to data */
    double *PtrY = nullptr; /* pointer to data */
    double *PtrZ = nullptr; /* pointer to data */
	double *PtrFinal = nullptr; /* pointer to data */
	double *PtrReturn = nullptr; /* pointer to data */
	double *PtrL = nullptr; /* pointer to data */
	std::string arr; /*name of variable*/
	std::string field = "arr"; // name of field
	mwSize nElements;       /* number of elements in array */
    mwIndex eIdx;           /* element index */
    const mxArray *fPtr;    /* field pointer */
	const mxArray *fPtrW;    /* field pointer */
	const mxArray *fPtrX;    /* field pointer */
	const mxArray *fPtrY;    /* field pointer */
	const mxArray *fPtrZ;    /* field pointer */
    size_t x,y,z; 			/* for index*/
	const char* name;		/* for getting variable names */
	char nameTemp;		/* for getting variable names */
	bool wReached = false;
	bool xReached = false;
	bool yReached = false;
	bool zReached = false;
	bool defaultReached = false;

	mfPtr = matOpen(file.c_str(), "r");
   	if (mfPtr == NULL) {
       	printf("Error opening file %s\n", file.c_str());
	}

	for (int i= 0; i < 4; i++)
	{
		aPtr = matGetNextVariableInfo(mfPtr, &name);
		std::string temp = name;
		if (temp.length() == 1)
		{
			nameTemp = temp[0];
		}
		else
		{
			if (temp.substr(1, 1) == "_")
			{
				nameTemp = temp[0];
			}
			else if ((temp.substr(0, 2) == "pc"))
			{
				nameTemp = 'V';
			}
			else if ((temp.substr(0, 2) == "al"))
			{
				nameTemp = 'K';
			}
			else
			{
				nameTemp = '!';
			}
		}
		switch (nameTemp)
		{
			case 'P':
				wReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrW = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrW != NULL) && (mxGetClassID(fPtrW) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrW)))
       		 				{
								PtrW = mxGetPr(fPtrW);
								//printf("%e \n", PtrW[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'R':
				xReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrX = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrX != NULL) && (mxGetClassID(fPtrX) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrX)))
       		 				{
								PtrX = mxGetPr(fPtrX);
								//printf("%e \n", PtrX[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'L':
				xReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrX = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrX != NULL) && (mxGetClassID(fPtrX) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrX)))
       		 				{
								PtrX = mxGetPr(fPtrX);
								PtrL = mxGetPr(fPtrX);
								//printf("%e \n", PtrX[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'V':
				yReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrY = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrY != NULL) && (mxGetClassID(fPtrY) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrY)))
       		 				{
								PtrY = mxGetPr(fPtrY);
								//printf("%e \n", PtrY[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'K':
				zReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrZ = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrZ != NULL) && (mxGetClassID(fPtrZ) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrZ)))
       		 				{
								PtrZ = mxGetPr(fPtrZ);
								//printf("%e \n", PtrZ[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			default:
				defaultReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtr = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtr != NULL) && (mxGetClassID(fPtr) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtr)))
       		 				{
								PtrFinal = mxGetPr(fPtr);
								//printf("%e \n", PtrFinal[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
		}
	}


	// To get the column number of the return
	if (wReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrW;
		}
	}
	if (xReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrX;
		}
	}
	if (yReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrY;
		}
	}
	if (zReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrZ;
		}
	}
	if (defaultReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrFinal;
		}
	}
	if (columnNumber > 0)
	{
		printf("column number too high \n");
		EXIT_FAILURE;
	}




	//sets the matrix array to be equal to an array of doubles
	for (x = 0; x < size_q1; x++) {
		for (y = 0; y < size_q2; y++) {
			for (z = 0; z < size_q3; z++) {
				matrix_array[x][y][z] = PtrReturn[z*(size_q1 * size_q2) + y*(size_q1) +  x];
			}
		}
	}

    if (matClose(mfPtr) != 0) {
        printf("Error closing file %s\n", file.c_str());
    }

	mxDestroyArray(aPtr);

#endif
}




/**
* Read matrix data from file with grid,
* Checks if the matrix data in the file is the same as the grids that were sent in with error < 1e-8,
* if not within error range will signal error and exit
*
* Overloaded readFromFile function
* \param filename - file to read grids from
* \param grids x,y,z - checks grids data against the file data
*/
template<class T>
void Matrix3D<T>::readFromFile(const std::string& filename, const Matrix3D<T>& grid_x, const Matrix3D<T>& grid_y, const Matrix3D<T>& grid_z) {
	std::string inBuf;
	double loaded_x, loaded_y, loaded_z;

	if (!initialized) {
		printf("MATRIX_ERROR: Using unitialized matrix");
		exit(EXIT_FAILURE);
	} else {
		this->name = filename;
		ifstream input(filename.c_str());
		//if (input != NULL && !input.eof()) {
		if (input.is_open() && !input.eof()) {
			// Skipping first two lines.
			input >> inBuf;
			// !!!			while (strcmp(_strupr((char *) inBuf.c_str()), "ZONE") != 0 && !input.eof() ) {
			// !!!			while (strcasecmp(inBuf.c_str(), "ZONE") != 0 && !input.eof() ) {
			while (strcmp(inBuf.c_str(), "ZONE") != 0 && !input.eof() ) {
				input >> inBuf;
			}
			// read to the end of the line with 'zone'
			input.ignore(9999, '\n');
			for (size_t i1 = 0; i1 < size_q1; i1++) {
				for (size_t i2 = 0; i2 < size_q2; i2++) {
					for (size_t i3 = 0; i3 < size_q3; i3++) {
						//input >> grid_x[i1][i2][i3] >> grid_y[i1][i2][i3] >> grid_z[i1][i2][i3] >> matrix_array[i1][i2][i3];
						// skip till the end of the line
						//input.ignore(9999, '\n');

						 input >> loaded_x >> loaded_y >> loaded_z;
						 // check if grid is the same
						 if (fabs(log10(loaded_x) - log10(grid_x[i1][i2][i3])) > err || fabs(log10(loaded_y) - log10(grid_y[i1][i2][i3])) > err || fabs(log10(loaded_z) - log10(grid_z[i1][i2][i3])) > err) {
							printf("MATRIX_LOAD_GRID_ERR: Loading %s: grid mismatch [%zu, %zu, %zu].\nLoaded: %e, %e, %e\nGrid: %e, %e, %e\n", filename.c_str(), i1, i2, i3, loaded_x, loaded_y, loaded_z, grid_x[i1][i2][i3], grid_y[i1][i2][i3], grid_z[i1][i2][i3]);
							exit(EXIT_FAILURE);
						} else {
							input >> matrix_array[i1][i2][i3];
							// skip till to the end of the line
							input.ignore(9999, '\n');
						}
					}
				}
			}
		} else {
			printf("MATRIX_LOAD_ERROR: Error reading file %s.\n", filename.c_str());
			exit(EXIT_FAILURE);
		}
		input.close();
	}
}







/**
* Function for reading from matlab file in 3-dimensions
* Will check the variables in the order they are saved in matlab, thus (P R V K Var) should be the standard
* The variables will be checked against the input grid parameters in order to make sure the right variables/values are being loaded
* This is the same as the readFromFile() function although only compatible with .mat files instead of .plt or other text files
*/
template<class T>
void Matrix3D<T>::readFromMatlabFile(const std::string& file , const Matrix3D<T>& grid_x, const Matrix3D<T>& grid_y, const Matrix3D<T>& grid_z)
{

#if (MATLAB_CAPABLE)

	Logger::message << "Reading " << file << ": " << endl;

	MATFile *mfPtr; /* MAT-file pointer */
    mxArray *aPtr;  /* mxArray pointer */
    double *realPtr = nullptr; /* pointer to data */
	double *PtrW = nullptr; /* pointer to data */
	double *PtrX = nullptr; /* pointer to data */
    double *PtrY = nullptr; /* pointer to data */
    double *PtrZ = nullptr; /* pointer to data */
	double *PtrFinal = nullptr; /* pointer to data */
	double *PtrL = nullptr; /* pointer to data */
	std::string arr; /*name of variable*/
	std::string field = "arr"; // name of field
	mwSize nElements;       /* number of elements in array */
    mwIndex eIdx;           /* element index */
    const mxArray *fPtr;    /* field pointer */
	const mxArray *fPtrW;    /* field pointer */
	const mxArray *fPtrX;    /* field pointer */
	const mxArray *fPtrY;    /* field pointer */
	const mxArray *fPtrZ;    /* field pointer */
    size_t x,y,z; 			/* for index*/
	const char* name;		/* for getting variable names */
	char nameTemp;		/* for getting variable names */
	bool wReached = false;
	bool xReached = false;
	bool yReached = false;
	bool zReached = false;
	bool defaultReached = false;

	mfPtr = matOpen(file.c_str(), "r");
   	if (mfPtr == NULL) {
       	printf("Error opening file %s\n", file.c_str());
	}

	for (int i= 0; i < 4; i++)
	{

		aPtr = matGetNextVariableInfo(mfPtr, &name);
		std::string temp = name;
		if (temp.length() == 1)
		{
			nameTemp = temp[0];
		}
		else
		{
			if (temp.substr(1, 1) == "_")
			{
				nameTemp = temp[0];
			}
			else if ((temp.substr(0, 2) == "pc"))
			{
				nameTemp = 'V';
			}
			else if ((temp.substr(0, 2) == "al"))
			{
				nameTemp = 'K';
			}
			else
			{
				nameTemp = '!';
			}
		}
		switch (nameTemp)
		{
			case 'P':
				wReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrW = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrW != NULL) && (mxGetClassID(fPtrW) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrW)))
       		 				{
								PtrW = mxGetPr(fPtrW);
								//printf("%e \n", PtrW[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'R':
				xReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrX = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrX != NULL) && (mxGetClassID(fPtrX) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrX)))
       		 				{
								PtrX = mxGetPr(fPtrX);
								//printf("%e \n", PtrX[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'L':
				xReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrX = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrX != NULL) && (mxGetClassID(fPtrX) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrX)))
       		 				{
								PtrX = mxGetPr(fPtrX);
								PtrL = mxGetPr(fPtrX);
								//printf("%e \n", PtrX[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'V':
				yReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrY = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrY != NULL) && (mxGetClassID(fPtrY) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrY)))
       		 				{
								PtrY = mxGetPr(fPtrY);
								//printf("%e \n", PtrY[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'K':
				zReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrZ = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrZ != NULL) && (mxGetClassID(fPtrZ) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrZ)))
       		 				{
								PtrZ = mxGetPr(fPtrZ);
								//printf("%e \n", PtrZ[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			default:
				defaultReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtr = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtr != NULL) && (mxGetClassID(fPtr) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtr)))
       		 				{
								PtrFinal = mxGetPr(fPtr);
								//printf("%e \n", PtrFinal[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
		}
	}

	if (!xReached)
	{
		PtrX = PtrW;
	}
	if (!yReached)
	{
		PtrY = PtrX;
		PtrX = PtrW;
	}
	if (!zReached)
	{
		PtrZ = PtrY;
		PtrY = PtrX;
		PtrX = PtrW;
	}


	//sets the matrix array to be equal to an array of doubles
	for (x = 0; x < size_q1; x++) {
		for (y = 0; y < size_q2; y++) {
			for (z = 0; z < size_q3; z++) {
				if (fabs(log10(PtrX[z*(size_q1 * size_q2) + y*(size_q1) +  x]) - log10(grid_x[x][y][z])) > err || fabs(log10(PtrY[z*(size_q1 * size_q2) + y*(size_q1) +  x]) - log10(grid_y[x][y][z])) > err || fabs(log10(PtrZ[z*(size_q1 * size_q2) + y*(size_q1) +  x]) - log10(grid_z[x][y][z])) > err) {
				 			printf("MATRIX_LOAD_GRID_ERR: Loading %s: grid mismatch [%zu, %zu, %zu].\nLoaded: %e, %e, %e\nGrid: %e, %e, %e\n", file.c_str(), x, y, z, PtrX[z*(size_q1 * size_q2) + y*(size_q1) +  x],PtrY[z*(size_q1 * size_q2) + y*(size_q1) +  x],PtrZ[z*(size_q1 * size_q2) + y*(size_q1) +  x], grid_x[x][y][z], grid_y[x][y][z], grid_z[x][y][z]);
							//printf("grid error \n");
							exit(EXIT_FAILURE);
				}
				matrix_array[x][y][z] = PtrFinal[z*(size_q1 * size_q2) + y*(size_q1) +  x];
			}
		}
	}

    if (matClose(mfPtr) != 0) {
        printf("Error closing file %s\n", file.c_str());
    }

	mxDestroyArray(aPtr);

#endif

}


/**
* Write matrix data to binary file

* \param filename - file to read grids from
*/
template<class T>
void Matrix3D<T>::writeToBinaryFile(const std::string& filename) const {
	std::ofstream outputFile(filename, std::ios::binary | std::ios::out);
    if (outputFile.is_open()) {
		int32_t size_array[3] =  { (int32_t)size_q1, (int32_t)size_q2, (int32_t)size_q3 };
		outputFile.write((char*)size_array, 3 * sizeof(int32_t));
		if (!outputFile.good()){
			printf("Writing error");
			exit(EXIT_FAILURE);
		}

		outputFile.write((char*)plane_array, num_elements * sizeof(T));
		if (!outputFile.good()){
			printf("Writing error");
			exit(EXIT_FAILURE);
		}
    } else {
        printf("MATRIX_WRITE_ERROR: Error writing file %s.\n", filename.c_str());
        exit(EXIT_FAILURE);
    }
    outputFile.close();
}




/**
* Read matrix data from binary file
*
* \param filename - file to read grids from
*/
template<class T>
void Matrix3D<T>::readFromBinaryFile(const std::string& filename) {
	if (!initialized) {
		printf("MATRIX_ERROR: Using un-itialized matrix");
		exit(EXIT_FAILURE);
	} else {
		this->name = filename;
		std::ifstream inputFile(filename, std::ios::binary | std::ios::in);
		if (inputFile.is_open()) {
			int32_t buffer_int32[3];
			inputFile.read((char*)buffer_int32, 3 * sizeof(int32_t));
			if (!inputFile.good()){
				printf("Reading error");
				exit(EXIT_FAILURE);
			}

			inputFile.read((char*)plane_array, num_elements * sizeof(T));
			if (!inputFile.good()){
				printf("Reading error");
				exit(EXIT_FAILURE);
			}
		} else {
			printf("MATRIX_LOAD_ERROR: Error reading file %s.\n", filename.c_str());
			exit(EXIT_FAILURE);
		}
		inputFile.close();
	}
}


/**
* Write 3D data to .plt, .pltb or .mat files
* WARNING: writing to 3D matlab array is not available at the moment
*/
template<class T>
void Matrix3D<T>::writeToAnyFile(const std::string& filename, const std::string& io_method, const std::string& info) const {
    std::string ext;
    if (io_method.compare("ascii") == 0){
        ext = ".plt";
        writeToFile(filename + ext, info);
    } else if (io_method.compare("binary") == 0) {
        ext = ".pltb";
        writeToBinaryFile(filename + ext);
    } else if (io_method.compare("matlab") == 0) {
        printf("Writing to 3D arrays to MATLAB filse is not available in the moment");
        //ext = ".mat";
        //writeToMatlabFile(filename + ext);
    } else {
        printf("IO error: unknown io_method");
        exit(EXIT_FAILURE);
    }
}

/**
* Read 3D data from .plt, .pltb or .mat files containing only one column
*/
template<class T>
void Matrix3D<T>::readFromAnyFile(const std::string& filename, const std::string& io_method) {
    std::string ext;
    if (io_method.compare("ascii") == 0){
        ext = ".plt";
        readFromFile(filename + ext, 1);
    } else if (io_method.compare("binary") == 0){
        ext = ".pltb";
        readFromBinaryFile(filename + ext);
    } else if (io_method.compare("matlab") == 0){
        ext = ".mat";
        readFromMatlabFile(filename + ext, 1);
    } else {
        printf("IO error: unknown io_method");
        exit(EXIT_FAILURE);
    }
}

/**
* Read 3D data from .plt, .pltb or .mat files with grids
* WARNING: if io_method == "binary", .pltb file does not contain grid. In this case the function works the same as Matrix3D<T>::readFromAnyFile(std::string, std::string, int)
*/
template<class T>
void Matrix3D<T>::readFromAnyFile(const std::string& filename, const std::string& io_method, const Matrix3D<T>& grid_q1, const Matrix3D<T>& grid_q2, const Matrix3D<T>& grid_q3){
    std::string ext;
    if (io_method.compare("ascii") == 0){
        ext = ".plt";
        readFromFile(filename + ext, grid_q1, grid_q2, grid_q3);
    } else if (io_method.compare("binary") == 0) {
        ext = ".pltb";
        readFromBinaryFile(filename + ext);
    } else if (io_method.compare("matlab") == 0) {
        ext = ".mat";
        readFromMatlabFile(filename + ext, grid_q1, grid_q2, grid_q3);
    } else {
        printf("IO error: unknown io_method");
        exit(EXIT_FAILURE);
    }
}


/**
* Returns corresponding index of 3d matrix if represented as a 1d array
* \param x,y,z - index of element in every dimension for the 3d matrix
*/
template<class T>
inline size_t Matrix3D<T>::index1d(size_t x, size_t y, size_t z) {
	return (x*size_q2 + y)*size_q3 + z;
}

/**
* Return minimum value of the 3d matrix.
* Default value set to 1e99
*/
template<class T>
T Matrix3D<T>::min() const {
	T tmp = std::numeric_limits<T>::infinity();
	for(size_t  i = 0; i < num_elements; i++)
	{
		tmp = std::min(tmp, plane_array[i]);
	}
	return tmp;
}

/**
* Return maximum value of the 3d matrix.
* Default value seet to 0
*/
template<class T>
T Matrix3D<T>::max() const {
	T tmp = -std::numeric_limits<T>::infinity();
	for(size_t  i = 0; i < num_elements; i++)
	{
		tmp = std::max(tmp, plane_array[i]);
	}
	return tmp;
}

/**
* Return absolute maximum value of the 3d matrix.
* Default value set to 0
*/
template<class T>
T Matrix3D<T>::maxabs() const {
	T tmp = 0;
	for(size_t  i = 0; i < num_elements; i++)
	{
		tmp = std::max(tmp, std::abs(plane_array[i]));
	}
	return tmp;
}

/**
* Return absolute value of the 3d matrix.
* Changes every element to a positive value with the same magnitude
*/
template<class T>
Matrix3D<T> Matrix3D<T>::abs() const {
	Matrix3D<T> tmp(size_q1, size_q2, size_q3);
	for(size_t  i = 0; i < num_elements; i++)
	{
		tmp.plane_array[i] = std::abs(plane_array[i]);
	}
	return tmp;
}

/**
* Replace all matrix malues smaller than val with val.
*
* \param val - new minimal value.
*/
template<class T>
Matrix3D<T>& Matrix3D<T>::max_of(T val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] = std::max(plane_array[i], val);
	}
	return *this;
}

/**
* Take x-slice of 3d matrix turning it into a 2d matrix.
* \param p_x - index at which to slice x dimension
*/
template<class T>
Matrix2D<T> Matrix3D<T>::xSlice(size_t p_x) const {
	Matrix2D<T> tmp(this->size_q2, this->size_q3);
	// !!! tmp.name = this->name + "_slice";
	for (size_t i2 = 0; i2 < size_q2; i2++) {
		for (size_t i3 = 0; i3 < size_q3; i3++) {
			tmp[i2][i3] = matrix_array[p_x][i2][i3];
		}
	}
	return tmp;
}

/**
* Take an x-slice of 3d matrix and write it to out
* \param p_x - index at which to slice x dimension
*/
template<class T>
void Matrix3D<T>::xSlice(Matrix2D<T>& out, size_t p_x) const {
	for (size_t i2 = 0; i2 < size_q2; i2++) {
		for (size_t i3 = 0; i3 < size_q3; i3++) {
			out[i2][i3] = matrix_array[p_x][i2][i3];
		}
	}
}

/**
* Take y-slice of 3d matrix turning it into a 2d matrix.
* \param p_y - index at which to slice y dimension
*/
template<class T>
Matrix2D<T> Matrix3D<T>::ySlice(size_t p_y) const {
	Matrix2D<T> tmp(this->size_q1, this->size_q3);
	// !!! tmp.name = this->name + "_slice";
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		for (size_t i3 = 0; i3 < size_q3; i3++) {
			tmp[i1][i3] = matrix_array[i1][p_y][i3];
		}
	}
	return tmp;
}

/**
* Take a y-slice of 3d matrix and write it to out
* \param p_x - index at which to slice x dimension
*/
template<class T>
void Matrix3D<T>::ySlice(Matrix2D<T>& out, size_t p_y) const {
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		for (size_t i3 = 0; i3 < size_q3; i3++) {
			out[i1][i3] = matrix_array[i1][p_y][i3];
		}
	}
}

/**
* Take z-slice of 3d matrix turning it into a 2d matrix.
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix2D<T> Matrix3D<T>::zSlice(size_t p_z) const {
	Matrix2D<T> tmp(this->size_q1, this->size_q2);
	// !!! tmp.name = this->name+ "_slice";
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		for (size_t i2 = 0; i2 < size_q2; i2++) {
			tmp[i1][i2] = matrix_array[i1][i2][p_z];
		}
	}
	return tmp;
}

/**
* Take a z-slice of 3d matrix and write it to out
* \param p_x - index at which to slice x dimension
*/
template<class T>
void Matrix3D<T>::zSlice(Matrix2D<T>& out, size_t p_z) const {
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		for (size_t i2 = 0; i2 < size_q2; i2++) {
			out[i1][i2] = matrix_array[i1][i2][p_z];
		}
	}
}

/**
* Take xy-slice of 3d matrix turning it into a 1d matrix.
* \param p_x - index at which to slice x dimension
* \param p_y - index at which to slice y dimension
*/
template<class T>
Matrix1D<T> Matrix3D<T>::xySlice(size_t p_x, size_t p_y) const {
	Matrix1D<T> tmp(this->size_q3);
	// !!! tmp.name = this->name + "_xySlice";
	for (size_t i3 = 0; i3 < size_q3; i3++) {
		tmp[i3] = matrix_array[p_x][p_y][i3];
	}
	return tmp;
}

/**
* Take an xy-slice of 3d matrix and write it to out.
* \param p_x - index at which to slice x dimension
* \param p_y - index at which to slice y dimension
*/
template<class T>
void Matrix3D<T>::xySlice(Matrix1D<T>& out, size_t p_x, size_t p_y) const {
	for (size_t i3 = 0; i3 < size_q3; i3++) {
		out[i3] = matrix_array[p_x][p_y][i3];
	}
}

/**
* Take yz-slice of 3d matrix turning it into a 1d matrix.
* \param p_y - index at which to slice y dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix1D<T> Matrix3D<T>::yzSlice(size_t p_y, size_t p_z) const {
	Matrix1D<T> tmp(this->size_q1);
	// !!! tmp.name = this->name + "_xzSlice";
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		tmp[i1] = matrix_array[i1][p_y][p_z];
	}
	return tmp;
}

/**
* Take a yz-slice of 3d matrix and write it to out.
* \param p_y - index at which to slice y dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
void Matrix3D<T>::yzSlice(Matrix1D<T>& out, size_t p_y, size_t p_z) const {
	for (size_t i1 = 0; i1 < size_q1; i1++) {
		out[i1] = matrix_array[i1][p_y][p_z];
	}
}

/**
* Take xz-slice of 3d matrix turning it into a 1d matrix.
* \param p_x - index at which to slice x dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix1D<T> Matrix3D<T>::xzSlice(size_t p_x, size_t p_z) const {
	Matrix1D<T> tmp(this->size_q2);
	// !!! tmp.name = this->name + "_xzSlice";
	for (size_t i2 = 0; i2 < size_q2; i2++) {
		tmp[i2] = matrix_array[p_x][i2][p_z];
	}
	return tmp;
}

/**
* Take an xz-slice of 3d matrix and write it to out.
* \param p_x - index at which to slice x dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
void Matrix3D<T>::xzSlice(Matrix1D<T>& out, size_t p_x, size_t p_z) const {
	for (size_t i2 = 0; i2 < size_q2; i2++) {
		out[i2] = matrix_array[p_x][i2][p_z];
	}
}


//////////////////////////////////////////////////
//
//	Matrix 4D
//
//////////////////////////////////////////////////

/**
* Constructor.
* Allocate memory.
*/
template<class T>
Matrix4D<T>::Matrix4D(size_t w_size, size_t x_size, size_t y_size, size_t z_size) {
	initialized = false;
	// allocating memory
	AllocateMemory(w_size, x_size, y_size, z_size);
}

/**
* Constructor.
* Create matrix equal to Matrix M.
*
* \param &M - Matrix M.
*/
template<class T>
Matrix4D<T>::Matrix4D( const Matrix4D<T> &M ) {
	initialized = false;
	this->operator = (M);
}

/**
* Destructor
*/
template<class T>
Matrix4D<T>::~Matrix4D() {
	if (initialized) free_matrix<T>(matrix_array, size_w, size_x, size_y);
	//		delete[] plane_array; - it's get deleted with matrix_array[0][0][0] or something?
}

/**
* Allocating memory and filling it with zero-values.
*/
template<class T>
void Matrix4D<T>::AllocateMemory(size_t w_size, size_t x_size, size_t y_size, size_t z_size) {
	this->size_w = w_size;
	this->size_x = x_size;
	this->size_y = y_size;
	this->size_z = z_size;
	num_elements = size_w * size_x * size_y * size_z;
	matrix_array = matrix<T>(w_size, x_size, y_size, z_size);
	plane_array = matrix_array[0][0][0];
	initialized = true;
#ifdef DEBUG_MODE
	// should not initialize matrix with zeros, it can slow the code greatly in some cases
	for (int w = 0; w < size_w; w++)
		for (int x = 0; x < size_x; x++)
			for (int y = 0; y < size_y; y++)
				for (int z = 0; z < size_z; z++)
					matrix_array[w][x][y][z] = 0;
#endif
}

/**
* Operator [i], returns pointer to 3D array. Next [j][k][l] can be applied to return value.
* If  DEBUG_MODE defined, check if matrix has been initialized.
*
* \param i - number of element to return
*/
template<class T>
inline T*** Matrix4D<T>::operator[] (size_t i) {
#ifdef DEBUG_MODE
	if (!initialized) {
		printf("MATRIX_ERROR: Using unitialized matrix");
		exit(EXIT_FAILURE);
	}
#endif
	// return i-th pointer
	return matrix_array[i];
}

/**
* Operator (w, x, y, z), returns value of element [w][x][y][z].
* If  DEBUG_MODE defined, check if matrix has been initialized.
*
*/
template<class T>
inline T& Matrix4D<T>::operator() (int w, int x, int y, int z) {
#ifdef DEBUG_MODE
	if (!initialized) {
		printf("MATRIX_ERROR: Using not initialized matrix");
		exit(EXIT_FAILURE);
	}
	if ((w < 0 || w > size_w-1) || (x < 0 || x > size_x-1) || (y < 0 || y > size_y-1) || (z < 0 || z > size_z-1)) {
		printf("MATRIX_ERROR: Index is out of bound");
		exit(EXIT_FAILURE);
	}
#endif
	// return (w,x,y,z) value
	return plane_array[((w*size_x + x)*size_y + y)*size_z + z];
}

/**
* Makes current matrix equal to Matrix M and returns the current matrix
*
* \param &M - Matrix M.
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::operator= (const Matrix4D<T> &M) {
	// check if not LHS and RHS matrix are the same
	if (this!=&M) {
		// check if RHS matrix is initialized
		if (M.initialized) {
			// free LHS matrix, if it is initialized
			if (initialized && (size_w != M.size_w || size_x != M.size_x || size_y != M.size_y || size_z != M.size_z)) {
				free_matrix<T>(matrix_array, this->size_w, this->size_x, this->size_y);
				initialized = false;
			}
			if (!initialized) {
				this->size_w = M.size_w;
				this->size_x = M.size_x;
				this->size_y = M.size_y;
				this->size_z = M.size_z;
				// allocating memory for LHS matrix
				this->AllocateMemory(this->size_w, this->size_x, this->size_y, this->size_z);
			}
			// !!! XXX this->name = M.name;
			this->name = std::string(M.name.length(), ' ');
			std::copy ( M.name.begin( ), M.name.end( ), this->name.begin( )) ;

			// fast values copying as memory range
			memcpy( matrix_array[0][0][0], M.matrix_array[0][0][0], this->size_w * this->size_x * this->size_y * this->size_z * sizeof( T ) );
		} else {
			this->initialized = false;
		}
	}
	return *this;
}

/**
* Makes every element in Matrix equal to Val and returns this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::operator= (const T Val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] = Val;
	}
	return *this;
}


/**
* Add each element of the current matrix by the corresponding element in M and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::operator+= (const Matrix4D<T> &M) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] += M.plane_array[i];
	}
	return *this;
}

/**
* Subtract each element of the current matrix by the corresponding element in M and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::operator-= (const Matrix4D<T> &M) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] -= M.plane_array[i];
	}
	return *this;
}

/**
* Multiply each element of the current matrix by Val and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::operator*= (const T Val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] *= Val;
	}
	return *this;
}

/**
* Divide each element of the current matrix by Val and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::operator/= (const T Val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] /= Val;
	}
	return *this;
}

/**
* Add each element of the current matrix by Val and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::operator+= (const T Val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] += Val;
	}
	return *this;
}

/**
* Subtract each element of the current matrix by Val and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::operator-= (const T Val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] -= Val;
	}
	return *this;
}

/**
* Multiply each element of the current matrix by the corresponding element of matrix M and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::times_equal (const Matrix4D<T> &M) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] *= M.plane_array[i];
	}
	return *this;
}

/**
* Divide each element of the current matrix by the corresponding element of matrix M and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::divide_equal (const Matrix4D<T> &M) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] /= M.plane_array[i];
	}
	return *this;
}

/**
* Add each element of the matrix by the corresponding element of matrix M and return new matrix
*/
template<class T>
inline Matrix4D<T> Matrix4D<T>::operator+ (const Matrix4D<T> &M) const {
	Matrix4D<T> Tmp(size_w, size_x, size_y, size_z);
	for(size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] + M.plane_array[i];
	}
	return Tmp;
}

/**
* Substract each element of the matrix by the corresponding element of matrix M and return new matrix
*/
template<class T>
inline Matrix4D<T> Matrix4D<T>::operator- (const Matrix4D<T> &M) const {
	Matrix4D<T> Tmp(size_w, size_x, size_y, size_z);
	for(size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] - M.plane_array[i];
	}
	return Tmp;
}


/**
* Multiply each element of the matrix by Val, save result to a new matrix.
*/
template<class T>
inline Matrix4D<T> Matrix4D<T>::operator* (const T Val) const {
	Matrix4D<T> Tmp(size_w, size_x, size_y, size_z);
	for(size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] * Val;
	}
	return Tmp;
}


/**
* Divide each element of the matrix by Val, save result to a new matrix.
*/
template<class T>
inline Matrix4D<T> Matrix4D<T>::operator/ (const T Val) const {
	Matrix4D<T> Tmp(size_w, size_x, size_y, size_z);
	for(size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] / Val;
	}
	return Tmp;
}

/**
* Multiply each element of the matrix by the corresponding element of matrix M and return new matrix
*/
template<class T>
inline Matrix4D<T> Matrix4D<T>::times (const Matrix4D<T> &M) const {
	Matrix4D<T> Tmp(size_w, size_x, size_y, size_z);
	for(size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] * M.plane_array[i];
	}
	return Tmp;
}

/**
* Divide each element of the current matrix by the corresponding element of matrix M and return new matrix
*/
template<class T>
inline Matrix4D<T> Matrix4D<T>::divide (const Matrix4D<T> &M) const {
	Matrix4D<T> Tmp(size_w, size_x, size_y, size_z);
	for(size_t i = 0; i < num_elements; i++)
	{
		Tmp.plane_array[i] = plane_array[i] / M.plane_array[i];
	}
	return Tmp;
}


/**
* Write matrix to file.
* File has two header lines.
*/
template<class T>
void Matrix4D<T>::writeToFile(const std::string& filename, const std::string& info) const {
	ofstream output(filename.c_str());
	output << "VARIABLES = \""<< ((this->name!="")?this->name:"f") <<"\" "<< endl;
	output << "ZONE T=\"" << ((info=="")?filename:info) << "\", I=" << size_z << ", J=" << size_y << ", K=" << size_x << ", L=" << size_w << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (size_t w = 0; w < size_w; w++)
		for (size_t x = 0; x < size_x; x++)
			for (size_t y = 0; y < size_y; y++)
				for (size_t z = 0; z < size_z; z++)
					output << matrix_array[w][x][y][z] << endl;
	output.close();
}








#if (MATLAB_CAPABLE)



/**
* Packaging function that turns a Matrix4D into a variable that can be stored into a .mat file.
* Used in conjunction with writeToMatlabFile() to save the 4 grid variables and val into a single .mat file
*/
template<class T>
mxArray* Matrix4D<T>::createStructMatrix(const std::string& filename, const std::string& info) const
{

	size_t w,x,y,z;
    // create a struct array with 7 fields
    const char *fieldnames[7] = {"arr", "time", "size", "size1", "size2" , "size3" , "size4"};
	// create a 1x1 struct that will hold the array of values and the time info, as well as the size of the dimensions
	mxArray *s = mxCreateStructMatrix(1, 1, 7, fieldnames);
	mwSize size_W = size_w;
	mwSize size_X = size_x;
	mwSize size_Y = size_y;
	mwSize size_Z = size_z;



    // fill struct fields
    for (mwIndex i=0; i<7; i++) {

		// For the arr variable in the 1x1 struct
		// Holds the contents of matrix_array
		if (i== 0)
		{
			// create array of doubles and give it the data in matrix_array
			const int ndim = 4;
			mwSize dims[ndim] = {size_W , size_X, size_Y, size_Z};
			mxArray *aPtr = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
			if (aPtr == NULL) {
    			printf("Unable to create double mxArray.\n");
    			exit(EXIT_FAILURE);
  			}
			double *data = mxGetPr(aPtr);
			for (w = 0; w < size_w; w++) {
				for (x = 0; x < size_x; x++) {
					for (y = 0; y < size_y; y++) {
						for (z = 0; z < size_z; z++) {
							data[z*(size_w * size_x * size_y) + y*(size_w * size_x) +  x*(size_w) + w] = matrix_array[w][x][y][z];
						}
					}
				}
			}
			mxSetField(s, 0, fieldnames[i], aPtr);
		}

		// For the time variable in the 1x1 struct
		if ( i == 1 )
		{
			// Create a std::string that has the inputted time
			mxArray *aPtr = mxCreateDoubleMatrix(1,1,mxREAL);
			double *data = mxGetPr(aPtr);
			double temp = 0.0;
			std::stringstream ss;
			ss << info;
			ss >> temp;
			data[0] = temp;
			mxSetField(s, 0, fieldnames[i], aPtr);
		}

		// For the size variable in the 1x1 struct
		if ( i == 2 )
		{
			mxArray *aPtr = mxCreateDoubleMatrix(1,1,mxREAL);
			double *data = mxGetPr(aPtr);
			data[0] = 1;
			mxSetField(s, 0, fieldnames[i], aPtr);
		}

		// size of dimension 1 or w
		if ( i == 3 )
		{
			mxArray *aPtr = mxCreateDoubleMatrix(1,1,mxREAL);
			double *data = mxGetPr(aPtr);
			data[0] = static_cast<double>(size_w);
			mxSetField(s, 0, fieldnames[i], aPtr);
		}

		// size of dimension 2 or x
		if ( i == 4 )
		{
			mxArray *aPtr = mxCreateDoubleMatrix(1,1,mxREAL);
			double *data = mxGetPr(aPtr);
			data[0] = static_cast<double>(size_x);
			mxSetField(s, 0, fieldnames[i], aPtr);
		}

		// size of dimension 3 or y
		if ( i == 5 )
		{
			mxArray *aPtr = mxCreateDoubleMatrix(1,1,mxREAL);
			double *data = mxGetPr(aPtr);
			data[0] = static_cast<double>(size_y);
			mxSetField(s, 0, fieldnames[i], aPtr);
		}

		// size of dimension 4 or z
		if ( i == 6 )
		{
			mxArray *aPtr = mxCreateDoubleMatrix(1,1,mxREAL);
			double *data = mxGetPr(aPtr);
			data[0] = static_cast<double>(size_z);
			mxSetField(s, 0, fieldnames[i], aPtr);
		}
    }
	return s;

}

#endif

/**
* Write matrix to .mat file.
* Creates a struct for the matrix.
* Struct has 7 fields in including - arr time size, size1, size2, size3, size4
*/
template<class T>
void Matrix4D<T>::writeToMatlabFile(const std::string& filename, const std::string& info) const {

#if (MATLAB_CAPABLE)


	Logger::message << "writing " << filename << ": " << endl;

	// Open file for writing
	MATFile *pmat = matOpen( filename.c_str(), "w");
	if (pmat == NULL) {
   		printf("Error creating file %s\n", filename.c_str());
    	exit(EXIT_FAILURE);
  	}



	mwSize size_W = size_w;
	mwSize size_X = size_x;
	mwSize size_Y = size_y;
	mwSize size_Z = size_z;
	size_t w,x,y,z;
	int status;
    // create a struct array with two fields
    const char *fieldnames[7] = {"arr", "time", "size", "size1", "size2" , "size3" , "size4"};
	// create a 1x1 struct that will hold the array of values and the time info
	mxArray *s = mxCreateStructMatrix(1, 1, 7, fieldnames);

    // fill struct fields
    for (mwIndex i=0; i<7; i++) {
        // For the arr variable in the 1x1 struct
		if (i== 0)
		{
			// create array of doubles and give it the data in matrix_array
			const int ndim = 4;
			mwSize dims[ndim] = {size_W , size_X, size_Y, size_Z};
			mxArray *aPtr = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
			if (aPtr == NULL) {
    			printf("Unable to create double mxArray.\n");
    			exit(EXIT_FAILURE);
  			}
			double *data = mxGetPr(aPtr);
			//int bytes_to_copy = size_w * size_x * size_y * size_z * mxGetElementSize(aPtr);
			for (w = 0; w < size_w; w++) {
				for (x = 0; x < size_x; x++) {
					for (y = 0; y < size_y; y++) {
						for (z = 0; z < size_z; z++) {
							data[z*(size_w * size_x * size_y) + y*(size_w * size_x) +  x*(size_w) + w] = matrix_array[w][x][y][z];
						}
					}
				}
			}
			mxSetField(s, 0, fieldnames[i], aPtr);
		}
		// For the time variable in the 1x1 struct
		if ( i == 1 )
		{
			// Create a std::string that has the inputted time
			mxArray *aPtr = mxCreateDoubleMatrix(1,1,mxREAL);
			double *data = mxGetPr(aPtr);
			double temp = 0.0;
			std::stringstream ss;
			ss << info;
			ss >> temp;
			data[0] = temp;
			mxSetField(s, 0, fieldnames[i], aPtr);
		}
		// For the size variable in the 1x1 struct
		if ( i == 2 )
		{
			// Create a std::string that has the inputted time
			mxArray *aPtr = mxCreateDoubleMatrix(1,1,mxREAL);
			double *data = mxGetPr(aPtr);
			data[0] = 1;
			mxSetField(s, 0, fieldnames[i], aPtr);
		}
		if ( i == 3 )
		{
			// Create a std::string that has the inputted time
			mxArray *aPtr = mxCreateDoubleMatrix(1,1,mxREAL);
			double *data = mxGetPr(aPtr);
			data[0] = static_cast<double>(size_w);
			mxSetField(s, 0, fieldnames[i], aPtr);
		}
		if ( i == 4 )
		{
			// Create a std::string that has the inputted time
			mxArray *aPtr = mxCreateDoubleMatrix(1,1,mxREAL);
			double *data = mxGetPr(aPtr);
			data[0] = static_cast<double>(size_x);
			mxSetField(s, 0, fieldnames[i], aPtr);
		}
		if ( i == 5 )
		{
			// Create a std::string that has the inputted time
			mxArray *aPtr = mxCreateDoubleMatrix(1,1,mxREAL);
			double *data = mxGetPr(aPtr);
			data[0] = static_cast<double>(size_y);
			mxSetField(s, 0, fieldnames[i], aPtr);
		}
		if ( i == 6 )
		{
			// Create a std::string that has the inputted time
			mxArray *aPtr = mxCreateDoubleMatrix(1,1,mxREAL);
			double *data = mxGetPr(aPtr);
			data[0] = static_cast<double>(size_z);
			mxSetField(s, 0, fieldnames[i], aPtr);
		}
    }
	// Save the struct into the .mat file
	// hardcoded to get rid of the "./output_folder/" and only get PSD_XXXXX
    size_t str_idx = filename.rfind("PSD_");
    status = matPutVariable(pmat, filename.substr(str_idx,9).c_str(), s);
    if (status != 0) {
    	printf("error matputvariable");
    	exit(EXIT_FAILURE);
 	}


	// cleanup
    mxDestroyArray(s);
    matClose(pmat);

#endif

}




/**
* Write matrix to file, using 4 other matrixes as a grid (simply - write all 5 matrixes to the file).
* File has two header lines.
*/
template<class T>
void Matrix4D<T>::writeToFile(const std::string& filename, const Matrix4D<T> &grid_w, const Matrix4D<T> &grid_x, const Matrix4D<T> &grid_y, const Matrix4D<T> &grid_z) const {
	ofstream output(filename.c_str());
	//if (output==NULL && (filename.find("Debug") == std::string::npos)) {
	if (!output.is_open() && (filename.find("Debug") == std::string::npos)) {
		printf("FILE: Unable to output file: %s", filename.c_str());
		exit(EXIT_FAILURE);
	}
	output << "VARIABLES = \"" << ((grid_w.name!="")?grid_w.name:"w") << "\", \"" << ((grid_x.name!="")?grid_x.name:"x") << "\", \"" << ((grid_y.name!="")?grid_y.name:"y") << "\", \"" << ((grid_z.name!="")?grid_z.name:"z") << "\", \"" << ((this->name!="")?this->name:"f") << "\" "<< endl;
	output << "ZONE T=\"" << filename << "\", I=" << size_z << ", J=" << size_y << ", K=" << size_x << ", L=" << size_w << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (size_t w = 0; w < size_w; w++) {
		for (size_t x = 0; x < size_x; x++) {
			for (size_t y = 0; y < size_y; y++) {
				for (size_t z = 0; z < size_z; z++) {
					output << "\t" << grid_w[w][x][y][z] << "\t" << grid_x[w][x][y][z] << "\t" << grid_y[w][x][y][z] << "\t" << grid_z[w][x][y][z] << "\t" << matrix_array[w][x][y][z] << endl;
				}
			}
		}
	}
	output.close();
}



/**
* Write matrix to file, using 4 other matrixes as a grid (simply - write all 5 matrixes to the file).
* Uses the createStructMatrix() function to pack all the grid dimensions into seperate variables and then combines these variables into a single matlab structure to save in .mat
*/
template<class T>
void Matrix4D<T>::writeToMatlabFile(const std::string& file, const Matrix4D<T> &grid_w, const Matrix4D<T> &grid_x, const Matrix4D<T> &grid_y, const Matrix4D<T> &grid_z) const {

#if (MATLAB_CAPABLE)


	Logger::message << "writing " << file << ": " << endl;
	MATFile *pmat = matOpen( file.c_str(), "w");


	mxArray* current = createStructMatrix(file);
	mxArray* w = grid_w.createStructMatrix(file);
	mxArray* x = grid_x.createStructMatrix(file);
	mxArray* y = grid_y.createStructMatrix(file);
	mxArray* z = grid_z.createStructMatrix(file);

	// Save the struct into the .mat file
    matPutVariable(pmat, this->name.c_str(), current);
	matPutVariable(pmat, "grid_w", w);
	matPutVariable(pmat, "grid_x", x);
	matPutVariable(pmat, "grid_y", y);
	matPutVariable(pmat, "grid_z", z);
    mxDestroyArray(current);
	mxDestroyArray(w);
	mxDestroyArray(x);
	mxDestroyArray(y);
	mxDestroyArray(z);
    matClose(pmat);

#endif

}





/**
* Read matrix data from file with grid, by column
*
* Overloaded readFromFile function
* \param filename - file to read grids from
* \param read_column - read up to this column from file
*/
template<class T>
void Matrix4D<T>::readFromFile(const std::string& filename, int read_column) {
	std::string inBuf;
	int column;
	if (!initialized) {
		printf("MATRIX_ERROR: Using un-itialized matrix");
		exit(EXIT_FAILURE);
	} else {
		this->name = filename;
		ifstream input(filename.c_str());
		//if (input != NULL && !input.eof()) {
		if (input.is_open() && !input.eof()) {
			// Skipping first two lines.
			input >> inBuf;
			// !!!			while (strcmp(strupr((char *) inBuf.c_str()), "ZONE") != 0 && !input.eof() ) {
			while (strcasecmp(inBuf.c_str(), "ZONE") != 0 && !input.eof() ) {

				if (strcasecmp(inBuf.c_str(), "VARIABLES") != 0) {
					// Read variable names
					getline(input, inBuf);
					std::stringstream linestream(inBuf);
					std::string word;
					int ic = 0; // column count
					while (linestream >> word) {
						if (++ic == read_column) {
							if (word.length() >= 4)
								this->name = word.substr(1, word.size()-3);
							break;
						}
					}
				}

				input >> inBuf;
			}
			// read to the end of the line with 'zone'
			input.ignore(9999, '\n'); // read till the end of the line
			//getline(input, inBuf);  // read till the end of the line
			for (size_t w = 0; w < size_w; w++) {
				for (size_t x = 0; x < size_x; x++) {
					for (size_t y = 0; y < size_y; y++) {
						for (size_t z = 0; z < size_z; z++) {
							for (column = 1; column < read_column; column++) input >> inBuf; // skip first columns
							input >> matrix_array[w][x][y][z];
							//getline(input, inBuf);  // read till the end of the line
							input.ignore(9999, '\n'); // read till the end of the line
						}
					}
				}
			}
		} else {
			printf("MATRIX_LOAD_ERROR: Error reading file %s.\n", filename.c_str());
			exit(EXIT_FAILURE);
		}
		input >> inBuf;
		if (!input.eof()) {
			printf("Check grid size - %s.\n", filename.c_str());
			exit(EXIT_FAILURE);
		}
		input.close();
	}
}


/**
* Function for reading from matlab file in 4-dimensions
* Will check the variables, order them in (P, R/L, V, K, Val) format and then set matrix_array to be the variable with the corresponding column number
* This is the same as the readFromFile() function although only compatible with .mat files instead of .plt or other text files
*/
template<class T>
void Matrix4D<T>::readFromMatlabFile(const std::string& file , int columnNumber)
{

#if (MATLAB_CAPABLE)

	Logger::message << "Reading " << file << ": " << endl;

	MATFile *mfPtr; /* MAT-file pointer */
    mxArray *aPtr;  /* mxArray pointer */
    double *realPtr = nullptr; /* pointer to data */
	double *PtrW = nullptr; /* pointer to data */
	double *PtrX = nullptr; /* pointer to data */
    double *PtrY = nullptr; /* pointer to data */
    double *PtrZ = nullptr; /* pointer to data */
	double *PtrFinal = nullptr; /* pointer to data */
	double *PtrReturn = nullptr; /* pointer to data */
	double *PtrL = nullptr; /* pointer to data */
	std::string arr; /*name of variable*/
	std::string field = "arr"; // name of field
	mwSize nElements;       /* number of elements in array */
    mwIndex eIdx;           /* element index */
    const mxArray *fPtr;    /* field pointer */
	const mxArray *fPtrW;    /* field pointer */
	const mxArray *fPtrL;    /* field pointer */
	const mxArray *fPtrX;    /* field pointer */
	const mxArray *fPtrY;    /* field pointer */
	const mxArray *fPtrZ;    /* field pointer */
    size_t w,x,y,z; 			/* for index*/
	const char* name;		/* for getting variable names */
	char nameTemp;	/* for getting variable names */
	bool wReached = false;	/* for checking which variables go to which coordinates */
	bool xReached = false;	/* for checking which variables go to which coordinates */
	bool yReached = false;	/* for checking which variables go to which coordinates */
	bool zReached = false;	/* for checking which variables go to which coordinates */
	bool defaultReached = false;	/* for checking which variables go to which coordinates */
	bool rReached = false;	/* for checking which variables go to which coordinates */

	// open matlab file
	mfPtr = matOpen(file.c_str(), "r");
   	if (mfPtr == NULL) {
       	printf("Error opening file %s\n", file.c_str());
	}

	// for each of the 5 columns
	for (int i= 0; i < 5; i++)
	{
		// get variable name
		aPtr = matGetNextVariableInfo(mfPtr, &name);
		std::string temp = name;
		// If the name is one char long - match on that char
		if (temp.length() == 1)
	    {
	        nameTemp = temp[0];
	    }
		else
		{
			// if the second char is underscore '_' then match on first char example V_BC
			if (temp.substr(1, 1) == "_")
			{
				nameTemp = temp[0];
			}
			else if ((temp.substr(0, 2) == "pc"))
			{
				nameTemp = 'V';
			}
			else if ((temp.substr(0, 2) == "al"))
			{
				nameTemp = 'K';
			}
			// If length is more than 1 and second char is not underscore will go to default
			else
			{
				nameTemp = '!';
			}
		}

		// switch on name to get P R(or L) V K Value in that order
		// The code is nearly the same for all switch statements - comments are found only for case 'P'
		switch (nameTemp)
		{
			// If variable is the array for 'P'
			case 'P':
				// coordinate w has been recorded
				wReached = true;
				// get variable
				aPtr = matGetVariable(mfPtr, name);
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				// if variable is matlab struct
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					// make sure it has a "arr" field
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						// get the number of elements in struct
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						// for each element in struct
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							// try to get the field "arr"
							fPtrW = mxGetField(aPtr, eIdx, field.c_str());
							// if successfully got the field "arr" and it is an array of doubles
       		 				if ((fPtrW != NULL) && (mxGetClassID(fPtrW) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrW)))
       		 				{
								// get the data from "arr" (array of doubles)
								PtrW = mxGetPr(fPtrW);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'R':
				rReached = true;
				xReached = true;
				aPtr = matGetVariable(mfPtr, name);
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrX = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrX != NULL) && (mxGetClassID(fPtrX) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrX)))
       		 				{
								PtrX = mxGetPr(fPtrX);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'L':
				xReached = true;
				aPtr = matGetVariable(mfPtr, name);
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrL = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrL != NULL) && (mxGetClassID(fPtrL) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrL)))
       		 				{
								PtrL = mxGetPr(fPtrL);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'V':
				yReached = true;
				aPtr = matGetVariable(mfPtr, name);
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrY = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrY != NULL) && (mxGetClassID(fPtrY) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrY)))
       		 				{
								PtrY = mxGetPr(fPtrY);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'K':
				zReached = true;
				aPtr = matGetVariable(mfPtr, name);
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrZ = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrZ != NULL) && (mxGetClassID(fPtrZ) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrZ)))
       		 				{
								PtrZ = mxGetPr(fPtrZ);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			default:
				defaultReached = true;
				aPtr = matGetVariable(mfPtr, name);
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtr = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtr != NULL) && (mxGetClassID(fPtr) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtr)))
       		 				{
								PtrFinal = mxGetPr(fPtr);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
		}
	}



	// if no R but there is L
	if (!rReached)
	{
		if (xReached)
		{
			PtrX = PtrL;
		}
	}
	// if no default - must have been L star (5th parameter = L which should be Val in "P R V K Val")
	else if (!defaultReached)
	{
		PtrFinal = PtrL;
	}



	// To get the column number of the return
	// since it is assumed the variables are in P - R/L - V - K - Val form iterate through until we reach the correct column
	// Not necesarily required for the 4D case where all 5 variables will be present - however in 3D and below might be missing arbitrary V or P etc.
	if (wReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrW;
		}
	}
	if (xReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrX;
		}
	}
	if (yReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrY;
		}
	}
	if (zReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrZ;
		}
	}
	if (defaultReached)
	{
		columnNumber--;
		if (columnNumber == 0)
		{
			PtrReturn = PtrFinal;
		}
	}
	if (columnNumber > 0)
	{
		printf("column number too high \n");
		EXIT_FAILURE;
	}


	// sets the matrix array to be equal to an array of doubles
	for (w = 0; w < size_w; w++) {
		for (x = 0; x < size_x; x++) {
			for (y = 0; y < size_y; y++) {
				for (z = 0; z < size_z; z++) {
					// Saving variable from matlab requires reversing indeces because matlab is stored in column major order while c++ is row major order
					matrix_array[w][x][y][z] = PtrReturn[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w];
				}
			}
		}
	}

	// close file
    if (matClose(mfPtr) != 0) {
        printf("Error closing file %s\n", file.c_str());
    }

	// free allocated memory
	mxDestroyArray(aPtr);

#endif

}




















/**
* Read matrix data from file with grid,
* Checks if the matrix data in the file is the same as the grids that were sent in with error < 1e-8,
* if not within error range will signal error and exit
*
* Overloaded readFromFile function
* \param filename - file to read grids from
* \param grids w,x,y,z - checks grids data against the file data
*/
template<class T>
void Matrix4D<T>::readFromFile(const std::string& filename, const Matrix4D<T>& grid_w, const Matrix4D<T>& grid_x, const Matrix4D<T>& grid_y, const Matrix4D<T>& grid_z) {
	std::string inBuf;
	double loaded_w, loaded_x, loaded_y, loaded_z;

	if (!initialized) {
		printf("MATRIX_ERROR: Using uninitialized matrix");
		exit(EXIT_FAILURE);
	} else {
		this->name = filename;
		ifstream input(filename.c_str());
		//if (input != NULL && !input.eof()) {
		if (input.is_open() && !input.eof()) {
			// Skipping first two lines.
			input >> inBuf;
			// !!!			while (strcmp(_strupr((char *) inBuf.c_str()), "ZONE") != 0 && !input.eof() ) {
			// !!!			while (strcasecmp(inBuf.c_str(), "ZONE") != 0 && !input.eof() ) {
			while (strcmp(inBuf.c_str(), "ZONE") != 0 && !input.eof() ) {

				if (strcasecmp(inBuf.c_str(), "VARIABLES") != 0) {
					// Read variable names
					getline(input, inBuf);
					std::stringstream linestream(inBuf);
					std::string word;
					int ic = 0; // column count
					int read_column = 5;
					while (linestream >> word) {
						if (++ic == read_column) {
							if (word.length() >= 4)
								this->name = word.substr(1, word.size()-3);
							break;
						}
					}
				}

				input >> inBuf;
			}
			// read to the end of the line with 'zone'
			input.ignore(9999, '\n');
			for (size_t w = 0; w < size_w; w++) {
				for (size_t x = 0; x < size_x; x++) {
					for (size_t y = 0; y < size_y; y++) {
						for (size_t z = 0; z < size_z; z++) {
							//input >> grid_x[x][y][z] >> grid_y[x][y][z] >> grid_z[x][y][z] >> matrix_array[x][y][z];
							// skip till the end of the line
							//input.ignore(9999, '\n');

							 input >> loaded_w >> loaded_x >> loaded_y >> loaded_z;
							 // check if grid is the same
							 if (fabs(log10(loaded_w) - log10(grid_w[w][x][y][z])) > err || fabs(log10(loaded_x) - log10(grid_x[w][x][y][z])) > err || fabs(log10(loaded_y) - log10(grid_y[w][x][y][z])) > err || fabs(log10(loaded_z) - log10(grid_z[w][x][y][z])) > err) {
								printf("MATRIX_LOAD_GRID_ERR: Loading %s: grid mismatch [%zu, %zu, %zu, %zu].\nLoaded: %e, %e, %e, %e\nGrid: %e, %e, %e, %e\n", filename.c_str(), w, x, y, z, loaded_w, loaded_x, loaded_y, loaded_z, grid_w[w][x][y][z], grid_x[w][x][y][z], grid_y[w][x][y][z], grid_z[w][x][y][z]);
								exit(EXIT_FAILURE);
							} else {
								input >> matrix_array[w][x][y][z];
								//ADDED
								//printf( "%e \n",  matrix_array[w][x][y][z]);
								//std::cout <<  input.rdbuf();
								// skip till to the end of the line
								input.ignore(9999, '\n');
							}
						}
					}
				}
			}
		} else {
			printf("MATRIX_LOAD_ERROR: Error reading file %s.\n", filename.c_str());
			exit(EXIT_FAILURE);
		}
		input.close();
	}
}




/**
* Function for reading from matlab file in 4-dimensions
* Will check the variables in the order they are saved in matlab, thus (P R V K Var) should be the standard
* The variables will be checked against the input grid parameters in order to make sure the right variables/values are being loaded
* This is the same as the readFromFile() function although only compatible with .mat files instead of .plt or other text files
*/
template<class T>
void Matrix4D<T>::readFromMatlabFile(const std::string& file , const Matrix4D<T>& grid_w, const Matrix4D<T>& grid_x, const Matrix4D<T>& grid_y, const Matrix4D<T>& grid_z)
{

#if (MATLAB_CAPABLE)

	Logger::message << "Reading " << file << ": " << endl;

	MATFile *mfPtr; /* MAT-file pointer */
    mxArray *aPtr;  /* mxArray pointer */
    double *realPtr = nullptr; /* pointer to data */
	double *PtrW = nullptr; /* pointer to data */
	double *PtrX = nullptr; /* pointer to data */
    double *PtrY = nullptr; /* pointer to data */
    double *PtrZ = nullptr; /* pointer to data */
	double *PtrFinal = nullptr; /* pointer to data */
	double *PtrL = nullptr; /* pointer to data */
	std::string arr; /*name of variable*/
	std::string field = "arr"; // name of field
	mwSize nElements;       /* number of elements in array */
    mwIndex eIdx;           /* element index */
    const mxArray *fPtr;    /* field pointer */
	const mxArray *fPtrW;    /* field pointer */
	const mxArray *fPtrL;    /* field pointer */
	const mxArray *fPtrX;    /* field pointer */
	const mxArray *fPtrY;    /* field pointer */
	const mxArray *fPtrZ;    /* field pointer */
    size_t w,x,y,z; 			/* for index*/
	const char* name;		/* for getting variable names */
	char nameTemp;	/* for getting variable names */
	bool defaultReached = false;
	bool rReached = false;
	bool lReached = false;


	mfPtr = matOpen(file.c_str(), "r");
   	if (mfPtr == NULL) {
       	printf("Error opening file %s\n", file.c_str());
	}

	for (int i= 0; i < 5; i++)
	{
		aPtr = matGetNextVariableInfo(mfPtr, &name);
		std::string temp = name;
		if (temp.length() == 1)
	    {
	        nameTemp = temp[0];
	    }
		else
		{
			if (temp.substr(1, 1) == "_")
			{
				nameTemp = temp[0];
			}
			else if ((temp.substr(0, 2) == "pc"))
			{
				nameTemp = 'V';
			}
			else if ((temp.substr(0, 2) == "al"))
			{
				nameTemp = 'K';
			}
			else
			{
				nameTemp = '!';
			}
		}
		switch (nameTemp)
		{
			case 'P':
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrW = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrW != NULL) && (mxGetClassID(fPtrW) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrW)))
       		 				{
								PtrW = mxGetPr(fPtrW);
								//printf("%e \n", PtrW[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'R':
				rReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrX = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrX != NULL) && (mxGetClassID(fPtrX) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrX)))
       		 				{
								PtrX = mxGetPr(fPtrX);
								//printf("%e \n", PtrX[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'L':
				lReached = true;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrL = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrL != NULL) && (mxGetClassID(fPtrL) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrL)))
       		 				{
								PtrL = mxGetPr(fPtrL);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'V':
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrY = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrY != NULL) && (mxGetClassID(fPtrY) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrY)))
       		 				{
								PtrY = mxGetPr(fPtrY);
								//printf("%e \n", PtrY[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			case 'K':
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtrZ = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtrZ != NULL) && (mxGetClassID(fPtrZ) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrZ)))
       		 				{
								PtrZ = mxGetPr(fPtrZ);
								//printf("%e \n", PtrZ[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
			default:
				defaultReached = true;
				this->name = name;
				aPtr = matGetVariable(mfPtr, name);
				//printf("%s \n" , name );
				if (aPtr == NULL) {
  					printf("mxArray not found: %s\n", name);
				}
				if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
					if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
   						printf("Field not found: %s\n", field.c_str());
	 				}
					else {
						nElements = (mwSize)mxGetNumberOfElements(aPtr);
						for (eIdx = 0; eIdx < nElements; eIdx++) {
							fPtr = mxGetField(aPtr, eIdx, field.c_str());
       		 				if ((fPtr != NULL) && (mxGetClassID(fPtr) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtr)))
       		 				{
								PtrFinal = mxGetPr(fPtr);
								//printf("%e \n", PtrFinal[0]);
							}
						}
					}
				}
				else
				{
					printf("%s is of unknown type\n", name);
			    }
				break;
		}
	}

	if (!rReached && lReached)
	{
		PtrX = PtrL;
	}
	else if (!defaultReached && lReached)
	{
		PtrFinal = PtrL;
	}

	// sets the matrix array to be equal to an array of doubles
	for (w = 0; w < size_w; w++) {
		for (x = 0; x < size_x; x++) {
			for (y = 0; y < size_y; y++) {
				for (z = 0; z < size_z; z++) {


					if (fabs(log10(PtrW[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w]) - log10(grid_w[w][x][y][z])) > err || fabs(log10(PtrX[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w]) - log10(grid_x[w][x][y][z])) > err || fabs(log10(PtrY[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w]) - log10(grid_y[w][x][y][z])) > err || fabs(log10(PtrZ[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w]) - log10(grid_z[w][x][y][z])) > err) {
					 			printf("MATRIX_LOAD_GRID_ERR: Loading %s: grid mismatch [%zu, %zu, %zu, %zu].\nLoaded: %e, %e, %e, %e\nGrid: %e, %e, %e, %e\n", file.c_str(), w, x, y, z, PtrW[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w], PtrX[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w],PtrY[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w],PtrZ[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w], grid_w[w][x][y][z], grid_x[w][x][y][z], grid_y[w][x][y][z], grid_z[w][x][y][z]);
								exit(EXIT_FAILURE);
					}
					matrix_array[w][x][y][z] = PtrFinal[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w];
				}
			}
		}
	}

    if (matClose(mfPtr) != 0) {
        printf("Error closing file %s\n", file.c_str());
    }

	mxDestroyArray(aPtr);

#endif

}




/**
* Write matrix data to binary file

* \param filename - file to read grids from
*/
template<class T>
void Matrix4D<T>::writeToBinaryFile(const std::string& filename) const {
	std::ofstream outputFile(filename, std::ios::binary | std::ios::out);
    if (outputFile.is_open()) {
		int32_t size_array[4] = { (int32_t)size_w, (int32_t)size_x, (int32_t)size_y, (int32_t)size_z };
		outputFile.write((char*)size_array, 4 * sizeof(int32_t));
		if (!outputFile.good()){
			printf("Writing error");
			exit(EXIT_FAILURE);
		}

		outputFile.write((char*)plane_array, num_elements * sizeof(T));
		if (!outputFile.good()){
			printf("Writing error");
			exit(EXIT_FAILURE);
		}
    } else {
        printf("MATRIX_WRITE_ERROR: Error writing file %s.\n", filename.c_str());
        exit(EXIT_FAILURE);
    }
    outputFile.close();
}

/**
* Write 4D data to .plt, .pltb or .mat files. No info is written for .pltb file
*/
template<class T>
void Matrix4D<T>::writeToAnyFile(const std::string& filename, const std::string& io_method, const std::string& info) const {
    std::string ext;
    if (io_method.compare("ascii") == 0){
        ext = ".plt";
        writeToFile(filename + ext, info);
    } else if (io_method.compare("binary") == 0) {
        ext = ".pltb";
        writeToBinaryFile(filename + ext);
    } else if (io_method.compare("matlab") == 0) {
        ext = ".mat";
        writeToMatlabFile(filename + ext, info);
    } else {
        printf("IO error: unknown io_method");
        exit(EXIT_FAILURE);
    }
}


/**
* Read matrix data from binary file
*
* \param filename - file to read grids from
*/
template<class T>
void Matrix4D<T>::readFromBinaryFile(const std::string& filename) {
	if (!initialized) {
		printf("MATRIX_ERROR: Using un-itialized matrix");
		exit(EXIT_FAILURE);
	} else {
		this->name = filename;
		std::ifstream inputFile(filename, std::ios::binary | std::ios::in);
		if (inputFile.is_open()) {
			int32_t buffer_int32[4];
			inputFile.read((char*)buffer_int32, 4 * sizeof(int32_t));
			if (!inputFile.good()){
				printf("Reading error");
				exit(EXIT_FAILURE);
			}

			inputFile.read((char*)plane_array, num_elements * sizeof(T));
			if (!inputFile.good()){
				printf("Reading error");
				exit(EXIT_FAILURE);
			}
		} else {
			printf("MATRIX_LOAD_ERROR: Error reading file %s.\n", filename.c_str());
			exit(EXIT_FAILURE);
		}
		inputFile.close();
	}
}


/**
* Read 4D data from .plt, .pltb or .mat files containing only one column
*/
template<class T>
void Matrix4D<T>::readFromAnyFile(const std::string& filename, const std::string& io_method) {
    std::string ext;
    if (io_method.compare("ascii") == 0){
        ext = ".plt";
        readFromFile(filename + ext, 1);
    } else if (io_method.compare("binary") == 0){
        ext = ".pltb";
        readFromBinaryFile(filename + ext);
    } else if (io_method.compare("matlab") == 0){
        ext = ".mat";
        readFromMatlabFile(filename + ext, 1);
    } else {
        printf("IO error: unknown io_method");
        exit(EXIT_FAILURE);
    }
}

/**
* Read 4D data from .plt, .pltb or .mat files with grids
* WARNING: if io_method == "binary", .pltb file does not contain grid. In this case the function works the same as Matrix4D<T>::readFromAnyFile(std::string, std::string, int)
*/
template<class T>
void Matrix4D<T>::readFromAnyFile(const std::string& filename, const std::string& io_method, const Matrix4D<T>& grid_w, const Matrix4D<T>& grid_x, const Matrix4D<T>& grid_y, const Matrix4D<T>& grid_z){
    std::string ext;
    if (io_method.compare("ascii") == 0){
        ext = ".plt";
        readFromFile(filename + ext, grid_w, grid_x, grid_y, grid_z);
    } else if (io_method.compare("binary") == 0) {
        ext = ".pltb";
        readFromBinaryFile(filename + ext);
    } else if (io_method.compare("matlab") == 0) {
        ext = ".mat";
        readFromMatlabFile(filename + ext, grid_w, grid_x, grid_y, grid_z);
    } else {
        printf("IO error: unknown io_method");
        exit(EXIT_FAILURE);
    }
}


/**
* Returns corresponding index of 4d matrix if represented as a 1d array
* \param w,x,y,z - index of element in every dimension for the 4d matrix
*/
template<class T>
inline size_t Matrix4D<T>::index1d(size_t w, size_t x, size_t y, size_t z) {
	return ((w*size_x + x)*size_y + y)*size_z;
}



/**
* Return minimum value of the 4d matrix.
* Default value set to 1e99
*/
template<class T>
T Matrix4D<T>::min() const {
	T tmp = std::numeric_limits<T>::infinity();
	for(size_t i = 0; i < num_elements; i++)
	{
		tmp = std::min(tmp, plane_array[i]);
	}
	return tmp;
}

/**
* Return maximum value of the 4d matrix.
* Default value set to 0
*/
template<class T>
T Matrix4D<T>::max() const {
	T tmp = -std::numeric_limits<T>::infinity();
	for(size_t i = 0; i < num_elements; i++)
	{
		tmp = std::max(tmp, plane_array[i]);
	}
	return tmp;
}

/**
* Return absolute maximum value of the 4d matrix.
* Default value set to 0
*/
template<class T>
T Matrix4D<T>::maxabs() const {
	T tmp = 0;
	for(size_t i = 0; i < num_elements; i++)
	{
		tmp = std::max(tmp, std::abs(plane_array[i]));
	}
	return tmp;
}

/**
* Return absolute value of the 4d matrix.
* Changes every element to a positive value with the same magnitude
*/
template<class T>
Matrix4D<T> Matrix4D<T>::abs() const {
	Matrix4D<T> tmp(size_w, size_x, size_y, size_z);
	for(size_t i = 0; i < num_elements; i++)
	{
		tmp.plane_array[i] = std::abs(plane_array[i]);
	}
	return tmp;
}

/**
* Replace all matrix malues smaller than val with val.
*
* \param val - new minimal value.
*/
template<class T>
Matrix4D<T>& Matrix4D<T>::max_of(T val) {
	for(size_t i = 0; i < num_elements; i++)
	{
		plane_array[i] = std::max(plane_array[i], val);
	}
	return *this;
}

/**
* Take w-slice of 4d matrix turning it into 3d matrix.
* squeeze matrix to 3d getting all the data when the w dimension is at index p_w
* \param p_w - index at which to slice w dimension
*/
template<class T>
Matrix3D<T> Matrix4D<T>::wSlice(size_t p_w) const {
	Matrix3D<T> tmp(this->size_x, this->size_y, this->size_z);
	// !!! tmp.name = this->name + "_slice";
	for (size_t x = 0; x < size_x; x++) {
		for (size_t y = 0; y < size_y; y++) {
			for (size_t z = 0; z < size_z; z++) {
				tmp[x][y][z] = matrix_array[p_w][x][y][z];
			}
		}
	}
	return tmp;
}

/**
* Take a w-slice of 4d matrix and write it to out.
* squeeze matrix to 3d getting all the data when the w dimension is at index p_w
* \param p_w - index at which to slice w dimension
*/
template<class T>
void Matrix4D<T>::wSlice(Matrix3D<T>& out, size_t p_w) const {
	for (size_t x = 0; x < size_x; x++) {
		for (size_t y = 0; y < size_y; y++) {
			for (size_t z = 0; z < size_z; z++) {
				out[x][y][z] = matrix_array[p_w][x][y][z];
			}
		}
	}
}

/**
* Take x-slice of 4d matrix turining it into 3d matrix.
* squeeze matrix to 3d getting all the data when the x dimension is at index p_x
* \param p_x - index at which to slice x dimension
*/
template<class T>
Matrix3D<T> Matrix4D<T>::xSlice(size_t p_x) const {
	Matrix3D<T> tmp(this->size_w, this->size_y, this->size_z);
	// !!! tmp.name = this->name + "_slice";
	for (size_t w = 0; w < size_w; w++) {
		for (size_t y = 0; y < size_y; y++) {
			for (size_t z = 0; z < size_z; z++) {
				tmp[w][y][z] = matrix_array[w][p_x][y][z];
			}
		}
	}
	return tmp;
}

/**
* Take an x-slice of 4d matrix and write it to out.
* squeeze matrix to 3d getting all the data when the x dimension is at index p_x
* \param p_x - index at which to slice x dimension
*/
template<class T>
void Matrix4D<T>::xSlice(Matrix3D<T>& out, size_t p_x) const {
	for (size_t w = 0; w < size_w; w++) {
		for (size_t y = 0; y < size_y; y++) {
			for (size_t z = 0; z < size_z; z++) {
				out[w][y][z] = matrix_array[w][p_x][y][z];
			}
		}
	}
}

/**
* Take y-slice of 4d matrix turning it into 3d matrix.
* squeeze matrix to 3d getting all the data when the y dimension is at index p_y
* \param p_y - index at which to slice y dimension
*/
template<class T>
Matrix3D<T> Matrix4D<T>::ySlice(size_t p_y) const {
	Matrix3D<T> tmp(this->size_w, this->size_x, this->size_z);
	// !!! tmp.name = this->name + "_slice";
	for (size_t w = 0; w < size_w; w++) {
		for (size_t x = 0; x < size_x; x++) {
			for (size_t z = 0; z < size_z; z++) {
				tmp[w][x][z] = matrix_array[w][x][p_y][z];
			}
		}
	}
	return tmp;
}

/**
* Take a y-slice of 4d matrix and write it to out.
* squeeze matrix to 3d getting all the data when the y dimension is at index p_y
* \param p_y - index at which to slice y dimension
*/
template<class T>
void Matrix4D<T>::ySlice(Matrix3D<T>& out, size_t p_y) const {
	for (size_t w = 0; w < size_w; w++) {
		for (size_t x = 0; x < size_x; x++) {
			for (size_t z = 0; z < size_z; z++) {
				out[w][x][z] = matrix_array[w][x][p_y][z];
			}
		}
	}
}

/**
* Take z-slice of 4d matrix turning it into 3d matrix.
* squeeze matrix to 3d getting all the data when the z dimension is at index p_z
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix3D<T> Matrix4D<T>::zSlice(size_t p_z) const {
	Matrix3D<T> tmp(this->size_w, this->size_x, this->size_y);
	// !!! tmp.name = this->name + "_slice";
	for (size_t w = 0; w < size_w; w++) {
		for (size_t x = 0; x < size_x; x++) {
			for (size_t y = 0; y < size_y; y++) {
				tmp[w][x][y] = matrix_array[w][x][y][p_z];
			}
		}
	}
	return tmp;
}

/**
* Take a z-slice of 4d matrix and write it to out.
* squeeze matrix to 3d getting all the data when the z dimension is at index p_z
* \param p_z - index at which to slice z dimension
*/
template<class T>
void Matrix4D<T>::zSlice(Matrix3D<T>& out, size_t p_z) const {
	for (size_t w = 0; w < size_w; w++) {
		for (size_t x = 0; x < size_x; x++) {
			for (size_t y = 0; y < size_y; y++) {
				out[w][x][y] = matrix_array[w][x][y][p_z];
			}
		}
	}
}

/**
* Take wx-slice of 4d matrix turning it into a 2d matrix.
* \param p_w - index at which to slice w dimension
* \param p_x - index at which to slice x dimension
*/
template<class T>
Matrix2D<T> Matrix4D<T>::wxSlice(size_t p_w, size_t p_x) const {
	Matrix2D<T> tmp(this->size_y, this->size_z);
	// !!! tmp.name = this->name + "_xySlice";
	for (size_t y = 0; y < size_y; y++) {
		for (size_t z = 0; z < size_z; z++) {
			tmp[y][z] = matrix_array[p_w][p_x][y][z];
		}
	}
	return tmp;
}

/**
* Take a wx-slice of 4d matrix and write it to out.
* \param p_w - index at which to slice w dimension
* \param p_x - index at which to slice x dimension
*/
template<class T>
void Matrix4D<T>::wxSlice(Matrix2D<T>& out, size_t p_w, size_t p_x) const {
	for (size_t y = 0; y < size_y; y++) {
		for (size_t z = 0; z < size_z; z++) {
			out[y][z] = matrix_array[p_w][p_x][y][z];
		}
	}
}

/**
* Take wy-slice of 4d matrix turning it into a 2d matrix.
* \param p_w - index at which to slice w dimension
* \param p_y - index at which to slice y dimension
*/
template<class T>
Matrix2D<T> Matrix4D<T>::wySlice(size_t p_w, size_t p_y) const {
	Matrix2D<T> tmp(this->size_x, this->size_z);
	// !!! tmp.name = this->name + "_xySlice";
	for (size_t x = 0; x < size_x; x++) {
		for (size_t z = 0; z < size_z; z++) {
			tmp[x][z] = matrix_array[p_w][x][p_y][z];
		}
	}
	return tmp;
}

/**
* Take a wy-slice of 4d matrix and write it to out.
* \param p_w - index at which to slice w dimension
* \param p_y - index at which to slice y dimension
*/
template<class T>
void Matrix4D<T>::wySlice(Matrix2D<T>& out, size_t p_w, size_t p_y) const {
	for (size_t x = 0; x < size_x; x++) {
		for (size_t z = 0; z < size_z; z++) {
			out[x][z] = matrix_array[p_w][x][p_y][z];
		}
	}
}

/**
* Take wz-slice of 4d matrix turning it into a 2d matrix.
* \param p_w - index at which to slice w dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix2D<T> Matrix4D<T>::wzSlice(size_t p_w, size_t p_z) const {
	Matrix2D<T> tmp(this->size_x, this->size_y);
	// !!! tmp.name = this->name + "_xySlice";
	for (size_t x = 0; x < size_x; x++) {
		for (size_t y = 0; y < size_y; y++) {
			tmp[x][y] = matrix_array[p_w][x][y][p_z];
		}
	}
	return tmp;
}

/**
* Take  wz-slice of 4d matrix and write it to out.
* \param p_w - index at which to slice w dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
void Matrix4D<T>::wzSlice(Matrix2D<T>& out, size_t p_w, size_t p_z) const {
	for (size_t x = 0; x < size_x; x++) {
		for (size_t y = 0; y < size_y; y++) {
			out[x][y] = matrix_array[p_w][x][y][p_z];
		}
	}
}

/**
* Take xy-slice of 4d matrix turning it into a 2d matrix.
* \param p_x - index at which to slice x dimension
* \param p_y - index at which to slice y dimension
*/
template<class T>
Matrix2D<T> Matrix4D<T>::xySlice(size_t p_x, size_t p_y) const {
	Matrix2D<T> tmp(this->size_w, this->size_z);
	// !!! tmp.name = this->name + "_xySlice";
	for (size_t w = 0; w < size_w; w++) {
		for (size_t z = 0; z < size_z; z++) {
			tmp[w][z] = matrix_array[w][p_x][p_y][z];
		}
	}
	return tmp;
}

/**
* Take a xy-slice of 4d matrix and write it to out.
* \param p_x - index at which to slice x dimension
* \param p_y - index at which to slice y dimension
*/
template<class T>
void Matrix4D<T>::xySlice(Matrix2D<T>& out, size_t p_x, size_t p_y) const {
	for (size_t w = 0; w < size_w; w++) {
		for (size_t z = 0; z < size_z; z++) {
			out[w][z] = matrix_array[w][p_x][p_y][z];
		}
	}
}

/**
* Take xz-slice of 4d matrix turning it into a 2d matrix.
* \param p_x - index at which to slice x dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix2D<T> Matrix4D<T>::xzSlice(size_t p_x, size_t p_z) const {
	Matrix2D<T> tmp(this->size_w, this->size_y);
	// !!! tmp.name = this->name + "_xySlice";
	for (size_t w = 0; w < size_w; w++) {
		for (size_t y = 0; y < size_y; y++) {
			tmp[w][y] = matrix_array[w][p_x][y][p_z];
		}
	}
	return tmp;
}

/**
* Take an xz-slice of 4d matrix and write it to out.
* \param p_x - index at which to slice x dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
void Matrix4D<T>::xzSlice(Matrix2D<T>& out, size_t p_x, size_t p_z) const {
	for (size_t w = 0; w < size_w; w++) {
		for (size_t y = 0; y < size_y; y++) {
			out[w][y] = matrix_array[w][p_x][y][p_z];
		}
	}
}

/**
* Take yz-slice of 4d matrix turning it into a 2d matrix.
* \param p_y - index at which to slice y dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix2D<T> Matrix4D<T>::yzSlice(size_t p_y, size_t p_z) const {
	Matrix2D<T> tmp(this->size_w, this->size_x);
	// !!! tmp.name = this->name + "_xySlice";
	for (size_t w = 0; w < size_w; w++) {
		for (size_t x = 0; x < size_x; x++) {
			tmp[w][x] = matrix_array[w][x][p_y][p_z];
		}
	}
	return tmp;
}

/**
* Take a yz-slice of 4d matrix and write it to out.
* \param p_y - index at which to slice y dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
void Matrix4D<T>::yzSlice(Matrix2D<T>& out, size_t p_y, size_t p_z) const {
	for (size_t w = 0; w < size_w; w++) {
		for (size_t x = 0; x < size_x; x++) {
			out[w][x] = matrix_array[w][x][p_y][p_z];
		}
	}
}

/**
* Take wxy-slice of 4d matrix turning it into a 1d matrix.
* \param p_w - index at which to slice w dimension
* \param p_x - index at which to slice x dimension
* \param p_y - index at which to slice y dimension
*/
template<class T>
Matrix1D<T> Matrix4D<T>::wxySlice(size_t p_w, size_t p_x, size_t p_y) const {
	Matrix1D<T> tmp(this->size_z);
	// !!! tmp.name = this->name + "_xySlice";
	for (size_t z = 0; z < size_z; z++) {
		tmp[z] = matrix_array[p_w][p_x][p_y][z];
	}
	return tmp;
}

/**
* Take a wxy-slice of 4d matrix and write it to out.
* \param p_w - index at which to slice w dimension
* \param p_x - index at which to slice x dimension
* \param p_y - index at which to slice y dimension
*/
template<class T>
void Matrix4D<T>::wxySlice(Matrix1D<T>& out, size_t p_w, size_t p_x, size_t p_y) const {
	for (size_t z = 0; z < size_z; z++) {
		out[z] = matrix_array[p_w][p_x][p_y][z];
	}
}

/**
* Take wxz-slice of 4d matrix turning it into a 1d matrix.
* \param p_w - index at which to slice w dimension
* \param p_x - index at which to slice x dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix1D<T> Matrix4D<T>::wxzSlice(size_t p_w, size_t p_x, size_t p_z) const {
	Matrix1D<T> tmp(this->size_y);
	// !!! tmp.name = this->name + "_xySlice";
	for (size_t y = 0; y < size_y; y++) {
		tmp[y] = matrix_array[p_w][p_x][y][p_z];
	}
	return tmp;
}

/**
* Take a wxz-slice of 4d matrix and write it to out.
* \param p_w - index at which to slice w dimension
* \param p_x - index at which to slice x dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
void Matrix4D<T>::wxzSlice(Matrix1D<T>& out, size_t p_w, size_t p_x, size_t p_z) const {
	for (size_t y = 0; y < size_y; y++) {
		out[y] = matrix_array[p_w][p_x][y][p_z];
	}
}

/**
* Take wyz-slice of 4d matrix turning it into a 1d matrix.
* \param p_w - index at which to slice w dimension
* \param p_y - index at which to slice y dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix1D<T> Matrix4D<T>::wyzSlice(size_t p_w, size_t p_y, size_t p_z) const {
	Matrix1D<T> tmp(this->size_x);
	// !!! tmp.name = this->name + "_xySlice";
	for (size_t x = 0; x < size_x; x++) {
		tmp[x] = matrix_array[p_w][x][p_y][p_z];
	}
	return tmp;
}

/**
* Take a wyz-slice of 4d matrix and write it to out.
* \param p_w - index at which to slice w dimension
* \param p_y - index at which to slice y dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
void Matrix4D<T>::wyzSlice(Matrix1D<T>& out, size_t p_w, size_t p_y, size_t p_z) const {
	for (size_t x = 0; x < size_x; x++) {
		out[x] = matrix_array[p_w][x][p_y][p_z];
	}
}

/**
* Take xyz-slice of 4d matrix turning it into a 1d matrix.
* \param p_x - index at which to slice x dimension
* \param p_y - index at which to slice y dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix1D<T> Matrix4D<T>::xyzSlice(size_t p_x, size_t p_y, size_t p_z) const {
	Matrix1D<T> tmp(this->size_w);
	// !!! tmp.name = this->name + "_xySlice";
	for (size_t w = 0; w < size_w; w++) {
		tmp[w] = matrix_array[w][p_x][p_y][p_z];
	}
	return tmp;
}

/**
* Take an xyz-slice of 4d matrix and write it to out.
* \param p_x - index at which to slice x dimension
* \param p_y - index at which to slice y dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
void Matrix4D<T>::xyzSlice(Matrix1D<T>& out, size_t p_x, size_t p_y, size_t p_z) const {
	for (size_t w = 0; w < size_w; w++) {
		out[w] = matrix_array[w][p_x][p_y][p_z];
	}
}

CalculationMatrix::CalculationMatrix()
{
	initialized = false;
	change_ind = 0;
	total_size = 0;
	x_size = 0;
	y_size = 0;
}

/**
* Constructor for CalculationMatrix class.
* \param x_size - x size
* \param y_size - y size
* \param z_size - z size
* \param n_of_diags - NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL (main diagonal is not counted)
*/
CalculationMatrix::CalculationMatrix(int x_size, int y_size, int z_size, int n_of_diags) {
	this->initialized = false;
	Initialize(x_size, y_size, z_size, n_of_diags);
}

/**
* Allocating memory for CalculationMatrix,
* Setting the diagonals to be 0
*/
void CalculationMatrix::Initialize(int x_size, int y_size, int z_size, int n_of_diags) {
	this->initialized = false;
	this->x_size = x_size;
	this->y_size = y_size;
	this->total_size = x_size;
	if (y_size > 0) this->total_size = this->total_size * y_size;
	if (z_size > 0) this->total_size = this->total_size * z_size;

	// initializing diagonals
	int id;
	// !!!	if (z_size > 0) { // means 3d
	if (z_size > 1) { // means 3d
		for (int x = -n_of_diags; x <= n_of_diags; x++) {
			for (int y = -n_of_diags; y <= n_of_diags; y++) {
				for (int z = -n_of_diags; z <= n_of_diags; z++) {
					// calculating a diagonal number (id)
					id = this->index1d(x, y, z);
					// allocating memory for each diagonal
					(*this)[id] = Matrix1D<double>(this->total_size);
					(*this)[id] = 0;
				}
			}
		}
		// !!!	} else if (y_size > 0) {
	} else if (y_size > 1) {
		for (int x = -n_of_diags; x <= n_of_diags; x++) {
			for (int y = -n_of_diags; y <= n_of_diags; y++) {
				// calculating a diagonal number (id)
				id = this->index1d(x, y);
				// allocating memory for each diagonal
				(*this)[id] = Matrix1D<double>(this->total_size);
				(*this)[id] = 0;
			}
		}
		// !!!	} else if (x_size > 0) {
	} else if (x_size > 1) {
		for (int x = -n_of_diags; x <= n_of_diags; x++) {
			// calculating a diagonal number (id)
			id = this->index1d(x);
			// allocating memory for each diagonal
			(*this)[id] = Matrix1D<double>(this->total_size);
			(*this)[id] = 0;
		}
	}

	// If im_size is not set, then it's 1d and ia_size should be zero for future
	if (y_size == 0) this->x_size = 0;
	// If il_size is not set, then it's 2d and im_size should be zero for future
	if (z_size == 0) this->y_size = 0;

	this->initialized = true;
}

/**
* Function returns 1d index for 2D-3D arrays
*/
int CalculationMatrix::index1d(int x, int y, int z) {
	return (z*y_size + y)*x_size + x;
}

/**
* Save matrix to file.
* Includes varaible names and sizes
*/
void CalculationMatrix::writeToFile(const std::string& filename) const {
	int in;

	ofstream output(filename.c_str());
	//if (output==NULL && (filename.find("Debug") == std::string::npos)) {
	if (!output.is_open() && (filename.find("Debug") == std::string::npos)) {
		printf("FILE: Unable to output file: %s", filename.c_str());
		exit(EXIT_FAILURE);
	}
	output << "VARIABLES = \"";
	for (const auto& element : *this) {
		output << "\"" << element.first << "\", ";
	}
	output << endl;
	output << "ZONE T=\"" << "\", I=" << total_size << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (in = 0; in < this->total_size; in++) {
		for (const auto& element : *this) {
			output << "\t" << element.second[in];
		}
		output << endl;
	}
	output.close();
}

//////////////////////////////////////////
// Implementations
//////////////////////////////////////////





template class Matrix1D<double>;
template class Matrix2D<double>;
template class Matrix3D<double>;
template class Matrix4D<double>;

#endif

