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
	
#include "Matrix.h"

using namespace std;

#if defined(_WIN32) || defined(_WIN64)
	#define strncasecmp _strnicmp
	#define strcasecmp _stricmp
#endif

const double err = 1e-6;

// #define DEBUG_MODE

// Memory related functions

/// Allocating memory for 1D matrix
template<class T>inline T* matrix(long Rows)
{
	T *m=new T[Rows];
	//	assert(m!=NULL);
	if (m == NULL) {
		printf("MEMORY_ERRROR: Memory can't be initialized: %ld size", Rows*sizeof(T));
		exit(EXIT_FAILURE);
	}
	return m;
}

/// Initilizing memory for 2D matrix
template<class T>inline T** matrix(long Rows, long Columns)
{
	// allocating memory for array of pinters
	T **m=new T*[Rows];
	//	assert(m!=NULL);
	if (m == NULL) {
		printf("MEMORY_ERRROR: Memory can't be initialized: %ld size", Rows * sizeof(T));
		exit(EXIT_FAILURE);
	}
	// allocating memory for data array
	m[0] = new T[Rows * Columns];
	//	assert(m[0]!=NULL);
	if (m[0] == NULL) {
		printf("MEMORY_ERRROR: Memory can't be initialized: %ld size", Rows * Columns * sizeof(T));
		exit(EXIT_FAILURE);
	}
	// assign pointers to data ranges
	for(long i=1; i<Rows; i++) m[i] = m[i-1] + Columns;
	return m;
}

/// Initializing memory for 3D matrix
template<class T>inline T*** matrix(int size_x, int size_y, int size_z)
{
	// allocating memory for array of pointers to pointers
	T ***m=new T**[size_x];
	//	assert(m!=NULL);
	if (m == NULL) {
		printf("MEMORY_ERRROR: Memory can't be initialized: %d size", size_x * sizeof(T));
		exit(EXIT_FAILURE);
	}
	for (int x = 0; x < size_x; x++) {
		// for each pointer allocating memory for array of pointers
		m[x] = new T*[size_y];
		//		assert(m[x]!=NULL);
		if (m[x] == NULL) {
			printf("MEMORY_ERRROR: Memory can't be initialized: %d size", size_y * sizeof(T));
			exit(EXIT_FAILURE);
		}
	}
	//allocating memory for data array
	m[0][0] = new T[size_x * size_y * size_z];
	//	assert(m[0][0]!=NULL);
	if (m[0][0] == NULL) {
		printf("MEMORY_ERRROR: Memory can't be initialized: %d size", size_x * size_y * size_z * sizeof(T));
		exit(EXIT_FAILURE);
	}
	for (int x = 0; x < size_x; x++) {
		for (int y = 0; y < size_y; y++) {
			// assign pointers to data ranges
			m[x][y] = m[0][0] + (x*size_y + y)*size_z;
		}
	}

	return m;
}

/// Initializing memory for 4D matrix
template<class T>inline T**** matrix(int size_w, int size_x, int size_y, int size_z)
{
	// allocating memory for array of pointers to pointers
	T ****m=new T***[size_w];
	//	assert(m!=NULL);

	if (m == NULL) {
		printf("MEMORY_ERRROR: Memory can't be initialized: %d size", size_w * sizeof(T));
		exit(EXIT_FAILURE);
	}

	for (int w = 0; w < size_w; w++) {
		// for each pointer allocating memory for array of pointers
		m[w] = new T**[size_x];
		//		assert(m[x]!=NULL);
		if (m[w] == NULL) {
			printf("MEMORY_ERRROR: Memory can't be initialized: %d size", size_x * sizeof(T));
			exit(EXIT_FAILURE);
		}

		for (int x = 0; x < size_x; x++) {
			// for each pointer allocating memory for array of pointers
			m[w][x] = new T*[size_y];
			//		assert(m[x]!=NULL);
			if (m[w][x] == NULL) {
				printf("MEMORY_ERRROR: Memory can't be initialized: %d size", size_y * sizeof(T));
				exit(EXIT_FAILURE);
			}
		}
	}

	//allocating memory for data array
	m[0][0][0] = new T[size_w * size_x * size_y * size_z];
	//	assert(m[0][0]!=NULL);
	if (m[0][0][0] == NULL) {
		printf("MEMORY_ERRROR: Memory can't be initialized: %d size", size_w * size_x * size_y * size_z * sizeof(T));
		exit(EXIT_FAILURE);
	}
	for (int w = 0; w < size_w; w++) {
		for (int x = 0; x < size_x; x++) {
			for (int y = 0; y < size_y; y++) {
				// assign pointers to data ranges
				m[w][x][y] = m[0][0][0] + ((w*size_x + x)*size_y + y)*size_z;
			}
		}
	}

	return m;
}


/// Freeing memory for 1D matrix
template<class T>inline void free_matrix(T* m) {
	delete m;
}

/// Freeing memory for 2D matrix
template<class T>inline void free_matrix(T** m) {
	delete[](m[0]);
	delete[](m);
}

/// Freeing memory for 3D matrix
template<class T>inline void free_matrix(T*** m, int size_x, int size_y) {
	delete[](m[0][0]);
	for (int x = 0; x < size_x; x++) {
		delete[](m[x]);
	}
	delete[](m);
}


/// Freeing memory for 4D matrix
template<class T>inline void free_matrix(T**** m, int size_w, int size_x, int size_y) {
	delete[](m[0][0][0]);
	for (int w = 0; w < size_w; w++) {
		for (int x = 0; x < size_x; x++) {
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
Matrix1D<T>::Matrix1D( int size_q1 , string name) {
	initialized = false;
	// !!! this->name = name;
	// allocating memory
	AllocateMemory(size_q1);
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
void Matrix1D<T>::AllocateMemory( int size_q1 ) {
	this->size_q1 = size_q1;
	// using inline template for memory allocation
	matrix_array = matrix<T>(size_q1);
	initialized = true;
}


/**
* Operator [i], returns value of element i.
* If  DEBUG_MODE defined, check if matrix has been initialized.
*
* \param i - number of element to return
*/
template<class T>
inline T& Matrix1D<T>::operator[](int i1) {
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
inline T& Matrix1D<T>::operator[](int i1) const {
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
	int i1;
	for (i1 = 0; i1 < this->size_q1; i1++)
		matrix_array[i1] += M.matrix_array[i1];
	return *this;
}

/**
* Matrix subtraction, result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix1D<T>& Matrix1D<T>::operator-= (const Matrix1D<T> &M) {
	int i1;
	for (i1 = 0; i1 < this->size_q1; i1++)
		matrix_array[i1] -= M.matrix_array[i1];
	return *this;
}

/**
* Multiplication to a value. Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix1D<T>& Matrix1D<T>::operator*= (const T Val) {
	int i1;
	for (i1 = 0; i1 < this->size_q1; i1++)
		matrix_array[i1] *= Val;
	return *this;
}

/**
* Division by a value. Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix1D<T>& Matrix1D<T>::operator/= (const T Val) {
	int i1;
	for (i1 = 0; i1 < this->size_q1; i1++)
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
	for (int i1 = 0; i1 < this->size_q1; i1++)
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
	for (int i1 = 0; i1 < this->size_q1; i1++)
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
	int i1;
	Matrix1D<T> Tmp(*this);
	for (i1 = 0; i1 < this->size_q1; i1++)
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
	int i1;
	Matrix1D<T> Tmp(*this);
	for (i1 = 0; i1 < this->size_q1; i1++)
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
	int i1;
	Matrix1D<T> Tmp(*this);
	for (i1 = 0; i1 < this->size_q1; i1++)
		Tmp[i1] = matrix_array[i1] / M.matrix_array[i1];
	return Tmp;
}

/**
* Norm of vector
*/
template<class T>
inline T Matrix1D<T>::norm() const {
	T res = 0;
	int i1;
	for (i1 = 0; i1 < this->size_q1; i1++) {
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
	int i1;
	for (i1 = 0; i1 < this->size_q1; i1++) {
		res += matrix_array[i1] * W[i1];
	}
	return res;
}


/**
* Write matrix data to file.
*/
template<class T>
void Matrix1D<T>::writeToFile(string filename) {
	int i1;
	ofstream output(filename.c_str());
	//if (output==NULL && (filename.find("Debug") == string::npos)) {
	if (!output.is_open() && (filename.find("Debug") == string::npos)) {
		printf("FILE: Unable to output file: %s", filename.c_str());
		exit(EXIT_FAILURE);
	}
	output << "VARIABLES = \"" << ((this->name!="")?this->name:"function") << "\" "<< endl;
	output << "ZONE T=\"" << filename << "\", I=" << size_q1 << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (i1 = 0; i1 < size_q1; i1++) {
		output << matrix_array[i1] << endl;
	}
	output.close();
}

/**
* Write matrix data to file with grid.
*/
template<class T>
void Matrix1D<T>::writeToFile(string filename, Matrix1D<T> &grid_q1) {
	int i1;
	ofstream output(filename.c_str());
	//if (output==NULL && (filename.find("Debug") == string::npos)) {
	if (!output.is_open() && (filename.find("Debug") == string::npos)) {
		printf("FILE: Unable to output file: %s", filename.c_str());
		exit(EXIT_FAILURE);
	}
	output << "VARIABLES = \"" << ((grid_q1.name!="")?grid_q1.name:"x") << "\", \""  << ((this->name!="")?this->name:"function") << "\" "<< endl;
	output << "ZONE T=\"" << filename << "\", I=" << size_q1 << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (i1 = 0; i1 < size_q1; i1++) {
		output << "\t" << grid_q1[i1] << "\t" << "\t" << matrix_array[i1] << endl;
	}
	output.close();
}


/**
* Read matrix data from file.
*/
template<class T>
void Matrix1D<T>::readFromFile(string filename) {
	int i1;
	string inBuf;
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
			for (i1 = 0; i1 < size_q1; i1++) {
				input >> matrix_array[i1];
			}
		} else {
			printf("MATRIX_LOAD_ERROR: Error reading file %s.\n", filename.c_str());
			exit(EXIT_FAILURE);
		}
		input.close();
	}
}





// /**
// * Read matrix data from MATLAB file.
// */
// template<class T>
// void Matrix1D<T>::readFromMatlabFile(string filename) {
	
// 	MATfile *pmat;
// 	MxArray *pa;
// 	pmat = matOpen(filename.c_str(), "r");
// 	const char* name;
// 	pa = matGetNextVariable(pmat, &name);
// 	matrix_array = mxGetData(pa);
	
// }





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
void Matrix1D<T>::readFromFile(string filename, const Matrix1D<T> grid_q1) {
	int i1;
	string inBuf;
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
			for (i1 = 0; i1 < size_q1; i1++) {
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






// ADDED
/* Analyze field FNAME in struct array SPTR. */
static void
analyzestructarray(const mxArray *sPtr, const char *fName)
{
    mwSize nElements;       /* number of elements in array */
    mwIndex eIdx;           /* element index */
    mwIndex fPTRIdx;           /* fptr index */
    const mxArray *fPtr;    /* field pointer */
    double *realPtr;        /* pointer to data */
    mwSize nElementsInRealData; /* number of elements in array */
    
	
	// Goes through all of the structs
	// I think all of them only have one struct consisting of a couple fields (one or more of which is the double array which we want)
    nElements = (mwSize)mxGetNumberOfElements(sPtr);
    for (eIdx = 0; eIdx < nElements; eIdx++) {
        fPtr = mxGetField(sPtr, eIdx, fName);
        nElementsInRealData = (mwSize)mxGetNumberOfElements(fPtr);
        printf("number of elements in %s: %.2d\n", fName, nElementsInRealData);
        if ((fPtr != NULL)
            && (mxGetClassID(fPtr) == mxDOUBLE_CLASS) 
            && (!mxIsComplex(fPtr))) 
        {
            realPtr = mxGetPr(fPtr);
            for (fPTRIdx = 0; fPTRIdx < nElementsInRealData; fPTRIdx++)
            {
                printf("%.2f \n", realPtr[fPTRIdx]);
            }
        }
    }
}

/* Find struct array ARR in MAT-file FILE.
 * Pass field name FIELD to analyzestructarray function. */
int findstructure(
        const char *file,
        const char *arr,
        const char *field) {
    
   
    
}















/**
* Return minimum value of the 1d matrix.
* Default value set to 1e99
*/
template<class T>
T Matrix1D<T>::min() {
	T tmp = 1e99;
	int i1;
	for (i1 = 0; i1 < size_q1; i1++) {
		tmp = (tmp<matrix_array[i1])?tmp:matrix_array[i1];
	}
	return tmp;
}

/**
* Return maximum value of the 1d matrix.
* Default value seet to 0
*/
template<class T>
T Matrix1D<T>::max() {
	T tmp = 0;
	int i1;
	for (i1 = 0; i1 < size_q1; i1++) {
		tmp = (tmp>matrix_array[i1])?tmp:matrix_array[i1];
	}
	return tmp;
}

/**
* Return absolute maximum value of the 1d matrix.
* Default value set to 0
*/
template<class T>
T Matrix1D<T>::maxabs() {
	T tmp = 0;
	int i1;
	for (i1 = 0; i1 < size_q1; i1++) {
		tmp = (tmp>fabs((double)matrix_array[i1]))?tmp:fabs((double)matrix_array[i1]);
	}
	return tmp;
}

/**
* Return absolute value of the 1d matrix.
* Changes every element to a positive value with the same magnitude
*/
template<class T>
Matrix1D<T> Matrix1D<T>::abs() {
	Matrix1D<T> tmp(this->size_q1);
	int i1;
	for (i1 = 0; i1 < size_q1; i1++) {
		tmp[i1] = (matrix_array[i1]>0)?matrix_array[i1]:-matrix_array[i1];
	}
	return tmp;
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
Matrix2D<T>::Matrix2D( int size_q1, int size_q2 ) {
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
void Matrix2D<T>::AllocateMemory( int size_q1, int size_q2 ) {
	this->size_q1 = size_q1;
	this->size_q2 = size_q2;
	// using matrix inline template to allocate memory
	matrix_array = matrix<T>(size_q1, size_q2);
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
	int i1, i2;
	if (initialized) {
		for (i1 = 0; i1 < this->size_q1; i1++)
			for (i2 = 0; i2 < this->size_q2; i2++)
				matrix_array[i1][i2] = val;
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
	int i1, i2;
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i1 = 0; i2 < this->size_q2; i2++)
				matrix_array[i1][i2] += M.matrix_array[i1][i2];
	return *this;
}

/**
* Matrix subtraction, result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix2D<T>& Matrix2D<T>::operator-= (const Matrix2D<T> &M) {
	int i1, i2;
	for (i1 = 0; i1 < this->size_q1; i2++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			matrix_array[i1][i2] -= M.matrix_array[i1][i2];
	return *this;
}


/**
* Multiplication to a value. Result is stored into the left hand side matrix
*/
template<class T>
inline Matrix2D<T>& Matrix2D<T>::operator*= (const T Val) {
	int i1, i2;
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i1 = 0; i1 < this->size_q2; i2++)
				matrix_array[i1][i2] *= Val;
	return *this;
}


/**
* Multiply matrix to Val.
*
* \param Val - value Val.
*/
template<class T>
inline Matrix2D<T> Matrix2D<T>::operator* (const T Val) const {
	int i1, i2;
	Matrix2D<T> Tmp(*this);
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			Tmp[i1][i2] = matrix_array[i1][i2] * Val;
	return Tmp;
}

/**
* Divide matrix to Val.
*
* \param Val - value Val.
*/
template<class T>
inline Matrix2D<T> Matrix2D<T>::operator/(const T Val) const {
	int i1, i2;
	Matrix2D<T> Tmp(*this);
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			Tmp[i1][i2] = matrix_array[i1][i2] / Val;
	return Tmp;
}

/**
* Divide all values of matrix to values of matrix M.
*
* \param &M - matrix M.
*/
template<class T>
inline Matrix2D<T> Matrix2D<T>::divide (const Matrix2D<T> &M) const {
	int i1, i2;
	Matrix2D<T> Tmp(*this);
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			Tmp[i1][i2] = matrix_array[i1][i2] / M[i1][i2];
	return Tmp;
}

/**
* Multiply all values of matrix to values of matrix M.
*
* \param &M - matrix M.
*/
template<class T>
inline Matrix2D<T> Matrix2D<T>::times (const Matrix2D<T> &M) const {
	int i1, i2;
	Matrix2D<T> Tmp(*this);
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			Tmp[i1][i2] = matrix_array[i1][i2] * M[i1][i2];
	return Tmp;
}


/**
* Return maximum value of the 2d matrix.
* Default value seet to 0
*/
template<class T>
T Matrix2D<T>::max() {
	T tmp = 0;
	int i1, i2;
	for (i1 = 0; i1 < size_q1; i1++) {
		for (i2 = 0; i2 < size_q2; i2++) {
			tmp = (tmp>matrix_array[i1][i2])?tmp:matrix_array[i1][i2];
		}
	}
	return tmp;
}

/**
* Return absolute maximum value of the 2d matrix.
* Default value set to 0
*/
template<class T>
T Matrix2D<T>::maxabs() {
	T tmp = 0;
	int i1, i2;
	for (i1 = 0; i1 < size_q1; i1++) {
		for (i2 = 0; i2 < size_q2; i2++) {
			tmp = (tmp>fabs((double)matrix_array[i1][i2]))?tmp:fabs((double)matrix_array[i1][i2]);
		}
	}
	return tmp;
}


/**
* Return minimum value of the 2d matrix.
* Default value set to 1e99
*/
template<class T>
T Matrix2D<T>::min() {
	T tmp = 1e99;
	int i1, i2;
	for (i1 = 0; i1 < size_q1; i1++) {
		for (i2 = 0; i2 < size_q2; i2++) {
			tmp = (tmp<matrix_array[i1][i2])?tmp:matrix_array[i1][i2];
		}
	}
	return tmp;
}

/**
* Return absolute value of the 2d matrix.
* Changes every element to a positive value with the same magnitude
*/
template<class T>
Matrix2D<T> Matrix2D<T>::abs() {
	Matrix2D<T> tmp(this->size_q1, this->size_q2);
	int i1, i2;
	for (i1 = 0; i1 < size_q1; i1++) {
		for (i2 = 0; i2 < size_q2; i2++) {
			tmp[i1][i2] = (matrix_array[i1][i2]>0)?matrix_array[i1][i2]:-matrix_array[i1][i2];
		}
	}
	return tmp;
}


/**
* Divide all values of matrix to value Val.
*
* \param val - value to divide.
*/
template<class T>
Matrix2D<T> Matrix2D<T>::max_of(T val) {
	int i1, i2;
	Matrix2D<T> Tmp(*this);
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			Tmp[i1][i2] = (matrix_array[i1][i2]>val)?matrix_array[i1][i2]:val;
	return Tmp;
}

/**
* Returns corresponding index of 2d matrix if represented as a 1d array
* \param x,y - index of element in every dimension for the 2d matrix
*/
template<class T>
inline int Matrix2D<T>::index1d(int x, int y) const {
	return x*size_q2 + y;
}



/**
* Writes the matrix to a file.
* File has two header lines.
*
* \param filename - file name
*/
template<class T>
void Matrix2D<T>::writeToFile(string filename) {
	int i1, i2;
	ofstream output(filename.c_str());
	//if (output==NULL && (filename.find("Debug") == string::npos)) {
	if (!output.is_open() && (filename.find("Debug") == string::npos)) {
		printf("FILE: Unable to output file: %s", filename.c_str());
		exit(EXIT_FAILURE);
	}
	output << "VARIABLES = \""<< ((this->name!="")?this->name:"f") << "\" "<< endl;
	output << "ZONE T=\"" << filename << "\", I=" << size_q2 << ", J= " << size_q1 << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (i1 = 0; i1 < size_q1; i1++) {
		for (i2 = 0; i2 < size_q2; i2++) {
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
void Matrix2D<T>::writeToFile(string filename, Matrix2D<T> &grid_x, Matrix2D<T> &grid_y) {
	int i1, i2;
	ofstream output(filename.c_str());
	output << "VARIABLES = \"" << ((grid_x.name!="")?grid_x.name:"x") << "\", \"" << ((grid_y.name!="")?grid_y.name:"y") << "\", \"" << ((this->name!="")?this->name:"f") << "\" "<< endl;
	output << "ZONE T=\"" << filename << "\", I=" << size_q2 << ", J=" << size_q1 << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (i1 = 0; i1 < size_q1; i1++) {
		for (i2 = 0; i2 < size_q2; i2++) {
			output << "\t" << grid_x[i1][i2] << "\t" << grid_y[i1][i2] << "\t" << matrix_array[i1][i2] << endl;
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
void Matrix2D<T>::readFromFile(string filename, int read_column) {
	int i1, i2, column;
	string inBuf;
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
			for (i1 = 0; i1 < size_q1; i1++) {
				for (i2 = 0; i2 < size_q2; i2++) {
					for (column = 1; column < read_column; column++) input >> inBuf; // skip first columns
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
* Read matrix data from file with grid, 
* Checks if the matrix data in the file is the same as the grids that were sent in with error < 1e-8,
* if not within error range will signal error and exit
* 
* Overloaded readFromFile function
* \param filename - file to read grids from
* \param grids x,y - checks grids data against the file data
*/
template<class T>
void Matrix2D<T>::readFromFile(string filename, const Matrix2D<T> grid_x, const Matrix2D<T> grid_y) {
	int i1, i2;
	string inBuf;
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
			for (i1 = 0; i1 < size_q1; i1++) {
				for (i2 = 0; i2 < size_q2; i2++) {
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
* Make x-slice of 2d matrix = 1d matrix.
*/
template<class T>
Matrix1D<T> Matrix2D<T>::xSlice(int p_x) const {
	if (p_x > this->size_q1) {
		printf("MATRIX_ERROR: xSlice outside or array boundary: %d > %d", p_x, this->size_q1);
		exit(EXIT_FAILURE);
	}

	int i2;
	Matrix1D<T> tmp(this->size_q2);
	// !!! tmp.name = this->name + "_slice";
	for (i2 = 0; i2 < size_q2; i2++) {
		tmp[i2] = matrix_array[p_x][i2];
	}
	return tmp;
}

/**
* Make y-slice of 2d matrix = 1d matrix.
*/
template<class T>
Matrix1D<T> Matrix2D<T>::ySlice(int p_y) const {
	if (p_y > this->size_q2) {
		printf("MATRIX_ERROR: ySlice outside or array boundary: %d > %d", p_y, this->size_q2);
		exit(EXIT_FAILURE);
	}

	int i1;
	Matrix1D<T> tmp(this->size_q1);
	// !!! tmp.name = this->name + "_slice";
	for (i1 = 0; i1 < size_q1; i1++) {
		tmp[i1] = matrix_array[i1][p_y];
	}
	return tmp;
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
Matrix3D<T>::Matrix3D( int size_q1, int size_q2, int size_q3 ) {
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
	//		delete[] plane_array; - it's get deleted with matrix_array[0][0][0] or something?
}

/**
* Allocating memory and filling it with zero-values.
*/
template<class T>
void Matrix3D<T>::AllocateMemory( int size_q1, int size_q2, int size_q3) {
	this->size_q1 = size_q1;
	this->size_q2 = size_q2;
	this->size_q3 = size_q3;
	matrix_array = matrix<T>(size_q1, size_q2, size_q3);
	plane_array = matrix_array[0][0];
	initialized = true;
#ifdef DEBUG_MODE
	// should not initialize matrix with zeros, it can slow the code greatly in some cases
	for (int i = 0; i < size_q1; i++)
		for (int j = 0; j < size_q2; j++)
			for (int k = 0; k < size_q3; k++)
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
inline T** Matrix3D<T>::operator[] (int i1) {
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
	// std::fill_n(matrix_array[0][0][0], size_q1*size*y*size*z, Val); // might be faster
	int i1, i2, i3;
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				this->matrix_array[i1][i2][i3] = Val;
	return *this;
}


/**
* Matrix summation, result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::operator+= (const Matrix3D<T> &M) {
	int i1, i2, i3;
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				matrix_array[i1][i2][i3] += M.matrix_array[i1][i2][i3];
	return *this;
}

/**
* Matrix subtraction, result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::operator-= (const Matrix3D<T> &M) {
	int i1, i2, i3;
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				matrix_array[i1][i2][i3] -= M.matrix_array[i1][i2][i3];
	return *this;
}

/**
* Multiplication to a value. Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::operator*= (const T Val) {
	int i1, i2, i3;
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				matrix_array[i1][i2][i3] *= Val;
	return *this;
}

/**
* Division by a value. Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::operator/= (const T Val) {
	int i1, i2, i3;
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				matrix_array[i1][i2][i3] /= Val;
	return *this;
}

/**
* Summation with a value. Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::operator+= (const T Val) {
	int i1, i2, i3;
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				matrix_array[i1][i2][i3] += Val;
	return *this;
}

/**
* Subtraction of a value. Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::operator-= (const T Val) {
	int i1, i2, i3;
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				matrix_array[i1][i2][i3] -= Val;
	return *this;
}

/**
* Multiplication between each element of the matrices (not a matrix multiplication). Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::times_equal (const Matrix3D<T> &M) {
	int i1, i2, i3;
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				matrix_array[i1][i2][i3] *= M.matrix_array[i1][i2][i3];
	return *this;
}

/**
* Division of each element of one matrices to the element of another. Result is stored into applied matrix (left hand side matrix)
*/
template<class T>
inline Matrix3D<T>& Matrix3D<T>::divide_equal (const Matrix3D<T> &M) {
	int i1, i2, i3;
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				matrix_array[i1][i2][i3] /= M.matrix_array[i1][i2][i3];
	return *this;
}

/**
* Add each element of the matrix to corresponds element of matrix M.
*/
template<class T>
inline Matrix3D<T> Matrix3D<T>::operator+ (const Matrix3D<T> &M) const {
	int i1, i2, i3;
	Matrix3D<T> Tmp(*this);
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				Tmp[i1][i2][i3] = matrix_array[i1][i2][i3] + M.matrix_array[i1][i2][i3];
	return Tmp;
}

/**
* Substract each element of the matrix to corresponds element of matrix M.
*/
template<class T>
inline Matrix3D<T> Matrix3D<T>::operator- (const Matrix3D<T> &M) const {
	int i1, i2, i3;
	Matrix3D<T> Tmp(*this);
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				Tmp[i1][i2][i3] = matrix_array[i1][i2][i3] - M.matrix_array[i1][i2][i3];
	return Tmp;
}


/**
* Multiply each element of the matrix to Val, save result to a new matrix.
*/
template<class T>
inline Matrix3D<T> Matrix3D<T>::operator* (const T Val) const {
	int i1, i2, i3;
	Matrix3D<T> Tmp(*this);
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				Tmp[i1][i2][i3] = matrix_array[i1][i2][i3] * Val;
	return Tmp;
}


/**
* Divide each element of the matrix to Val, save result to a new matrix.
*/
template<class T>
inline Matrix3D<T> Matrix3D<T>::operator/ (const T Val) const {
	int i1, i2, i3;
	Matrix3D<T> Tmp(*this);
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				Tmp[i1][i2][i3] = matrix_array[i1][i2][i3] / Val;
	return Tmp;
}

/**
* Multiply each element of the matrix to corresponds element of matrix M.
*/
template<class T>
inline Matrix3D<T> Matrix3D<T>::times (const Matrix3D<T> &M) const {
	int i1, i2, i3;
	Matrix3D<T> Tmp(*this);
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				Tmp[i1][i2][i3] = matrix_array[i1][i2][i3] * M.matrix_array[i1][i2][i3];
	return Tmp;
}

/**
* Divide each element of the matrix to corresponds element of matrix M.
*/
template<class T>
inline Matrix3D<T> Matrix3D<T>::divide (const Matrix3D<T> &M) const {
	int i1, i2, i3;
	Matrix3D<T> Tmp(*this);
	for (i1 = 0; i1 < this->size_q1; i1++)
		for (i2 = 0; i2 < this->size_q2; i2++)
			for (i3 = 0; i3 < this->size_q3; i3++)
				Tmp[i1][i2][i3] = matrix_array[i1][i2][i3] / M.matrix_array[i1][i2][i3];
	return Tmp;
}


/**
* Write matrix to file.
* File has two header lines.
*/
template<class T>
void Matrix3D<T>::writeToFile(string filename, string info) {
	int i1, i2, i3;
	ofstream output(filename.c_str());
	output << "VARIABLES = \""<< ((this->name!="")?this->name:"f") <<"\" "<< endl;
	output << "ZONE T=\"" << ((info=="")?filename:info) << "\", I=" << size_q3 << ", J=" << size_q2 << ", K=" << size_q1 << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (i1 = 0; i1 < size_q1; i1++) {
		for (i2 = 0; i2 < size_q2; i2++) {
			for (i3 = 0; i3 < size_q3; i3++) {
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
void Matrix3D<T>::writeToFile(string filename, Matrix3D<T> &grid_x, Matrix3D<T> &grid_y, Matrix3D<T> &grid_z) {
	int i1, i2, i3;
	ofstream output(filename.c_str());
	//if (output==NULL && (filename.find("Debug") == string::npos)) {
	  if (!output.is_open() && (filename.find("Debug") == string::npos)) {
		printf("FILE: Unable to output file: %s", filename.c_str());
		exit(EXIT_FAILURE);
	}
	output << "VARIABLES = \"" << ((grid_x.name!="")?grid_x.name:"x") << "\", \"" << ((grid_y.name!="")?grid_y.name:"y") << "\", \"" << ((grid_z.name!="")?grid_z.name:"z") << "\", \"" << ((this->name!="")?this->name:"f") << "\" "<< endl;
	output << "ZONE T=\"" << filename << "\", I=" << size_q3 << ", J=" << size_q2 << ", K=" << size_q1 << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (i1 = 0; i1 < size_q1; i1++) {
		for (i2 = 0; i2 < size_q2; i2++) {
			for (i3 = 0; i3 < size_q3; i3++) {
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
void Matrix3D<T>::readFromFile(string filename, int read_column) {
	int i1, i2, i3;
	string inBuf;
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
			for (i1 = 0; i1 < size_q1; i1++) {
				for (i2 = 0; i2 < size_q2; i2++) {
					for (i3 = 0; i3 < size_q3; i3++) {
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
* Read matrix data from file with grid, 
* Checks if the matrix data in the file is the same as the grids that were sent in with error < 1e-8,
* if not within error range will signal error and exit
* 
* Overloaded readFromFile function
* \param filename - file to read grids from
* \param grids x,y,z - checks grids data against the file data
*/
template<class T>
void Matrix3D<T>::readFromFile(string filename, const Matrix3D<T> grid_x, const Matrix3D<T> grid_y, const Matrix3D<T> grid_z) {
	int i1, i2, i3;
	string inBuf;
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
			for (i1 = 0; i1 < size_q1; i1++) {
				for (i2 = 0; i2 < size_q2; i2++) {
					for (i3 = 0; i3 < size_q3; i3++) {
						//input >> grid_x[i1][i2][i3] >> grid_y[i1][i2][i3] >> grid_z[i1][i2][i3] >> matrix_array[i1][i2][i3];
						// skip till the end of the line
						//input.ignore(9999, '\n');

						 input >> loaded_x >> loaded_y >> loaded_z;
						 // check if grid is the same
						 if (fabs(log10(loaded_x) - log10(grid_x[i1][i2][i3])) > err || fabs(log10(loaded_y) - log10(grid_y[i1][i2][i3])) > err || fabs(log10(loaded_z) - log10(grid_z[i1][i2][i3])) > err) {
							printf("MATRIX_LOAD_GRID_ERR: Loading %s: grid mismatch [%d, %d, %d].\nLoaded: %e, %e, %e\nGrid: %e, %e, %e\n", filename.c_str(), i1, i2, i3, loaded_x, loaded_y, loaded_z, grid_x[i1][i2][i3], grid_y[i1][i2][i3], grid_z[i1][i2][i3]);
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
* Returns corresponding index of 3d matrix if represented as a 1d array
* \param x,y,z - index of element in every dimension for the 3d matrix
*/
template<class T>
inline int Matrix3D<T>::index1d(int x, int y, int z) {
	return (x*size_q2 + y)*size_q3 + z;
}

/**
* Return minimum value of the 3d matrix.
* Default value set to 1e99
*/
template<class T>
T Matrix3D<T>::min() {
	T tmp = 1e99;
	int i1, i2, i3;
	for (i1 = 0; i1 < size_q1; i1++) {
		for (i2 = 0; i2 < size_q2; i2++) {
			for (i3 = 0; i3 < size_q3; i3++) {
				tmp = (tmp<matrix_array[i1][i2][i3])?tmp:matrix_array[i1][i2][i3];
			}
		}
	}
	return tmp;
}

/**
* Return maximum value of the 3d matrix.
* Default value seet to 0
*/
template<class T>
T Matrix3D<T>::max() {
	T tmp = 0;
	int i1, i2, i3;
	for (i1 = 0; i1 < size_q1; i1++) {
		for (i2 = 0; i2 < size_q2; i2++) {
			for (i3 = 0; i3 < size_q3; i3++) {
				tmp = (tmp>matrix_array[i1][i2][i3])?tmp:matrix_array[i1][i2][i3];
			}
		}
	}
	return tmp;
}

/**
* Return absolute maximum value of the 3d matrix.
* Default value set to 0
*/
template<class T>
T Matrix3D<T>::maxabs() {
	T tmp = 0;
	int i1, i2, i3;
	for (i1 = 0; i1 < size_q1; i1++) {
		for (i2 = 0; i2 < size_q2; i2++) {
			for (i3 = 0; i3 < size_q3; i3++) {
				tmp = (tmp>fabs((double)matrix_array[i1][i2][i3]))?tmp:fabs((double)matrix_array[i1][i2][i3]);
			}
		}
	}
	return tmp;
}

/**
* Return absolute value of the 3d matrix.
* Changes every element to a positive value with the same magnitude
*/
template<class T>
Matrix3D<T> Matrix3D<T>::abs() {
	Matrix3D<T> tmp(this->size_q1, this->size_q2, this->size_q3);
	int i1, i2, i3;
	for (i1 = 0; i1 < size_q1; i1++) {
		for (i2 = 0; i2 < size_q2; i2++) {
			for (i3 = 0; i3 < size_q3; i3++) {
				tmp[i1][i2][i3] = (matrix_array[i1][i2][i3]>0)?matrix_array[i1][i2][i3]:-matrix_array[i1][i2][i3];
			}
		}
	}
	return tmp;
}


/**
* Take x-slice of 3d matrix turning it into a 2d matrix.
* \param p_x - index at which to slice x dimension
*/
template<class T>
Matrix2D<T> Matrix3D<T>::xSlice(int p_x) const {
	int i2, i3;
	Matrix2D<T> tmp(this->size_q2, this->size_q3);
	// !!! tmp.name = this->name + "_slice";
	for (i2 = 0; i2 < size_q2; i2++) {
		for (i3 = 0; i3 < size_q3; i3++) {
			tmp[i2][i3] = matrix_array[p_x][i2][i3];
		}
	}
	return tmp;
}

/**
* Take y-slice of 3d matrix turning it into a 2d matrix.
* \param p_y - index at which to slice y dimension
*/
template<class T>
Matrix2D<T> Matrix3D<T>::ySlice(int p_y) const {
	int i1, i3;
	Matrix2D<T> tmp(this->size_q1, this->size_q3);
	// !!! tmp.name = this->name + "_slice";
	for (i1 = 0; i1 < size_q1; i1++) {
		for (i3 = 0; i3 < size_q3; i3++) {
			tmp[i1][i3] = matrix_array[i1][p_y][i3];
		}
	}
	return tmp;
}

/**
* Take z-slice of 3d matrix turning it into a 2d matrix.
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix2D<T> Matrix3D<T>::zSlice(int p_z) const {
	int i1, i2;
	Matrix2D<T> tmp(this->size_q1, this->size_q2);
	// !!! tmp.name = this->name+ "_slice";
	for (i1 = 0; i1 < size_q1; i1++) {
		for (i2 = 0; i2 < size_q2; i2++) {
			tmp[i1][i2] = matrix_array[i1][i2][p_z];
		}
	}
	return tmp;
}

/**
* Take xy-slice of 3d matrix turning it into a 1d matrix.
* \param p_x - index at which to slice x dimension
* \param p_y - index at which to slice y dimension
*/
template<class T>
Matrix1D<T> Matrix3D<T>::xySlice(int p_x, int p_y) const {
	int i3;
	Matrix1D<T> tmp(this->size_q3);
	// !!! tmp.name = this->name + "_xySlice";
	for (i3 = 0; i3 < size_q3; i3++) {
		tmp[i3] = matrix_array[p_x][p_y][i3];
	}
	return tmp;
}

/**
* Take yz-slice of 3d matrix turning it into a 1d matrix.
* \param p_y - index at which to slice y dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix1D<T> Matrix3D<T>::yzSlice(int p_y, int p_z) const {
	int i1;
	Matrix1D<T> tmp(this->size_q1);
	// !!! tmp.name = this->name + "_xzSlice";
	for (i1 = 0; i1 < size_q1; i1++) {
		tmp[i1] = matrix_array[i1][p_y][p_z];
	}
	return tmp;
}

/**
* Take xz-slice of 3d matrix turning it into a 1d matrix.
* \param p_x - index at which to slice x dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix1D<T> Matrix3D<T>::xzSlice(int p_x, int p_z) const {
	int i2;
	Matrix1D<T> tmp(this->size_q2);
	// !!! tmp.name = this->name + "_xzSlice";
	for (i2 = 0; i2 < size_q2; i2++) {
		tmp[i2] = matrix_array[p_x][i2][p_z];
	}
	return tmp;
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
Matrix4D<T>::Matrix4D( int w_size, int x_size, int y_size, int z_size ) {
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
void Matrix4D<T>::AllocateMemory( int w_size, int x_size, int y_size, int z_size) {
	this->size_w = w_size;
	this->size_x = x_size;
	this->size_y = y_size;
	this->size_z = z_size;
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
inline T*** Matrix4D<T>::operator[] (int i) {
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
	return plane_array[((w*size_x + x)*size_y + y)*size_z];
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
			this->name = string(M.name.length(), ' ');
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
	// std::fill_n(matrix_array[0][0][0], size_x*size*y*size*z, Val); // might be faster
	for (int w = 0; w < this->size_w; w++)
		for (int x = 0; x < this->size_x; x++)
			for (int y = 0; y < this->size_y; y++)
				for (int z = 0; z < this->size_z; z++)
					this->matrix_array[w][x][y][z] = Val;
	return *this;
}


/**
* Add each element of the current matrix by the corresponding element in M and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::operator+= (const Matrix4D<T> &M) {
	int w, x, y, z;
	for (w = 0; w < this->size_w; w++)
		for (x = 0; x < this->size_x; x++)
			for (y = 0; y < this->size_y; y++)
				for (z = 0; z < this->size_z; z++)
					matrix_array[w][x][y][z] += M.matrix_array[w][x][y][z];
	return *this;
}

/**
* Subtract each element of the current matrix by the corresponding element in M and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::operator-= (const Matrix4D<T> &M) {
	int w, x, y, z;
	for (w = 0; w < this->size_w; w++)
		for (x = 0; x < this->size_x; x++)
			for (y = 0; y < this->size_y; y++)
				for (z = 0; z < this->size_z; z++)
					matrix_array[w][x][y][z] -= M.matrix_array[w][x][y][z];
	return *this;
}

/**
* Multiply each element of the current matrix by Val and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::operator*= (const T Val) {
	int w, x, y, z;
	for (w = 0; w < this->size_w; w++)
		for (x = 0; x < this->size_x; x++)
			for (y = 0; y < this->size_y; y++)
				for (z = 0; z < this->size_z; z++)
					matrix_array[w][x][y][z] *= Val;
	return *this;
}

/**
* Divide each element of the current matrix by Val and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::operator/= (const T Val) {
	int w, x, y, z;
	for (w = 0; w < this->size_w; w++)
		for (x = 0; x < this->size_x; x++)
			for (y = 0; y < this->size_y; y++)
				for (z = 0; z < this->size_z; z++)
					matrix_array[w][x][y][z] /= Val;
	return *this;
}

/**
* Add each element of the current matrix by Val and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::operator+= (const T Val) {
	int w, x, y, z;
	for (w = 0; w < this->size_w; w++)
		for (x = 0; x < this->size_x; x++)
			for (y = 0; y < this->size_y; y++)
				for (z = 0; z < this->size_z; z++)
					matrix_array[w][x][y][z] += Val;
	return *this;
}

/**
* Subtract each element of the current matrix by Val and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::operator-= (const T Val) {
	int w, x, y, z;
	for (w = 0; w < this->size_w; w++)
		for (x = 0; x < this->size_x; x++)
			for (y = 0; y < this->size_y; y++)
				for (z = 0; z < this->size_z; z++)
					matrix_array[w][x][y][z] -= Val;
	return *this;
}

/**
* Multiply each element of the current matrix by the corresponding element of matrix M and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::times_equal (const Matrix4D<T> &M) {
	int w, x, y, z;
	for (w = 0; w < this->size_w; w++)
		for (x = 0; x < this->size_x; x++)
			for (y = 0; y < this->size_y; y++)
				for (z = 0; z < this->size_z; z++)
					matrix_array[w][x][y][z] *= M.matrix_array[w][x][y][z];
	return *this;
}

/**
* Divide each element of the current matrix by the corresponding element of matrix M and return this matrix
*/
template<class T>
inline Matrix4D<T>& Matrix4D<T>::divide_equal (const Matrix4D<T> &M) {
	int w, x, y, z;
	for (w = 0; w < this->size_w; w++)
		for (x = 0; x < this->size_x; x++)
			for (y = 0; y < this->size_y; y++)
				for (z = 0; z < this->size_z; z++)
					matrix_array[w][x][y][z] /= M.matrix_array[w][x][y][z];
	return *this;
}

/**
* Add each element of the matrix by the corresponding element of matrix M and return new matrix
*/
template<class T>
inline Matrix4D<T> Matrix4D<T>::operator+ (const Matrix4D<T> &M) const {
	int w, x, y, z;
	Matrix4D<T> Tmp(*this);
	for (w = 0; w < this->size_w; w++)
		for (x = 0; x < this->size_x; x++)
			for (y = 0; y < this->size_y; y++)
				for (z = 0; z < this->size_z; z++)
					Tmp[w][x][y][z] = matrix_array[w][x][y][z] + M.matrix_array[w][x][y][z];
	return Tmp;
}

/**
* Substract each element of the matrix by the corresponding element of matrix M and return new matrix
*/
template<class T>
inline Matrix4D<T> Matrix4D<T>::operator- (const Matrix4D<T> &M) const {
	int w, x, y, z;
	Matrix4D<T> Tmp(*this);
	for (w = 0; w < this->size_w; w++)
		for (x = 0; x < this->size_x; x++)
			for (y = 0; y < this->size_y; y++)
				for (z = 0; z < this->size_z; z++)
					Tmp[w][x][y][z] = matrix_array[w][x][y][z] - M.matrix_array[w][x][y][z];
	return Tmp;
}


/**
* Multiply each element of the matrix by Val, save result to a new matrix.
*/
template<class T>
inline Matrix4D<T> Matrix4D<T>::operator* (const T Val) const {
	int w, x, y, z;
	Matrix4D<T> Tmp(*this);
	for (w = 0; w < this->size_w; w++)
		for (x = 0; x < this->size_x; x++)
			for (y = 0; y < this->size_y; y++)
				for (z = 0; z < this->size_z; z++)
					Tmp[w][x][y][z] = matrix_array[w][x][y][z] * Val;
	return Tmp;
}


/**
* Divide each element of the matrix by Val, save result to a new matrix.
*/
template<class T>
inline Matrix4D<T> Matrix4D<T>::operator/ (const T Val) const {
	int w, x, y, z;
	Matrix4D<T> Tmp(*this);
	for (w = 0; w < this->size_w; w++)
		for (x = 0; x < this->size_x; x++)
			for (y = 0; y < this->size_y; y++)
				for (z = 0; z < this->size_z; z++)
					Tmp[w][x][y][z] = matrix_array[w][x][y][z] / Val;
	return Tmp;
}

/**
* Multiply each element of the matrix by the corresponding element of matrix M and return new matrix
*/
template<class T>
inline Matrix4D<T> Matrix4D<T>::times (const Matrix4D<T> &M) const {
	int w, x, y, z;
	Matrix4D<T> Tmp(*this);
	for (w = 0; w < this->size_w; w++)
		for (x = 0; x < this->size_x; x++)
			for (y = 0; y < this->size_y; y++)
				for (z = 0; z < this->size_z; z++)
					Tmp[w][x][y][z] = matrix_array[w][x][y][z] * M.matrix_array[w][x][y][z];
	return Tmp;
}

/**
* Divide each element of the current matrix by the corresponding element of matrix M and return new matrix
*/
template<class T>
inline Matrix4D<T> Matrix4D<T>::divide (const Matrix4D<T> &M) const {
	int w, x, y, z;
	Matrix4D<T> Tmp(*this);
	for (w = 0; w < this->size_w; w++)
		for (x = 0; x < this->size_x; x++)
			for (y = 0; y < this->size_y; y++)
				for (z = 0; z < this->size_z; z++)
					Tmp[w][x][y][z] = matrix_array[w][x][y][z] / M.matrix_array[w][x][y][z];
	return Tmp;
}


/**
* Write matrix to file.
* File has two header lines.
*/
template<class T>
void Matrix4D<T>::writeToFile(string filename, string info) {
	int w, x, y, z;
	ofstream output(filename.c_str());
	output << "VARIABLES = \""<< ((this->name!="")?this->name:"f") <<"\" "<< endl;
	output << "ZONE T=\"" << ((info=="")?filename:info) << "\", I=" << size_z << ", J=" << size_y << ", K=" << size_x << ", L=" << size_w << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (w = 0; w < size_w; w++)
		for (x = 0; x < size_x; x++)
			for (y = 0; y < size_y; y++)
				for (z = 0; z < size_z; z++)
					output << matrix_array[w][x][y][z] << endl;
	output.close();
}


/**
* Write matrix to file, using 3 other matrixes as a grid (simply - write all 4 matrixes to the file).
* File has two header lines.
*/
template<class T>
void Matrix4D<T>::writeToFile(string filename, Matrix4D<T> &grid_w, Matrix4D<T> &grid_x, Matrix4D<T> &grid_y, Matrix4D<T> &grid_z) {
	int w, x, y, z;
	ofstream output(filename.c_str());
	//if (output==NULL && (filename.find("Debug") == string::npos)) {
	if (!output.is_open() && (filename.find("Debug") == string::npos)) {
		printf("FILE: Unable to output file: %s", filename.c_str());
		exit(EXIT_FAILURE);
	}
	output << "VARIABLES = \"" << ((grid_w.name!="")?grid_w.name:"w") << "\", \"" << ((grid_x.name!="")?grid_x.name:"x") << "\", \"" << ((grid_y.name!="")?grid_y.name:"y") << "\", \"" << ((grid_z.name!="")?grid_z.name:"z") << "\", \"" << ((this->name!="")?this->name:"f") << "\" "<< endl;
	output << "ZONE T=\"" << filename << "\", W=" << size_w << ", I=" << size_z << ", J=" << size_y << ", K=" << size_x << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (w = 0; w < size_w; w++) {
		for (x = 0; x < size_x; x++) {
			for (y = 0; y < size_y; y++) {
				for (z = 0; z < size_z; z++) {
					output << "\t" << grid_x[w][x][y][z] << "\t" << grid_y[w][x][y][z] << "\t" << grid_z[w][x][y][z] << "\t" << matrix_array[w][x][y][z] << endl;
				}
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
void Matrix4D<T>::readFromFile(string filename, int read_column) {
	int w, x, y, z;
	string inBuf;
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
					std::istringstream linestream(inBuf);
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
			for (w = 0; w < size_w; w++) {
				for (x = 0; x < size_x; x++) {
					for (y = 0; y < size_y; y++) {
						for (z = 0; z < size_z; z++) {
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


// Function for reading from matlab file in 1-dimension
// Will always be stored in PRVK format so 1 = P, 2 = R ... 
template<class T>
void Matrix4D<T>::readFromMatlabFile(string file ,int columnNumber)
{
	MATFile *mfPtr; /* MAT-file pointer */
    mxArray *aPtr;  /* mxArray pointer */
    double *realPtr; /* pointer to data */
	string arr; /*name of variable*/
	string field = "arr"; // name of field
	mwSize nElements;       /* number of elements in array */
    mwIndex eIdx;           /* element index */
    const mxArray *fPtr;    /* field pointer */
    int w,x,y,z; 			/* for index*/
	
	
	if (columnNumber == 1)
		arr =  "P";
	if (columnNumber == 2)
		arr =  "R";
	if (columnNumber == 3)
		arr =  "V";
	if (columnNumber == 4)
		arr =  "K";
	
	this->name = arr;
		
	mfPtr = matOpen(file.c_str(), "r");
    if (mfPtr == NULL) {
        printf("Error opening file %s\n", file);
    }
    
    aPtr = matGetVariable(mfPtr, arr.c_str());
    if (aPtr == NULL) {
        printf("mxArray not found: %s\n", arr);
    }
    
    if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
        if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
            printf("Field not found: %s\n", field);
        }
        else {
            // Goes through all of the structs
			// I think all of them only have one struct consisting of a couple fields (one or more of which is the double array which we want)
    		nElements = (mwSize)mxGetNumberOfElements(aPtr);
   			 for (eIdx = 0; eIdx < nElements; eIdx++) {
       		 	fPtr = mxGetField(aPtr, eIdx, field.c_str()); // field was previously fname
       		 	if ((fPtr != NULL) && (mxGetClassID(fPtr) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtr))) 
       		 	{
          			realPtr = mxGetPr(fPtr);
       		 	}
   			 }
        }
    } 
   // ADDED
   // change to call function if it is just an array of doubles
   // realPtr should hold a pointer to the double array
    else if (mxGetClassID(aPtr) == mxDOUBLE_CLASS) 
	{
        realPtr = mxGetPr(aPtr);
	}
	else 
	{
		printf("%s is of unknown type\n", arr);
    }

	// sets the matrix array to be equal to an array of doubles
	for (w = 0; w < size_w; w++) {
		for (x = 0; x < size_x; x++) {
			for (y = 0; y < size_y; y++) {
				for (z = 0; z < size_z; z++) {
					matrix_array[w][x][y][z] = realPtr[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w];
				}
			}
		}
	}
	
    mxDestroyArray(aPtr);
    
    if (matClose(mfPtr) != 0) {
        printf("Error closing file %s\n", file);
    }
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
void Matrix4D<T>::readFromFile(string filename, const Matrix4D<T> grid_w, const Matrix4D<T> grid_x, const Matrix4D<T> grid_y, const Matrix4D<T> grid_z) {
	int w, x, y, z;
	string inBuf;
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
					std::istringstream linestream(inBuf);
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
			for (w = 0; w < size_w; w++) {
				for (x = 0; x < size_x; x++) {
					for (y = 0; y < size_y; y++) {
						for (z = 0; z < size_z; z++) {
							//input >> grid_x[x][y][z] >> grid_y[x][y][z] >> grid_z[x][y][z] >> matrix_array[x][y][z];
							// skip till the end of the line
							//input.ignore(9999, '\n');

							 input >> loaded_w >> loaded_x >> loaded_y >> loaded_z;
							 // check if grid is the same
							 if (fabs(log10(loaded_w) - log10(grid_w[w][x][y][z])) > err || fabs(log10(loaded_x) - log10(grid_x[w][x][y][z])) > err || fabs(log10(loaded_y) - log10(grid_y[w][x][y][z])) > err || fabs(log10(loaded_z) - log10(grid_z[w][x][y][z])) > err) {
								printf("MATRIX_LOAD_GRID_ERR: Loading %s: grid mismatch [%d, %d, %d, %d].\nLoaded: %e, %e, %e, %e\nGrid: %e, %e, %e, %e\n", filename.c_str(), w, x, y, z, loaded_w, loaded_x, loaded_y, loaded_z, grid_w[w][x][y][z], grid_x[w][x][y][z], grid_y[w][x][y][z], grid_z[w][x][y][z]);
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





// Function for reading from matlab file in 1-dimension
// Will always be stored in PRVK format so 1 = P, 2 = R ... 
template<class T>
void Matrix4D<T>::readFromMatlabFile(string file , const Matrix4D<T> grid_w, const Matrix4D<T> grid_x, const Matrix4D<T> grid_y, const Matrix4D<T> grid_z)
{
	MATFile *mfPtr; /* MAT-file pointer */
    mxArray *aPtr1;  /* mxArray pointer */
	mxArray *aPtr2;  /* mxArray pointer */
	mxArray *aPtr3;  /* mxArray pointer */
	mxArray *aPtr4;  /* mxArray pointer */
	mxArray *aPtr5;  /* mxArray pointer */
	mxArray *aPtrfinished;  /* mxArray pointer */
    double *realPtr; /* pointer to data */
	double *realPtrLoadedW; /* loaded W grid */
	double *realPtrLoadedX; /* loaded X grid */
	double *realPtrLoadedY; /* loaded Y grid */
	double *realPtrLoadedZ; /* loaded Z grid */
	string arr; /*name of variable*/
	string field = "arr"; // name of field
	mwSize nElements;       /* number of elements in array */
    mwIndex eIdx;           /* element index */
    const mxArray *fPtr;    /* field pointer */
	const mxArray *fPtrw;    /* field pointer */
	const mxArray *fPtrx;    /* field pointer */
	const mxArray *fPtry;    /* field pointer */
	const mxArray *fPtrz;    /* field pointer */
    int w,x,y,z; 			/* for index*/
	const char* name;		/* for getting variable names */		
	printf("reached here3");
			
	mfPtr = matOpen(file.c_str(), "r");
    if (mfPtr == NULL) {
        printf("Error opening file %s\n", file);
    }
    
	printf("reached here2");
	aPtr1 = matGetNextVariableInfo(mfPtr, &name);
	aPtr2 = matGetNextVariableInfo(mfPtr, &name);
	aPtr3 = matGetNextVariableInfo(mfPtr, &name);
	aPtr4 = matGetNextVariableInfo(mfPtr, &name);
	aPtr5 = matGetNextVariableInfo(mfPtr, &name);
	this->name = name;
    
	aPtrfinished = matGetVariable(mfPtr, name);
    if (aPtrfinished == NULL) {
        printf("mxArray not found: %s\n", arr);
    }
	
	
     printf("reached here1");	
    if (mxGetClassID(aPtrfinished) == mxSTRUCT_CLASS) {
        if (mxGetFieldNumber(aPtrfinished, field.c_str()) == -1) {
            printf("Field not found: %s\n", field);
        }
        else {
            // Goes through all of the structs
			// I think all of them only have one struct consisting of a couple fields (one or more of which is the double array which we want)
    		nElements = (mwSize)mxGetNumberOfElements(aPtrfinished);
   			 for (eIdx = 0; eIdx < nElements; eIdx++) {
       		 	fPtr = mxGetField(aPtrfinished, eIdx, field.c_str()); // field was previously fname
       		 	if ((fPtr != NULL) && (mxGetClassID(fPtr) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtr))) 
       		 	{
          			printf("size of realptr = %d \n" , (int)mxGetNumberOfElements(fPtr) );  
					realPtr = mxGetPr(fPtr);  
       		 	}
   			 }
        }
    } 
   // ADDED
   // change to call function if it is just an array of doubles
   // realPtr should hold a pointer to the double array
    else if (mxGetClassID(aPtrfinished) == mxDOUBLE_CLASS) 
	{
        realPtr = mxGetPr(aPtrfinished);
	}
	else 
	{
		printf("%s is of unknown type\n", arr);
    }


	// // // NOW GET ALL THE LOADED GRIDS TO CHECK ERROR
		
	// 	arr = grid_w.name;
    // 	aPtr = matGetVariable(mfPtr, arr.c_str());
    // 	if (aPtr == NULL) {
    //     	printf("mxArray not found: %s\n", arr);
   	// 	}   
		   
	// 	if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
    //     	if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
    //         	printf("Field not found: %s\n", field);
    //     	}
    //     	else 
	// 		{
	// 			nElements = (mwSize)mxGetNumberOfElements(aPtr);
   	// 			for (eIdx = 0; eIdx < nElements; eIdx++) {
	// 				fPtrw = mxGetField(aPtr, eIdx, field.c_str());
	//        		 	if ((fPtrw != NULL) && (mxGetClassID(fPtrw) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrw))) 
    //    		 		{
    //       				realPtrLoadedW = mxGetPr(fPtrw);
    //    		 		}
	// 			   }
	// 		}
    // 	}
		
	// 	arr = grid_x.name;
    // 	aPtr = matGetVariable(mfPtr, arr.c_str());
    // 	if (aPtr == NULL) {
    //     	printf("mxArray not found: %s\n", arr);
   	// 	}   
		   
	// 	if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
    //     	if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
    //         	printf("Field not found: %s\n", field);
    //     	}
    //     	else 
	// 		{
	// 			nElements = (mwSize)mxGetNumberOfElements(aPtr);
   	// 			 for (eIdx = 0; eIdx < nElements; eIdx++) {
	// 				fPtrx = mxGetField(aPtr, eIdx, field.c_str());
    //    		 		if ((fPtrx != NULL) && (mxGetClassID(fPtrx) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrx))) 
    //    		 		{
    //       				realPtrLoadedX = mxGetPr(fPtrx);
    //    		 		}
	// 			}
	// 		}
    // 	}
		
	// 	arr = grid_y.name;
    // 	aPtr = matGetVariable(mfPtr, arr.c_str());
    // 	if (aPtr == NULL) {
    //     	printf("mxArray not found: %s\n", arr);
   	// 	}   
		   
	// 	if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
    //     	if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
    //         	printf("Field not found: %s\n", field);
    //     	}
    //     	else 
	// 		{
	// 			nElements = (mwSize)mxGetNumberOfElements(aPtr);
   	// 		 	for (eIdx = 0; eIdx < nElements; eIdx++) {
	// 				fPtry = mxGetField(aPtr, eIdx, field.c_str());
    //    		 		if ((fPtry != NULL) && (mxGetClassID(fPtry) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtry))) 
    //    		 		{
    //       				realPtrLoadedY = mxGetPr(fPtry);
    //    		 		}
	// 			}
	// 		}
    // 	}
		
	// 	arr = grid_z.name;
    // 	aPtr = matGetVariable(mfPtr, arr.c_str());
    // 	if (aPtr == NULL) {
    //     	printf("mxArray not found: %s\n", arr);
   	// 	}   
		   
	// 	if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
    //     	if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
    //         	printf("Field not found: %s\n", field);
    //     	}
    //     	else 
	// 		{
	// 			nElements = (mwSize)mxGetNumberOfElements(aPtr);
   	// 		 	for (eIdx = 0; eIdx < nElements; eIdx++) {
	// 				fPtrz = mxGetField(aPtr, eIdx, field.c_str());
    //    		 		if ((fPtrz != NULL) && (mxGetClassID(fPtrz) == mxDOUBLE_CLASS) && (!mxIsComplex(fPtrz))) 
    //    		 		{
    //       				realPtrLoadedZ = mxGetPr(fPtrz);
    //    		 		}
	// 			}
	// 		}
    // 	}
		
		
	
printf("w: %d x: %d y: %d z: %d \n", size_w, size_x, size_y, size_z);

	// sets the matrix array to be equal to an array of doubles
	for (w = 0; w < size_w; w++) {
		for (x = 0; x < size_x; x++) {
			for (y = 0; y < size_y; y++) {
				for (z = 0; z < size_z; z++) {
					
					
					// if (fabs(log10(realPtrLoadedW[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w]) - log10(grid_w[w][x][y][z])) > err || fabs(log10(realPtrLoadedX[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w]) - log10(grid_x[w][x][y][z])) > err || fabs(log10(realPtrLoadedY[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w]) - log10(grid_y[w][x][y][z])) > err || fabs(log10(realPtrLoadedZ[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w]) - log10(grid_z[w][x][y][z])) > err) {
	// 				// 			printf("MATRIX_LOAD_GRID_ERR: Loading %s: grid mismatch [%d, %d, %d, %d].\nLoaded: %e\nGrid: %e, %e, %e, %e\n", file.c_str(), w, x, y, z, realPtr[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w], grid_w[w][x][y][z], grid_x[w][x][y][z], grid_y[w][x][y][z], grid_z[w][x][y][z]);
					matrix_array[w][x][y][z] = realPtr[z*(size_x * size_y * size_w) + y*(size_x * size_w) +  x*(size_w) + w];
					// }
				}
			}
		}
	}
	
    mxDestroyArray(aPtr1);
	mxDestroyArray(aPtr2);
	mxDestroyArray(aPtr3);	
	mxDestroyArray(aPtr4);
	mxDestroyArray(aPtr5);
	mxDestroyArray(aPtrfinished);

    
    if (matClose(mfPtr) != 0) {
        printf("Error closing file %s\n", file);
    }
}









/**
* Returns corresponding index of 4d matrix if represented as a 1d array
* \param w,x,y,z - index of element in every dimension for the 4d matrix
*/
template<class T>
inline int Matrix4D<T>::index1d(int w, int x, int y, int z) {
	return ((w*size_x + x)*size_y + y)*size_z;
}



/**
* Return minimum value of the 4d matrix.
* Default value set to 1e99
*/
template<class T>
T Matrix4D<T>::min() {
	T tmp = 1e99;
	int w, x, y, z;
	for (w = 0; w < size_w; w++) {
		for (x = 0; x < size_x; x++) {
			for (y = 0; y < size_y; y++) {
				for (z = 0; z < size_z; z++) {
					tmp = (tmp<matrix_array[w][x][y][z])?tmp:matrix_array[w][x][y][z];
				}
			}
		}
	}
	return tmp;
}

/**
* Return maximum value of the 4d matrix.
* Default value set to 0
*/
template<class T>
T Matrix4D<T>::max() {
	T tmp = 0;
	int w, x, y, z;
	for (w = 0; w < size_w; w++) {
		for (x = 0; x < size_x; x++) {
			for (y = 0; y < size_y; y++) {
				for (z = 0; z < size_z; z++) {
					tmp = (tmp>matrix_array[w][x][y][z])?tmp:matrix_array[w][x][y][z];
				}
			}
		}
	}
	return tmp;
}

/**
* Return absolute maximum value of the 4d matrix.
* Default value set to 0
*/
template<class T>
T Matrix4D<T>::maxabs() {
	T tmp = 0;
	int w, x, y, z;
	for (w = 0; w < size_w; w++) {
		for (x = 0; x < size_x; x++) {
			for (y = 0; y < size_y; y++) {
				for (z = 0; z < size_z; z++) {
					tmp = (tmp>fabs((double)matrix_array[w][x][y][z]))?tmp:fabs((double)matrix_array[w][x][y][z]);
				}
			}
		}
	}
	return tmp;
}

/**
* Return absolute value of the 4d matrix.
* Changes every element to a positive value with the same magnitude
*/
template<class T>
Matrix4D<T> Matrix4D<T>::abs() {
	Matrix4D<T> tmp(this->size_w, this->size_x, this->size_y, this->size_z);
	int w, x, y, z;
	for (w = 0; w < size_w; w++) {
		for (x = 0; x < size_x; x++) {
			for (y = 0; y < size_y; y++) {
				for (z = 0; z < size_z; z++) {
					tmp[w][x][y][z] = (matrix_array[w][x][y][z]>0)?matrix_array[w][x][y][z]:-matrix_array[w][x][y][z];
				}
			}
		}
	}
	return tmp;
}

/**
* Take w-slice of 4d matrix turning it into 3d matrix.
* squeeze matrix to 3d getting all the data when the w dimension is at index p_w
* \param p_w - index at which to slice w dimension
*/
template<class T>
Matrix3D<T> Matrix4D<T>::wSlice(int p_w) const {
	int x, y, z;
	Matrix3D<T> tmp(this->size_x, this->size_y, this->size_z);
	// !!! tmp.name = this->name + "_slice";
	for (x = 0; x < size_x; x++) {
		for (y = 0; y < size_y; y++) {
			for (z = 0; z < size_z; z++) {
				tmp[x][y][z] = matrix_array[p_w][x][y][z];
			}
		}
	}
	return tmp;
}

/**
* Take x-slice of 4d matrix turining it into 3d matrix.
* squeeze matrix to 3d getting all the data when the x dimension is at index p_x
* \param p_x - index at which to slice x dimension
*/
template<class T>
Matrix3D<T> Matrix4D<T>::xSlice(int p_x) const {
	int w, y, z;
	Matrix3D<T> tmp(this->size_w, this->size_y, this->size_z);
	// !!! tmp.name = this->name + "_slice";
	for (w = 0; w < size_w; w++) {
		for (y = 0; y < size_y; y++) {
			for (z = 0; z < size_z; z++) {
				tmp[w][y][z] = matrix_array[w][p_x][y][z];
			}
		}
	}
	return tmp;
}

/**
* Take y-slice of 4d matrix turning it into 3d matrix.
* squeeze matrix to 3d getting all the data when the y dimension is at index p_y
* \param p_y - index at which to slice y dimension
*/
template<class T>
Matrix3D<T> Matrix4D<T>::ySlice(int p_y) const {
	int w, x, z;
	Matrix3D<T> tmp(this->size_w, this->size_x, this->size_z);
	// !!! tmp.name = this->name + "_slice";
	for (w = 0; w < size_w; w++) {
		for (x = 0; x < size_x; x++) {
			for (z = 0; z < size_z; z++) {
				tmp[w][x][z] = matrix_array[w][x][p_y][z];
			}
		}
	}
	return tmp;
}

/**
* Take z-slice of 4d matrix turning it into 3d matrix.
* squeeze matrix to 3d getting all the data when the z dimension is at index p_z
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix3D<T> Matrix4D<T>::zSlice(int p_z) const {
	int w, x, y;
	Matrix3D<T> tmp(this->size_w, this->size_x, this->size_y);
	// !!! tmp.name = this->name + "_slice";
	for (w = 0; w < size_w; w++) {
		for (x = 0; x < size_x; x++) {
			for (y = 0; y < size_y; y++) {
				tmp[w][x][y] = matrix_array[w][x][y][p_z];
			}
		}
	}
	return tmp;
}

/**
* Take wx-slice of 4d matrix turning it into a 2d matrix.
* \param p_w - index at which to slice w dimension
* \param p_x - index at which to slice x dimension
*/
template<class T>
Matrix2D<T> Matrix4D<T>::wxSlice(int p_w, int p_x) const {
	int y, z;
	Matrix2D<T> tmp(this->size_y, this->size_z);
	// !!! tmp.name = this->name + "_xySlice";
	for (y = 0; y < size_y; y++) {
		for (z = 0; z < size_z; z++) {
			tmp[y][z] = matrix_array[p_w][p_x][y][z];
		}
	}
	return tmp;
}

/**
* Take wy-slice of 4d matrix turning it into a 2d matrix.
* \param p_w - index at which to slice w dimension
* \param p_y - index at which to slice y dimension
*/
template<class T>
Matrix2D<T> Matrix4D<T>::wySlice(int p_w, int p_y) const {
	int x, z;
	Matrix2D<T> tmp(this->size_x, this->size_z);
	// !!! tmp.name = this->name + "_xySlice";
	for (x = 0; x < size_x; x++) {
		for (z = 0; z < size_z; z++) {
			tmp[x][z] = matrix_array[p_w][x][p_y][z];
		}
	}
	return tmp;
}

/**
* Take wz-slice of 4d matrix turning it into a 2d matrix.
* \param p_w - index at which to slice w dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix2D<T> Matrix4D<T>::wzSlice(int p_w, int p_z) const {
	int x, y;
	Matrix2D<T> tmp(this->size_x, this->size_y);
	// !!! tmp.name = this->name + "_xySlice";
	for (x = 0; x < size_x; x++) {
		for (y = 0; y < size_y; y++) {
			tmp[x][y] = matrix_array[p_w][x][y][p_z];
		}
	}
	return tmp;
}

/**
* Take xy-slice of 4d matrix turning it into a 2d matrix.
* \param p_x - index at which to slice x dimension
* \param p_y - index at which to slice y dimension
*/
template<class T>
Matrix2D<T> Matrix4D<T>::xySlice(int p_x, int p_y) const {
	int w, z;
	Matrix2D<T> tmp(this->size_w, this->size_z);
	// !!! tmp.name = this->name + "_xySlice";
	for (w = 0; w < size_w; w++) {
		for (z = 0; z < size_z; z++) {
			tmp[w][z] = matrix_array[w][p_x][p_y][z];
		}
	}
	return tmp;
}

/**
* Take xz-slice of 4d matrix turning it into a 2d matrix.
* \param p_x - index at which to slice x dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix2D<T> Matrix4D<T>::xzSlice(int p_x, int p_z) const {
	int w, y;
	Matrix2D<T> tmp(this->size_w, this->size_y);
	// !!! tmp.name = this->name + "_xySlice";
	for (w = 0; w < size_w; w++) {
		for (y = 0; y < size_y; y++) {
			tmp[w][y] = matrix_array[w][p_x][y][p_z];
		}
	}
	return tmp;
}

/**
* Take yz-slice of 4d matrix turning it into a 2d matrix.
* \param p_y - index at which to slice y dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix2D<T> Matrix4D<T>::yzSlice(int p_y, int p_z) const {
	int w, x;
	Matrix2D<T> tmp(this->size_w, this->size_x);
	// !!! tmp.name = this->name + "_xySlice";
	for (w = 0; w < size_w; w++) {
		for (x = 0; x < size_x; x++) {
			tmp[w][x] = matrix_array[w][x][p_y][p_z];
		}
	}
	return tmp;
}

/**
* Take wxy-slice of 4d matrix turning it into a 1d matrix.
* \param p_w - index at which to slice w dimension
* \param p_x - index at which to slice x dimension
* \param p_y - index at which to slice y dimension
*/
template<class T>
Matrix1D<T> Matrix4D<T>::wxySlice(int p_w, int p_x, int p_y) const {
	int z;
	Matrix1D<T> tmp(this->size_z);
	// !!! tmp.name = this->name + "_xySlice";
	for (z = 0; z < size_z; z++) {
		tmp[z] = matrix_array[p_w][p_x][p_y][z];
	}
	return tmp;
}

/**
* Take wxz-slice of 4d matrix turning it into a 1d matrix.
* \param p_w - index at which to slice w dimension
* \param p_x - index at which to slice x dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix1D<T> Matrix4D<T>::wxzSlice(int p_w, int p_x, int p_z) const {
	int y;
	Matrix1D<T> tmp(this->size_y);
	// !!! tmp.name = this->name + "_xySlice";
	for (y = 0; y < size_y; y++) {
		tmp[y] = matrix_array[p_w][p_x][y][p_z];
	}
	return tmp;
}

/**
* Take wyz-slice of 4d matrix turning it into a 1d matrix.
* \param p_w - index at which to slice w dimension
* \param p_y - index at which to slice y dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix1D<T> Matrix4D<T>::wyzSlice(int p_w, int p_y, int p_z) const {
	int x;
	Matrix1D<T> tmp(this->size_x);
	// !!! tmp.name = this->name + "_xySlice";
	for (x = 0; x < size_x; x++) {
		tmp[x] = matrix_array[p_w][x][p_y][p_z];
	}
	return tmp;
}

/**
* Take xyz-slice of 4d matrix turning it into a 1d matrix.
* \param p_x - index at which to slice x dimension
* \param p_y - index at which to slice y dimension
* \param p_z - index at which to slice z dimension
*/
template<class T>
Matrix1D<T> Matrix4D<T>::xyzSlice(int p_x, int p_y, int p_z) const {
	int w;
	Matrix1D<T> tmp(this->size_w);
	// !!! tmp.name = this->name + "_xySlice";
	for (w = 0; w < size_w; w++) {
		tmp[w] = matrix_array[w][p_x][p_y][p_z];
	}
	return tmp;
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
	int x, y, z, id;
	// !!!	if (z_size > 0) { // means 3d
	if (z_size > 1) { // means 3d
		for (x = -n_of_diags; x <= n_of_diags; x++) {
			for (y = -n_of_diags; y <= n_of_diags; y++) {
				for (z = -n_of_diags; z <= n_of_diags; z++) {
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
		for (x = -n_of_diags; x <= n_of_diags; x++) {
			for (y = -n_of_diags; y <= n_of_diags; y++) {
				// calculating a diagonal number (id)
				id = this->index1d(x, y);
				// allocating memory for each diagonal
				(*this)[id] = Matrix1D<double>(this->total_size);
				(*this)[id] = 0;
			}
		}
		// !!!	} else if (x_size > 0) {
	} else if (x_size > 1) {
		for (x = -n_of_diags; x <= n_of_diags; x++) {
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
void CalculationMatrix::writeToFile(string filename) {
	int in;

	ofstream output(filename.c_str());
	//if (output==NULL && (filename.find("Debug") == string::npos)) {
	if (!output.is_open() && (filename.find("Debug") == string::npos)) {
		printf("FILE: Unable to output file: %s", filename.c_str());
		exit(EXIT_FAILURE);
	}
	output << "VARIABLES = \"";
	for (DiagMatrix::iterator it = (*this).begin(); it != (*this).end(); it++) {
		output << "\"" << it->first << "\", ";
	}
	output << endl;
	output << "ZONE T=\"" << "\", I=" << total_size << endl;
	output.setf(ios_base::scientific, ios_base::floatfield);
	for (in = 0; in < this->total_size; in++) {
		for (DiagMatrix::iterator it = (*this).begin(); it != (*this).end(); it++) {
			output << "\t" << it->second[in];
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

