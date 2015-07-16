/**
 * \file UpdatableMatrix.h
 * \brief Same functionality as matrices found in Matrix.h but can also be updated from ini-files
 */

#ifndef UPDATABLEMATRIX_H_
#define UPDATABLEMATRIX_H_

#include <string>
#include <istream>
#include <sstream>
#include <vector>

#include "Matrix.h"
#include "Logger.h"

/** @class MatrixND
 * 
 * @brief A normal matrix that can be defined to have 1, 2, 3, or 4 dimensions which is templated for UpdatableMatrix and UpdatableListMatrix
 */


/**
 * Updatable matrix - 1,2,3,4 dimensions.
 *
 * It's a template because it's easier to create, maintain, and debug one class, instead of 4 classes with the same functionality
 *
 * It's just like a normal matrix (and can be used just like one),
 * but it can also be updated according to rules from an ini-file
 *
 * TODO: delete this class, all it's functionality is inside UpdatableListMatrix
 *
 * @brief A matrix that can be created in 1, 2, 3, or 4 dimensions with the ability to be updated
 */
template <typename MatrixND>
class UpdatableMatrix : public MatrixND { // current array is the parent array - whenever we use UpdatableMatrix as an array, the parent class is used
private:
	// updated matrix
	// Matrix3D<double> updated_arr;

	//
	// string last_data_filename; /// data-filename
	string scale_string; /// Filename for matrix scaling
	string update_filename; /// Filename for Matrix data updating (e.g. to load a completely different data)
	string Q1_from_string; /// Filename for limiting Q1 (usually it's MLT, lower limit)
	string Q1_to_string; /// Filename for limiting Q1 (usually it's MLT, upper limit)
	string Q2_from_string; /// Filename for limiting Q2 (usually it's the radial distance range, e.g. above or below plasmapause)
	string Q2_to_string; /// Filename for limiting Q2 (usually it's the radial distance range, e.g. above or below plasmapause)
	string Q3_from_string; /// Filename for limiting Q3 (usually it's V, lower limit)
	string Q3_to_string; /// Filename for limiting Q3 (usually it's V, upper limit)
	string Q4_from_string; /// Filename for limiting Q4 (usually it's K, lower limit)
	string Q4_to_string; /// Filename for limiting Q4 (usually it's K, upper limit)

	double last_update_time; /// Indicates when the Matrix was updated last time

public:
	// original matrix
	MatrixND original_arr;

	// to simplify things
	inline MatrixND& operator= (const MatrixND &M);
	inline MatrixND& operator= (const double Val);

	void saveCurrent(); // Save current matrix as original to be scaled/limited later

	//void AllocateMemory( int size_Q1 );
	//void AllocateMemory( int size_Q1, int size_Q2 );
	//void AllocateMemory( int size_Q1, int size_Q2, int size_Q3 );
	//void AllocateMemory( int size_Q1, int size_Q2, int size_Q3, int size_Q4 );

	//void readFromFile(string filename, const MatrixND grid_w); ///< Load matrix to a file
	//void readFromFile(string filename, const MatrixND grid_w, const MatrixND grid_x); ///< Load matrix to a file
	//void readFromFile(string filename, const MatrixND grid_w, const MatrixND grid_x, const MatrixND grid_y); ///< Load matrix to a file
	//void readFromFile(string filename, const MatrixND grid_w, const MatrixND grid_x, const MatrixND grid_y, const MatrixND grid_z); ///< Load matrix to a file

	bool readFromIniFile(string ini_filename, MatrixND Q1, MatrixND Q2, MatrixND Q3, MatrixND Q4 = MatrixND());
	bool readFromString(string file_line_string, MatrixND &Q1, MatrixND &Q2, MatrixND &Q3, MatrixND &Q4);
	bool update(double time, MatrixND Q1, MatrixND Q2, MatrixND Q3, MatrixND Q4 = MatrixND());

};

/**
 * Updatable list of matrix - 1,2,3,4 dimensions.
 *
 * It's a template because it's easier to create, maintain, and debug one class, instead of 4 classes with the same functionality
 *
 * It's just like a normal matrix (and can be used just like one),
 * but it can also be updated according to rules from an ini-file
 *
 * The ini-files have specify several matrices that are added to each other,
 * e.g. diffusion coefficients for different waves to get one diffusion coefficient
 * 
 * @brief A matrix that can be created in 1, 2, 3, or 4 dimensions with the ability to be updated
 */
template <typename MatrixND>
class UpdatableListMatrix : public MatrixND {
private:
	/// These are the matrices that need to be combined to get the final, updated matrix
	vector < UpdatableMatrix<MatrixND> > matricesList;

public:

	// to simplify things
	inline MatrixND& operator= (const MatrixND &M);
	inline MatrixND& operator= (const double Val);

	bool readFromIniFile(string ini_filename, MatrixND &Q1, MatrixND &Q2, MatrixND &Q3, MatrixND &Q4);
	void update(double time, MatrixND Q1, MatrixND Q2, MatrixND Q3, MatrixND Q4 = MatrixND());

};

// couple of helpful functions
string GetCurrentTimeValue(string filename, double current_time, double &update_time);
bool is_number(const std::string& s);
/// FUNCTION NOT IMPLEMENTED
double stringToValue(string string_value, double current_time);

#endif /* UPDATABLEMATRIX_H_ */
