/**
 *  \file UpdatableMatrix.cpp
 *
 *  These can act just like matrices, but have the ability to be updated from ini-files at any point in time.
 *
 *  Templetes are used so we don't repeat exactly the same code 4 times - for 1D, 2D, 3D, and 4D matrices
 *
 * \brief Same functionality as matrices found in Matrix.h but can also be updated from ini-files
 */

#include "UpdatableMatrix.h"

#define UPDATE_EXT ".lst"

using namespace std;



// Functions that make using temlates possible:
// We need to create functions with the same number of arguments for 1D-4D so they can be called in a similar manner

// Allocating memory

/**  Allocates memory using the allocate memory function for Matrix1D
*/
void inline MatrixAllocateMemory(Matrix1D<double> &M, const Matrix1D<double>& Q1, const Matrix1D<double>& Q2, const Matrix1D<double>& Q3, const Matrix1D<double>& Q4) {
	M.Matrix1D<double>::AllocateMemory(Q1.size_q1);
}
/**  Allocates memory using the allocate memory function for Matrix2D
*/
void inline MatrixAllocateMemory(Matrix2D<double> &M, const Matrix2D<double>& Q1, const Matrix2D<double>& Q2, const Matrix2D<double>& Q3, const Matrix2D<double>& Q4) {
	M.Matrix2D<double>::AllocateMemory(Q1.size_q1, Q2.size_q2);
}
/**  Allocates memory using the allocate memory function for Matrix3D
*/
void inline MatrixAllocateMemory(Matrix3D<double> &M, const Matrix3D<double>& Q1, const Matrix3D<double>& Q2, const Matrix3D<double>& Q3, const Matrix3D<double>& Q4) {
	M.Matrix3D<double>::AllocateMemory(Q1.size_q1, Q2.size_q2, Q3.size_q3);
}
/**  Allocates memory using the allocate memory function for Matrix4D
*/
void inline MatrixAllocateMemory(Matrix4D<double> &M, const Matrix4D<double>& Q1, const Matrix4D<double>& Q2, const Matrix4D<double>& Q3, const Matrix4D<double>& Q4) {
	M.Matrix4D<double>::AllocateMemory(Q1.size_w, Q2.size_x, Q3.size_y, Q4.size_z);
}




/*
void MatrixAllocateMemory(Matrix1D<double> &M, int size_Q1, int size_Q2, int size_Q3, int size_Q4) {
	M.Matrix1D<double>::AllocateMemory(size_Q1);
}
void MatrixAllocateMemory(Matrix2D<double> &M, int size_Q1, int size_Q2, int size_Q3, int size_Q4) {
	M.Matrix2D<double>::AllocateMemory(size_Q1, size_Q2);
}
void MatrixAllocateMemory(Matrix3D<double> &M, int size_Q1, int size_Q2, int size_Q3, int size_Q4) {
	M.Matrix3D<double>::AllocateMemory(size_Q1, size_Q2, size_Q3);
}
void MatrixAllocateMemory(Matrix4D<double> &M, int size_Q1, int size_Q2, int size_Q3, int size_Q4) {
	M.Matrix4D<double>::AllocateMemory(size_Q1, size_Q2, size_Q3, size_Q4);
}*/




// Reading from a file

/**
* Reading from a file using Matrix1D::readFromFile(string filename, const Matrix1D< T > grid_x)
*
* Only using Q1
*/
void MatrixReadFromFile(Matrix1D<double> &M, string data_filename, const Matrix1D<double>& Q1, const Matrix1D<double>& Q2, const Matrix1D<double>& Q3, const Matrix1D<double>& Q4) {
	if (data_filename.substr(data_filename.length() - 4 ,4) == ".mat")
		M.Matrix1D<double>::readFromMatlabFile(data_filename, Q1);
	else
		M.Matrix1D<double>::readFromFile(data_filename, Q1);
}
/**  Reading from a file using Matrix2D::readFromFile(string filename, const Matrix2D< T > grid_x,const Matrix2D< T > grid_y)
*
* Only using Q1, Q2
*/
void MatrixReadFromFile(Matrix2D<double> &M, string data_filename, const Matrix2D<double>& Q1, const Matrix2D<double>& Q2, const Matrix2D<double>& Q3, const Matrix2D<double>& Q4) {
	if (data_filename.substr(data_filename.length() - 4 ,4) == ".mat")
		M.Matrix2D<double>::readFromMatlabFile(data_filename, Q1,Q2);
	else
		M.Matrix2D<double>::readFromFile(data_filename, Q1, Q2);
}
/**  Reading from a file using Matrix3D::readFromFile(string filename, const Matrix3D< T > grid_x,const Matrix3D< T > grid_y, const Matrix3D< T > grid_z)
*
* Only using Q1, Q2, Q3
*/
void MatrixReadFromFile(Matrix3D<double> &M, string data_filename, const Matrix3D<double>& Q1, const Matrix3D<double>& Q2, const Matrix3D<double>& Q3, const Matrix3D<double>& Q4) {
	if (data_filename.substr(data_filename.length() - 4 ,4) == ".mat")
		M.Matrix3D<double>::readFromMatlabFile(data_filename, Q1, Q2, Q3);
	else if (data_filename.substr(data_filename.length() - 4 ,4) == ".plt")
		M.Matrix3D<double>::readFromFile(data_filename, Q1, Q2, Q3);
    else if (data_filename.substr(data_filename.length() - 5 ,5) == ".pltb")
        M.Matrix3D<double>::readFromBinaryFile(data_filename);
    else {
        printf("Unknown data format");
        exit(EXIT_FAILURE);
    }
}
/**  Reading from a file using Matrix4D::readFromFile(string filename, const Matrix4D< T > grid_w, const Matrix4D< T > grid_x,const Matrix4D< T > grid_y, const Matrix4D< T > grid_z)
*
* Using Q1, Q2, Q3, Q4
*/
void MatrixReadFromFile(Matrix4D<double> &M, string data_filename, const Matrix4D<double>& Q1, const Matrix4D<double>& Q2, const Matrix4D<double>& Q3, const Matrix4D<double>& Q4) {
	if (data_filename.substr(data_filename.length() - 4 ,4) == ".mat")
		M.Matrix4D<double>::readFromMatlabFile(data_filename, Q1, Q2, Q3, Q4);
	else if (data_filename.substr(data_filename.length() - 4 ,4) == ".plt")
		M.Matrix4D<double>::readFromFile(data_filename, Q1, Q2, Q3, Q4);
    else if (data_filename.substr(data_filename.length() - 5 ,5) == ".pltb")
        M.Matrix4D<double>::readFromBinaryFile(data_filename);
    else {
        printf("Unknown data format");
        exit(EXIT_FAILURE);
    }
}

/**
 * UpdatableMatrix has 2 arrays: original and current (=parent).
 * We need to allocate memory for both of them when the memory is allocated for the aprent array.
 * So we use polymorphism concept to override parent method AllocateMemory and to allocate memory for both arrays
 */
/*template <typename MatrixND>
void UpdatableMatrix<MatrixND>::AllocateMemory( int size_Q1, int size_Q2, int size_Q3, int size_Q4 ) {
	MatrixAllocateMemory(this->original_arr, size_Q1, size_Q2, size_Q3, size_Q4);
	MatrixAllocateMemory(*this, size_Q1, size_Q2, size_Q3, size_Q4);
}
template <typename MatrixND>
void UpdatableMatrix<MatrixND>::AllocateMemory( int size_Q1, int size_Q2, int size_Q3 ) {
	MatrixAllocateMemory(this->original_arr, size_Q1, size_Q2, size_Q3, -1);
	MatrixAllocateMemory(*this, size_Q1, size_Q2, size_Q3, -1);
}
template <typename MatrixND>
void UpdatableMatrix<MatrixND>::AllocateMemory( int size_Q1, int size_Q2 ) {
	MatrixAllocateMemory(this->original_arr, size_Q1, size_Q2, -1, -1);
	MatrixAllocateMemory(*this, size_Q1, size_Q2, -1, -1);
}
template <typename MatrixND>
void UpdatableMatrix<MatrixND>::AllocateMemory( int size_Q1 ) {
	MatrixAllocateMemory(this->original_arr, size_Q1, -1, -1, -1);
	MatrixAllocateMemory(*this, size_Q1, -1, -1, -1);
}*/


/**
 * Similar to AllocateMomory, we need to update both arrays: original and current (parent) when we read values from a file
 */
/*template <typename MatrixND>
void UpdatableMatrix<MatrixND>::readFromFile(string data_filename, const MatrixND grid_Q1, const MatrixND grid_Q2, const MatrixND grid_Q3, const MatrixND grid_Q4) {
	MatrixReadFromFile(*this, data_filename, grid_Q1, grid_Q2, grid_Q3, grid_Q4);
	original_arr = *this;
}
template <typename MatrixND>
void UpdatableMatrix<MatrixND>::readFromFile(string data_filename, const MatrixND grid_Q1, const MatrixND grid_Q2, const MatrixND grid_Q3) {
	MatrixND unused_matrix;
	MatrixReadFromFile(*this, data_filename, grid_Q1, grid_Q2, grid_Q3, unused_matrix);
	original_arr = *this;
}
template <typename MatrixND>
void UpdatableMatrix<MatrixND>::readFromFile(string data_filename, const MatrixND grid_Q1, const MatrixND grid_Q2) {
	MatrixND unused_matrix;
	MatrixReadFromFile(*this, data_filename, grid_Q1, grid_Q2, unused_matrix, unused_matrix);
	original_arr = *this;
}
template <typename MatrixND>
void UpdatableMatrix<MatrixND>::readFromFile(string data_filename, const MatrixND grid_Q1) {
	MatrixND unused_matrix;
	MatrixReadFromFile(*this, data_filename, grid_Q1, unused_matrix, unused_matrix, unused_matrix);
	original_arr = *this;
}*/


/**
 * Save current matrix as 'original' to be scaled/limited later
 */
template <typename MatrixND>
void UpdatableMatrix<MatrixND>::saveCurrent() {
	this->original_arr = *this;
}

/**
 * Function to limit the data on any direction (e.g. only above plasmapause location, or only on day/night side)
 *
 * Sets all values to 0 which are not in the range of "from" to "to"
 *
 * \param M - The matrix that will be updated(limited)
 */
void MatrixLimit(UpdatableMatrix< Matrix1D<double> > &M, const Matrix1D<double> &Q1, const Matrix1D<double> &Q2, const Matrix1D<double> &Q3, const Matrix1D<double> &Q4,
		double Q1_from, double Q1_to, double Q2_from, double Q2_to, double Q3_from, double Q3_to, double Q4_from, double Q4_to) {

	// Set zeros everywhere where we shouldn't have the values (this is the limiting)
	if (Q1_from > Q1.min() || Q1_to < Q1.max()) {
		int iQ1;
		for (iQ1 = 0; iQ1 < Q1.size_q1; iQ1++) {
			if (
					Q1[iQ1] < Q1_from || Q1[iQ1] > Q1_to
				) {
					M[iQ1] = 0;
			}
		}
	}
}

/**
 * Function to limit the data on any direction (e.g. only above plasmapause location, or only on day/night side)
 *
 * Sets all values to 0 which are not in the range of "from" to "to"
 *
 * \param M - The matrix that will be updated(limited)
 */
void MatrixLimit(UpdatableMatrix<Matrix2D<double>> &M, const Matrix2D<double> &Q1, const Matrix2D<double> &Q2, const Matrix2D<double> &Q3, const Matrix2D<double> &Q4,
		double Q1_from, double Q1_to, double Q2_from, double Q2_to, double Q3_from, double Q3_to, double Q4_from, double Q4_to) {
	// Set zeros everywhere where we shouldn't have the values (this is the limiting)
	if (Q1_from > Q1.min() || Q1_to < Q1.max() || Q2_from > Q2.min() || Q2_to < Q2.max()) {
		int iQ1, iQ2;
		for (iQ1 = 0; iQ1 < Q1.size_q1; iQ1++) {
			for (iQ2 = 0; iQ2 < Q1.size_q2; iQ2++) {
				if (
						Q1[iQ1][iQ2] < Q1_from || Q1[iQ1][iQ2] > Q1_to ||
						Q2[iQ1][iQ2] < Q2_from || Q2[iQ1][iQ2] > Q2_to
					) {
						M[iQ1][iQ2] = 0;
				}
			}
		}
	}
}

/**
 * Function to limit the data on any direction (e.g. only above plasmapause location, or only on day/night side)
 *
 * Sets all values to 0 which are not in the range of "from" to "to"
 *
 * \param M - The matrix that will be updated(limited)
 */
void MatrixLimit(UpdatableMatrix< Matrix3D<double> > &M, const Matrix3D<double> &Q1, const Matrix3D<double> &Q2, const Matrix3D<double> &Q3, const Matrix3D<double> &Q4,
		double Q1_from, double Q1_to, double Q2_from, double Q2_to, double Q3_from, double Q3_to, double Q4_from, double Q4_to) {
	// Set zeros everywhere where we shouldn't have the values (this is the limiting)
	if (Q1_from > Q1.min() || Q1_to < Q1.max() || Q2_from > Q2.min() || Q2_to < Q2.max() || Q3_from > Q3.min() || Q3_to < Q3.max()) {
		int iQ1, iQ2, iQ3;
		for (iQ1 = 0; iQ1 < Q1.size_q1; iQ1++) {
			for (iQ2 = 0; iQ2 < Q1.size_q2; iQ2++) {
				for (iQ3 = 0; iQ3 < Q1.size_q3; iQ3++) {
					if (
							Q1[iQ1][iQ2][iQ3] < Q1_from || Q1[iQ1][iQ2][iQ3] > Q1_to ||
							Q2[iQ1][iQ2][iQ3] < Q2_from || Q2[iQ1][iQ2][iQ3] > Q2_to ||
							Q3[iQ1][iQ2][iQ3] < Q3_from || Q3[iQ1][iQ2][iQ3] > Q3_to
						) {
							M[iQ1][iQ2][iQ3] = 0;
					}
				}
			}
		}
	}
}

/**
 * Function to limit the data on any direction (e.g. only above plasmapause location, or only on day/night side)
 *
 * Sets all values to 0 which are not in the range of "from" to "to"
 *
 * \param M - The matrix that will be updated(limited)
 */
void MatrixLimit(UpdatableMatrix< Matrix4D<double> > &M, const Matrix4D<double> &Q1, const Matrix4D<double> &Q2, const Matrix4D<double> &Q3, const Matrix4D<double> &Q4,
		double Q1_from, double Q1_to, double Q2_from, double Q2_to, double Q3_from, double Q3_to, double Q4_from, double Q4_to) {
	// Set zeros everywhere where we shouldn't have the values (this is the limiting)
	if (Q1_from > Q1.min() || Q1_to < Q1.max() || Q2_from > Q2.min() || Q2_to < Q2.max() || Q3_from > Q3.min() || Q3_to < Q3.max() || Q4_from > Q4.min() || Q4_to < Q4.max()) {
		int iQ1, iQ2, iQ3, iQ4;
		for (iQ1 = 0; iQ1 < Q1.size_w; iQ1++) {
			for (iQ2 = 0; iQ2 < Q1.size_x; iQ2++) {
				for (iQ3 = 0; iQ3 < Q1.size_y; iQ3++) {
					for (iQ4 = 0; iQ4 < Q1.size_z; iQ4++) {
						if (
								Q1[iQ1][iQ2][iQ3][iQ4] < Q1_from || Q1[iQ1][iQ2][iQ3][iQ4] > Q1_to ||
								Q2[iQ1][iQ2][iQ3][iQ4] < Q2_from || Q2[iQ1][iQ2][iQ3][iQ4] > Q2_to ||
								Q3[iQ1][iQ2][iQ3][iQ4] < Q3_from || Q3[iQ1][iQ2][iQ3][iQ4] > Q3_to ||
								Q4[iQ1][iQ2][iQ3][iQ4] < Q4_from || Q4[iQ1][iQ2][iQ3][iQ4] > Q4_to
							) {
								M[iQ1][iQ2][iQ3][iQ4] = 0;
						}
					}
				}
			}
		}
	}
}

/**
 * Read UpdatableMatrix rules from ini-file
 *
 * Gets parameters from the first line using readFromString()
 *
 * File will either be data file (such as .plt which will store data immediately) or update file (such as .lst which will update matrix values at corresponding timestep from main function)
 *
 * @return True if success, False if failure
 */
template <typename MatrixND>
bool UpdatableMatrix<MatrixND>::readFromIniFile(string ini_filename, const MatrixND& q1, const MatrixND& q2, const MatrixND& q3, const MatrixND& q4) {
	fstream input;

	// Open the ini-file for read
	input.open(ini_filename.c_str());
	//if (input == NULL) {
	if (!input.is_open() ) {
		// if file not found - return false
		Logger::warning << ini_filename << " not found." << endl;
		return false;
	}
	Logger::message << "Reading " << ini_filename << ": " << endl;

	string file_line_string;

	// if the file is not empty
	if(!input.eof()) {
		file_line_string = "";

		// read one line from the file into a string
		getline(input, file_line_string);

		// check that the line is not empty and is not a comment line
		if (!file_line_string.empty() && file_line_string.find("#") == string::npos) {

			// get all the parameters from one line
			this->readFromString(file_line_string, q1, q2, q3, q4);

			// output - next line
			Logger::message << endl;

		}
	}

	return true;
}



/**
 * Read UpdateMatrix parameters from a string
 *
 * Will traverse the string and get the entire lines values -
 * starting with parameter q1 and getting values until reaching q4 or running out of values for the line
 *
 * If a data file (such as .plt) is selected then the contents will be stored into the matrix.
 * If a update-file (such as .lst) is selected then no matrix values will be saved, however the update_filename field will be saved.
 * Upon every timestep the main file will call update() on this matrix, and inside update will look for the update_filename field. If present it will update the matrix to the corresponding timestep
 *
 * @param file_line_string - string which contains one line of values spanning from q1 to q4
 *
 * @return True if success, False if failure
 */
template <typename MatrixND>
bool UpdatableMatrix<MatrixND>::readFromString(string file_line_string, const MatrixND &q1, const MatrixND &q2, const MatrixND &q3, const MatrixND &q4) {
	// This defines an "empty marker" - it indicates that we need to skip that parameter.
	const string empty_marker = "-";

	// Convert input string into input stream
	// It's much easier to work with streams
	stringstream file_line_stream;
	file_line_stream << file_line_string;

	// Read data filename from the stream (first parameter)
	string data_filename;
	file_line_stream >> data_filename;
	Logger::message << "	Data: " << data_filename << endl;

	// Check if what we've read is not empty (or is not an "empty marker")
	if (!data_filename.empty() && data_filename != empty_marker) {

		// The first parameter is either data-filename or update-filename, check it (simply by file extension)
		// if (data_filename.size() >= strlen(UPDATE_EXT) && data_filename.substr(data_filename.size() - strlen(UPDATE_EXT), data_filename.size()) == UPDATE_EXT) {
		if (data_filename.rfind(UPDATE_EXT) == data_filename.length() - strlen(UPDATE_EXT)) { // string::npos) {
			// This update-file will tell us when and how to update the array
			this->update_filename = data_filename;
			Logger::message << "	Update from: " << this->update_filename << endl;
			this->update_pos = 0;

		} else {
			// If it's not an update-file, no updates
			this->update_filename = empty_marker;

			// So, it's data-file, read the data-file
			// TODO: check if it's just a number, not a filename
			MatrixReadFromFile(*this, data_filename, q1, q2, q3, q4);

			// save the last name that we've read, for user-friendly output
			// this->last_data_filename = data_filename;
		}

	}

	// Check if memory was allocated for original_arr, allocate it if not
	if (!this->original_arr.initialized) {
		MatrixAllocateMemory(this->original_arr, q1, q2, q3, q4);
	}

	// Copy current array into original array
	this->original_arr = *this;

	// Save the last update timestep - it means the file will be updated after the code starts
	this->last_update_time = -1;

	// Check for the end of line
	if (!file_line_stream.eof()) {
		// If not the end - read the next parameter, which is scaling coefficient
		file_line_stream >> this->scale_string;
		Logger::message << "	Scale: " << this->scale_string << endl;
		this->scale_pos = 0;
	} else {
		// If there is nothing - no scaling
		this->scale_string  = "";
	}

	// Check for the end of line
	if (!file_line_stream.eof()) {
		// If not the end - read the next parameters, which are 'from' and 'to' boundaries for the array
		// (like above Lpp or below Lpp or something like that)
		file_line_stream >> this->Q1_from_string;
		file_line_stream >> this->Q1_to_string;
		Logger::message << "	Limits: " << this->Q1_from_string << " to " << this->Q1_to_string << endl;
		this->Q1_from_pos = 0;
		this->Q1_to_pos = 0;
	} else {
		// If there is nothing - no limiting
		this->Q1_from_string = "";
		this->Q1_to_string = "";
	}

	// Check for the end of line
	if (!file_line_stream.eof()) {
		// If not the end - read the next parameters, which are 'from' and 'to' boundaries for the array
		// (like above Lpp or belowe Lpp or something like that)
		file_line_stream >> this->Q2_from_string;
		file_line_stream >> this->Q2_to_string;
		Logger::message << "	Limits: " << this->Q2_from_string << " to " << this->Q2_to_string << endl;
		this->Q2_from_pos = 0;
		this->Q2_to_pos = 0;
	} else {
		// If there is nothing - no limiting
		this->Q2_from_string = "";
		this->Q2_to_string = "";
	}

	// Check for the end of line
	if (!file_line_stream.eof()) {
		// If not the end - read the next parameters, which are 'from' and 'to' boundaries for the array
		// (like above Lpp or belowe Lpp or something like that)
		file_line_stream >> this->Q3_from_string;
		file_line_stream >> this->Q3_to_string;
		Logger::message << "	Limits: " << this->Q3_from_string << " to " << this->Q3_to_string << endl;
		this->Q3_from_pos = 0;
		this->Q3_to_pos = 0;
	} else {
		// If there is nothing - no limiting
		this->Q3_from_string = "";
		this->Q3_to_string = "";
	}

	// Check for the end of line
	if (!file_line_stream.eof()) {
		// If not the end - read the next parameters, which are 'from' and 'to' boundaries for the array
		// (like above Lpp or belowe Lpp or something like that)
		file_line_stream >> this->Q4_from_string;
		file_line_stream >> this->Q4_to_string;
		Logger::message << "	Limits: " << this->Q4_from_string << " to " << this->Q4_to_string << endl;
		this->Q4_from_pos = 0;
		this->Q4_to_pos = 0;
	} else {
		// If there is nothing - no limiting
		this->Q4_from_string = "";
		this->Q4_to_string = "";
	}

	return true;
}


/**
 * Update the Matrix if last update time is less than current update time
 *
 * Scaling coefficients are added
 *
 * Limits are set for the coefficients for q1 through q4
 *
 * The limits are applied using MatrixLimit()
 *
 * Record which limits were applied
 *
 * This function is called at every time step from the main function. It will look to see if this matrix has the update_filename field set.
 * If so at any time where the update_time is less than current_time it will use the update_filename to load the new matrix. This pertains to .lst files mainly
 *
 * @return True if updated, False if not
 */
template <typename MatrixND>
bool UpdatableMatrix<MatrixND>::update(double current_time, const MatrixND& q1, const MatrixND& q2, const MatrixND& q3, const MatrixND& q4) {
	// Indicator if updated
	bool updated = false;

	// This is an "empty marker" - it indicates that we need to skip that parameter.
	// We can't just leave it as a space, cause multiple spaces in a row are considered as one space
	const string empty_marker = "-";

	ostringstream user_output("");

	// current computation time
	//double current_time = iteration * dt;

	// A time from the scaling, limiting, and updating file, will be used later
	double update_time;

	// Check if the original array was initialized,
	// if not - allocate memory and fill it with current array data
	if (!this->original_arr.initialized) {
		// Check if the current array was even initialized
		if (!this->initialized) {
			Logger::error << "Error - using uninitialized matrix" << endl;
			exit(EXIT_FAILURE);
		}
		MatrixAllocateMemory(this->original_arr, q1, q2, q3, q4);
		// Save current array as original array
		this->original_arr = *this;
	}

	// ////////////////////////////////////////////
	// Updating

	// update the data, if update-filename is not empty
	if (!this->update_filename.empty() && this->update_filename != empty_marker) {
		// filename for update
		string data_filename;

		// Search for current time step line (or closest to it) in the update-file,
		// The second column will be data-filename
		data_filename = GetCurrentTimeValue(this->update_filename, this->update_pos, current_time, update_time);

		// Check if data-filename is not an empty-file-marker
		if (!data_filename.empty()
				&& data_filename != empty_marker
				// && data_filename != this->last_data_filename
				) {
			if (current_time >= update_time && this->last_update_time < update_time) {
				// Update the Matrix data from the file
				user_output << data_filename;
				MatrixReadFromFile(this->original_arr, data_filename, q1, q2, q3, q4);
				this->last_update_time = current_time;
				// this->last_data_filename = data_filename;
				updated = true;


				// // ADDED FOR TESTING
				// original_arr.writeToFile(to_string(current_time) + "LUBC_Update.plt" );

			}
		}
	}

	// Copy the original array into the current array, so we can apply scaling, limiting, etc
	// The following line calls "operator=" from base (parent) Matrix2D class
	// and makes that class to be equal to the original_arr
	// MatrixND::operator = (original_arr);
	*this = original_arr;


	// ////////////////////////////////////////////
	// Scaling

	if (!this->scale_string.empty() && this->scale_string != empty_marker) {
		// the scaling coefficient
		double scale_coefficient;

		// the scaling coefficient in the string form
		string scale_coefficient_string;

		// this->scale_string contains either filename or a coefficient itself
		// 'GetCurrentTimeValue' will check it and return the coefficient in both cases
		scale_coefficient_string = GetCurrentTimeValue(this->scale_string, this->scale_pos, current_time, update_time);

		// convert from string value into double value
		stringstream(scale_coefficient_string) >> scale_coefficient;

		if (scale_coefficient != 1) {
			// Output for the user what we just did - scaled the array (array-name) by (scaling value)
			// cout << this->last_data_filename << " was scaled as " << scale_coefficient << endl;
			user_output << "*" << scale_coefficient;

			// Scale the array
			// The following line multiplies the current (base, parent) array by scale_coefficient
			// The result goes into base (parent) array
			(*this) *= scale_coefficient;
			// Another way to do the same:
			//this->MatrixND::operator*= (scale_coefficient);

			updated = true;
		}
	}

	// ////////////////////////////////////////////
	// Limiting

	// values for limiting - from, and to
	double Q1_from = -1e99, Q1_to = 1e99;
	double Q2_from = -1e99, Q2_to = 1e99;
	double Q3_from = -1e99, Q3_to = 1e99;
	double Q4_from = -1e99, Q4_to = 1e99;


	// limiting coefficient, as a string
	string limit_coefficient_string;
	bool was_limited = false;

	// check, if limiting rule is empty
	if (!Q1_from_string.empty() && Q1_from_string != empty_marker) {
		limit_coefficient_string = GetCurrentTimeValue(this->Q1_from_string, this->Q1_from_pos, current_time, update_time);

		// convert the limiting value from string to double
		stringstream(limit_coefficient_string) >> Q1_from;

		was_limited = true;
	}
	// Exactly the same, as for limiting-from: this is limiting-to
	if (!Q1_to_string.empty() && Q1_to_string != empty_marker) {
		// else - get the limiting value
		limit_coefficient_string = GetCurrentTimeValue(this->Q1_to_string, this->Q1_to_pos, current_time, update_time);

		// convert the limiting value from string to double
		stringstream(limit_coefficient_string) >> Q1_to;

		was_limited = true;
	}

	// check, if limiting rule is empty
	if (!Q2_from_string.empty() && Q2_from_string != empty_marker) {
		limit_coefficient_string = GetCurrentTimeValue(this->Q2_from_string, this->Q2_from_pos, current_time, update_time);

		// convert the limiting value from string to double
		stringstream(limit_coefficient_string) >> Q2_from;

		was_limited = true;
	}
	// Exactly the same, as for limiting-from: this is limiting-to
	if (!Q2_to_string.empty() && Q2_to_string != empty_marker) {
		// else - get the limiting value
		limit_coefficient_string = GetCurrentTimeValue(this->Q2_to_string, this->Q2_to_pos, current_time, update_time);

		// convert the limiting value from string to double
		stringstream(limit_coefficient_string) >> Q2_to;

		was_limited = true;
	}

	// check, if limiting rule is empty
	if (!Q3_from_string.empty() && Q3_from_string != empty_marker) {
		limit_coefficient_string = GetCurrentTimeValue(this->Q3_from_string, this->Q3_from_pos, current_time, update_time);

		// convert the limiting value from string to double
		stringstream(limit_coefficient_string) >> Q3_from;

		was_limited = true;
	}
	// Exactly the same, as for limiting-from: this is limiting-to
	if (!Q3_to_string.empty() && Q3_to_string != empty_marker) {
		// else - get the limiting value
		limit_coefficient_string = GetCurrentTimeValue(this->Q3_to_string, this->Q3_to_pos, current_time, update_time);

		// convert the limiting value from string to double
		stringstream(limit_coefficient_string) >> Q3_to;

		was_limited = true;
	}

	// check, if limiting rule is empty
	if (!Q4_from_string.empty() && Q4_from_string != empty_marker) {
		limit_coefficient_string = GetCurrentTimeValue(this->Q4_from_string, this->Q4_from_pos, current_time, update_time);

		// convert the limiting value from string to double
		stringstream(limit_coefficient_string) >> Q4_from;

		was_limited = true;
	}
	// Exactly the same, as for limiting-from: this is limiting-to
	if (!Q4_to_string.empty() && Q4_to_string != empty_marker) {
		// else - get the limiting value
		limit_coefficient_string = GetCurrentTimeValue(this->Q4_to_string, this->Q4_to_pos, current_time, update_time);

		// convert the limiting value from string to double
		stringstream(limit_coefficient_string) >> Q4_to;

		was_limited = true;
	}

	// Set zeros everywhere where we shouldn't have the values (this is the limiting)
	MatrixLimit(*this, q1, q2, q3, q4, Q1_from, Q1_to, Q2_from, Q2_to, Q3_from, Q3_to, Q4_from, Q4_to);

	// Tell user what we've done - what we have limited, from where and to where
	if (was_limited) {
		user_output << " ";
		//if (Q1_from > -1e99) user_output << " from " << Q1_from;
		//if (Q1_to < 1e99) user_output << " to " << Q1_to;
		if (Q1_from > -1e99 || Q1_to < 1e99) {
			user_output << "[" << q1.name << ": " << Q1_from << "-" << Q1_to << "] ";
		}

		//if (Q2_from > -1e99) user_output << " from " << Q2_from;
		//if (Q2_to < 1e99) user_output << " to " << Q2_to;
		if (Q2_from > -1e99 || Q2_to < 1e99) {
			user_output << "[" << q2.name << ": " << Q2_from << "-" << Q2_to << "] ";
		}

		//if (Q3_from > -1e99) user_output << " from " << Q3_from;
		//if (Q3_to < 1e99) user_output << " to " << Q3_to;
		if (Q3_from > -1e99 || Q3_to < 1e99) {
			user_output << "[" << q3.name << ": " << Q3_from << "-" << Q3_to << "] ";
		}

		//if (Q4_from > -1e99) user_output << " from " << Q4_from;
		//if (Q4_to < 1e99) user_output << " to " << Q4_to;
		if (Q4_from > -1e99 || Q4_to < 1e99)  {
			user_output << "[" << q4.name << ": " << Q4_from << "-" << Q4_to << "] ";
		}
	}

	if (user_output.str().length() > 0) {
		Logger::message << this->name << ": " << user_output.str() << endl;
	}

	return updated;

}


/**
 * Read UpdatableListMatrix rules from ini-file
 * UpdatableListMatrix is a list (vector) or UpdatableMatrices
 * They are combined together to form one matrix (e.g. diffusion coefficients for different waves)
 *
 * The inputted file can be of any format. If it is a data file (such as .plt) the data will be stored into the matrix.
 * If the file is an update file then the matrix will only be updated on the corresponding timestep from the main file.
 * The inputted file may be a .tab file that has a combination of data and update files.
 * This function goes through the inputted file and stores the contents individually into a vector of UpdatableMatrix, which in turn will be evaluated from the calls to the UpdatableMatrix functions
 *
 * @param ini_filename - ini filename with appropriate structure
 * @param Q1 - corresponding grid coordinate
 * @param Q2 - corresponding grid coordinate
 * @param Q3 - corresponding grid coordinate
 * @param Q4 - corresponding grid coordinate
 * @return True if success, False if failure
 */
template <typename MatrixND>
bool UpdatableListMatrix<MatrixND>::readFromIniFile(string ini_filename, const MatrixND &Q1, const MatrixND &Q2, const MatrixND &Q3, const MatrixND &Q4) {
	fstream input;

	// Open the ini-file for read
	input.open(ini_filename.c_str());
	//if (input == NULL) {
	if (!input.is_open()) {
		// if file not found - diff coeff is zero
		Logger::warning << ini_filename << " not found." << endl;
		return false;
	}
	Logger::message << "Reading " << ini_filename << ": " << endl;

	string file_line_string;// = "";

	// read while the ini-file is not empty
	while(!input.eof()) {
		file_line_string = "";

		// Read one line from the file into a string
		// we do it this way cause we don't know how many values are in one line,
		// since some of the values are optional
		getline(input, file_line_string);

		// check that the line is not a comment line
		if (!file_line_string.empty() && file_line_string.find("#") == string::npos) {

			// Create new UpdatableMatrix
			UpdatableMatrix<MatrixND> newUpdatableMatrix;
			MatrixAllocateMemory(newUpdatableMatrix, Q1, Q2, Q3, Q4);
			MatrixAllocateMemory(newUpdatableMatrix.original_arr, Q1, Q2, Q3, Q4);

			// read all the UpdatableMatrix's parameters from the line
			newUpdatableMatrix.readFromString(file_line_string, Q1, Q2, Q3, Q4);

			// add the matrix to the list
			matricesList.push_back(newUpdatableMatrix);

			// move output cursor to the next line
			Logger::message << endl;

		}
	}

	return true;

}



/**
 * Update the UpdatableListMatrix - call update() for each UpdatableMatrix in the list, then sum all the UpdatableMatrices together
 *
 * @param iteration - current iteration number
 * @param dt - time step
 * @param Q1 - corresponding grid coordinate
 * @param Q2 - corresponding grid coordinate
 * @param Q3 - corresponding grid coordinate
 * @param Q4 - corresponding grid coordinate
 */
template <typename MatrixND>
void UpdatableListMatrix<MatrixND>::update(double current_time, const MatrixND& Q1, const MatrixND& Q2, const MatrixND& Q3, const MatrixND& Q4) {
	unsigned int d_it;

	// check if there are any coefficients in the list at all
	if (matricesList.size() > 0) {
		// if there are,

		// set base (parent) class to be zero
		MatrixND::operator = (0);

		// update each matrix in the list
		for (d_it = 0; d_it < matricesList.size(); d_it++)
			matricesList[d_it].update(current_time, Q1, Q2, Q3, Q4);


		// combine the UpdatableMatrices from the list together into one final updated matrix
		for (d_it = 0; d_it < matricesList.size(); d_it++)
			// Applying "operator+=" of the base class to each of the matrices
			// this will add (sum) all the matrices together one by one
			(*this) += matricesList[d_it];
	} else {
		// if there is nothing in the list - leave everything the way it is
	}
}


/**
 * Search for current time-step in an update-file and return the corresponding value
 *
 * @param filename - filename to check. Also can be just a value - then we just return this value!
 * @param current_time - the time we're searching for
 * @param update_time - this is a return value! It is the time we've found!
 * @return
 */
string GetCurrentTimeValue(const string &filename, long &pos, const double &current_time, double &update_time) {

	// check - if the input is not a filename, but just a number - in that case just return it
	if (is_number(filename)) {
		return filename;
	}

	// We are going to read the file and search for the current time step

	ifstream input;
	string return_value = "";
	long current_line = 0; // current line of the file

	// open the file
	input.open(filename.c_str(), ifstream::in);
	//if (input == NULL || input.eof()) {
	if (!input.is_open() || input.eof()) {
		Logger::error << "Error reading file " << filename.c_str() << endl;
		exit(EXIT_FAILURE);
	}

	// //////////////////////////////////////////////////////
	// Looking for the correct line in the file
	// - the one which is closest to the current time step in the past, not in the future!
	string string_tmp;
	double time_tmp;
	double last_time = -1;

	// check if we position is > 0. Then we need to skip lines
	if (pos > 0)
	{
		while (!input.eof() && (current_line < pos))
		{
			input.ignore(999, '\n'); // This is fastest way to skip lines that I found
			current_line++;
		}

	}

	while (!input.eof()) {
		// read time
		input >> string_tmp;

		// check if it's not a comment or something else not important
		if (string_tmp.find("#") == string::npos && string_tmp.find("ZONE") == string::npos && string_tmp.find("VAR") == string::npos) {

			stringstream(string_tmp) >> time_tmp;

			// check if the time_tmp is ahead of current time
			if (time_tmp <= current_time) {
				// If we still behind the current time - save the values
				// save current time into last_time
				last_time = time_tmp;
				// save current value for later use, in case this will be the last value we read
				input >> return_value;

				current_line++;
				pos = current_line; // save the last line position
			} else {
				// We've just passed current_time
				// Whatever we've saved before - is what we need, just exit
				break;
			}

		}
		// read to the end of the line
		input.ignore(9999, '\n');

	}

	if (return_value == "") {
		Logger::error << std::setprecision(15) << "ERROR: Scale factor for " << current_time << " wasn't found in " << filename << endl;
		exit(EXIT_FAILURE);
	}

	// Use the last time we've read as update_time
	update_time = last_time;
	// save the position for the future search
	input.close();

	// Move position one line back in case we have to use the scaling again
	pos--;
	if (pos < 0)
		pos = 0; // Make sure that pos is not negative

	return return_value;

}

/*
* Overloaded GetCurrentTimeValue to maintain possible compatability
*/
string GetCurrentTimeValue(const string &filename, const double &current_time, double &update_time) {
	long pos = 0;
	return GetCurrentTimeValue(filename, pos, current_time, update_time);
}

/**
 * Check if a string is number
 *
 * Note: Will accept any string (even badly formatted ones) consisting entirely of the following: Numbers, . , e , - , +
 *
 * For example 1.2.3-e.4 will be accepted
 * @param s - the string we need to check
 * @return True if it's a number, false if it's not
 */
bool is_number(const std::string& s) {
	// Loop through each character and see if it's one of the allowed characters
	// This is an example of s: "-1.123000e+15"
    string::const_iterator it = s.begin();
    while (it != s.end() && (
    		isdigit(*it) ||
    		(*it == '.') ||
    		(*it == 'e') ||
    		(*it == '+') ||
    		(*it == '-')
    		)) ++it;
    return !s.empty() && it == s.end();
}

/**
* Overloaded = operator for any size Updatable Matrix.
* Copies every element from M and returns current matrix
*/
template <typename MatrixND>
inline MatrixND& UpdatableMatrix<MatrixND>::operator= (const MatrixND &M) {
	return MatrixND::operator= (M);
}

/**
* Overloaded = operator for any size Updatable Matrix.
* Sets every elemenet to Val and returns current matrix
*/
template <typename MatrixND>
inline MatrixND& UpdatableMatrix<MatrixND>::operator= (const double Val) {
	return MatrixND::operator= (Val);
}

/**
* Overloaded = operator for any size Updatable List Matrix.
* Copies every element from M and returns current matrix
*/
template <typename MatrixND>
inline MatrixND& UpdatableListMatrix<MatrixND>::operator= (const MatrixND &M) {
	return MatrixND::operator= (M);
}

/**
* Overloaded = operator for any size Updatable List Matrix.
* Sets every elemenet to Val and returns current matrix
*/
template <typename MatrixND>
inline MatrixND& UpdatableListMatrix<MatrixND>::operator= (const double Val) {
	return MatrixND::operator= (Val);
}

// ///////////////////////////////////////////////////////////
// Implementations - Enumerates all possibilities for MatrixND
// ///////////////////////////////////////////////////////////


template class UpdatableMatrix< Matrix1D<double> >;
template class UpdatableMatrix< Matrix2D<double> >;
template class UpdatableMatrix< Matrix3D<double> >;
template class UpdatableMatrix< Matrix4D<double> >;

template class UpdatableListMatrix< Matrix1D<double> >;
template class UpdatableListMatrix< Matrix2D<double> >;
template class UpdatableListMatrix< Matrix3D<double> >;
template class UpdatableListMatrix< Matrix4D<double> >;
