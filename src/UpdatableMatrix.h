/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * \file UpdatableMatrix.h
 * \brief Same functionality as matrices found in Matrix.h but can also be updated from ini-files
 */

#pragma once

#include <string>
#include <istream>
#include <sstream>
#include <vector>
#include <queue>
#include <iomanip> // std::setprecision

#include "Matrix.h"
#include "Logger.h"

template <typename T>
struct QueueElement {
    double time;
    T value;
};

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
 * The readFromIniFile function is used to load data into an UpdatableMatrix by storing the contents of a .tab file.
 * The .tab file will either list a .plt file to load or a .lst file which in turn has both time steps and .plt files.
 * The update function is called at every time step in the main function and will execute any updates that can be found from these files that match the designated timestep.
 *
 * @brief A matrix that can be created in 1, 2, 3, or 4 dimensions with the ability to be updated
 */
template <typename MatrixND>
class UpdatableMatrix : public MatrixND { // current array is the parent array - whenever we use UpdatableMatrix as an array, the parent class is used
private:
	// updated matrix
	// Matrix3D<double> updated_arr;

	//
	// std::string last_data_filename; /// data-filename
	std::string scale_string; /// Filename for matrix scaling
	std::string update_filename; /// Filename for Matrix data updating (e.g. to load a completely different data)
	std::string Q1_from_string; /// Filename for limiting Q1 (usually it's MLT, lower limit)
	std::string Q1_to_string; /// Filename for limiting Q1 (usually it's MLT, upper limit)
	std::string Q2_from_string; /// Filename for limiting Q2 (usually it's the radial distance range, e.g. above or below plasmapause)
	std::string Q2_to_string; /// Filename for limiting Q2 (usually it's the radial distance range, e.g. above or below plasmapause)
	std::string Q3_from_string; /// Filename for limiting Q3 (usually it's V, lower limit)
	std::string Q3_to_string; /// Filename for limiting Q3 (usually it's V, upper limit)
	std::string Q4_from_string; /// Filename for limiting Q4 (usually it's K, lower limit)
	std::string Q4_to_string; /// Filename for limiting Q4 (usually it's K, upper limit)

	double last_update_time; /// Indicates when the Matrix was updated last time
	bool is_limited = false;

	/// Last line position in the matrix scaling and limiting file
	//long scale_pos=0, update_pos=0, Q1_from_pos=0, Q1_to_pos=0, Q2_from_pos=0, Q2_to_pos=0, Q3_from_pos=0, Q3_to_pos=0, Q4_from_pos=0, Q4_to_pos=0; 

	// storing limiting and scale coefficients so we know if values have changed
	double latest_scale_coefficient = 1;
	double latest_Q1_from = -1e99;
	double latest_Q1_to = 1e99;
	double latest_Q2_from = -1e99;
	double latest_Q2_to = 1e99;
	double latest_Q3_from = -1e99;
	double latest_Q3_to = 1e99;
	double latest_Q4_from = -1e99;
	double latest_Q4_to = 1e99;

    std::queue<QueueElement<std::string>> lst_update_queue;
    std::queue<QueueElement<double>> scale_update_queue;

    std::queue<QueueElement<double>> Q1_from_update_queue;
    std::queue<QueueElement<double>> Q1_to_update_queue;
    std::queue<QueueElement<double>> Q2_from_update_queue;
    std::queue<QueueElement<double>> Q2_to_update_queue;
    std::queue<QueueElement<double>> Q3_from_update_queue;
    std::queue<QueueElement<double>> Q3_to_update_queue;
    std::queue<QueueElement<double>> Q4_from_update_queue;
    std::queue<QueueElement<double>> Q4_to_update_queue;

public:
	// original matrix
	MatrixND original_arr;

	// to simplify things
	MatrixND& operator= (const MatrixND &M);
	MatrixND& operator= (const double Val);

	void saveCurrent(); // Save current matrix as original to be scaled/limited later

	//void AllocateMemory( int size_Q1 );
	//void AllocateMemory( int size_Q1, int size_Q2 );
	//void AllocateMemory( int size_Q1, int size_Q2, int size_Q3 );
	//void AllocateMemory( int size_Q1, int size_Q2, int size_Q3, int size_Q4 );

	//void readFromFile(std::string filename, const MatrixND grid_w); ///< Load matrix to a file
	//void readFromFile(std::string filename, const MatrixND grid_w, const MatrixND grid_x); ///< Load matrix to a file
	//void readFromFile(std::string filename, const MatrixND grid_w, const MatrixND grid_x, const MatrixND grid_y); ///< Load matrix to a file
	//void readFromFile(std::string filename, const MatrixND grid_w, const MatrixND grid_x, const MatrixND grid_y, const MatrixND grid_z); ///< Load matrix to a file

	bool readFromIniFile(const std::string& ini_filename, const MatrixND& Q1, const MatrixND& Q2, const MatrixND& Q3, const MatrixND& Q4 = MatrixND());
	bool readFromString(const std::string& file_line_string, const MatrixND &Q1, const MatrixND &Q2, const MatrixND &Q3, const MatrixND &Q4);
	bool update(double time, const MatrixND& Q1, const MatrixND& Q2, const MatrixND& Q3, const MatrixND& Q4 = MatrixND());
};

/**
 * Updatable list of matrix - 1,2,3,4 dimensions. Currently set up to be used for diffusion coefficients, Sources, and Losses only.
 * All other matrices use UpdatableMatrix
 *
 * It's a template because it's easier to create, maintain, and debug one class, instead of 4 classes with the same functionality
 *
 * It's just like a normal matrix (and can be used just like one),
 * but it can also be updated according to rules from an ini-file
 *
 * The ini-files specify several matrices that are added to each other,
 * e.g. diffusion coefficients for different waves to get one diffusion coefficient
 * The readFromIniFile function is used to load data into an UpdatableListMatrix by storing line by line the contents of the ini file into a vector of UpdatableMatrix
 * At every time step the update function is called from the main file. The vector is then traversed and any valid updates are completed.
 * The standard method is to save diffusion coefficients, Sources, or Losses in a .tab file which is comprised of either .plt or .lst files or both.
 * 
 * 
 * @brief A matrix that can be created in 1, 2, 3, or 4 dimensions with the ability to be updated
 */
template <typename MatrixND>
class UpdatableListMatrix : public MatrixND {
private:
	/// These are the matrices that need to be combined to get the final, updated matrix
	std::vector <UpdatableMatrix<MatrixND>> matricesList;

public:

	enum MERGE_TYPE{SUM, MEAN};

	UpdatableListMatrix() {};
	UpdatableListMatrix(MERGE_TYPE merge_type);

	// to simplify things
	MatrixND& operator= (const MatrixND &M);
	MatrixND& operator= (const double Val);

	bool readFromIniFile(std::string ini_filename, const MatrixND &Q1, const MatrixND &Q2, const MatrixND &Q3, const MatrixND &Q4);
	bool update(double time, const MatrixND& Q1, const MatrixND& Q2, const MatrixND& Q3, const MatrixND& Q4 = MatrixND());
	void clearMatricesList();

private:
	MERGE_TYPE merge_type = MERGE_TYPE::SUM;

};

// couple of helpful functions
std::string GetCurrentTimeValue(const std::string &filename, const double &current_time, double &update_time);
std::string GetCurrentTimeValue(const std::string &filename, long &pos, const double &current_time, double &update_time);
bool is_number(const std::string& s);
/// FUNCTION NOT IMPLEMENTED
double stringToValue(std::string string_value, double current_time);
