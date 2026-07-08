/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * \file Parameters.h
 * \brief Parameter file parsing and command line argument handling
 *
 * This file implements a parameter management system that reads configuration
 * parameters from files and command line arguments. It supports type-safe
 * parameter extraction with default values and validation.
 *
 * Created on: Nov 17, 2012
 * Author: dimath
 */

#pragma once

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cctype>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>

#include "Logger.h"

/**
 * @brief Parameter file and command line argument parser
 * 
 * @details This class provides functionality to read configuration parameters
 * from parameter files and command line arguments. It inherits from std::stringstream
 * to provide stream-based parameter extraction capabilities.
 * 
 * **Parameter Resolution Order:**
 * 1. Command line arguments (highest priority)
 * 2. Parameter file values
 * 3. Default values (if provided)
 * 
 * **Parameter File Format:**
 * - Lines starting with '#' are treated as comments
 * - Parameters use format: parameter_name = value
 * - Empty lines are ignored
 * 
 * **Command Line Format:**
 * - Arguments use format: parameter_name=value
 * - No spaces around the '=' sign
 */
class Parameters : public std::stringstream {
private:
	/** @brief Input file stream for reading parameter file */
	std::ifstream parametersFile;
	
	/** @brief Vector storing command line arguments */
	std::vector<std::string> argv;

public:
	/**
	 * @brief Constructs Parameters object from file and command line arguments
	 * 
	 * @param filename [in] Path to parameter file to read
	 * @param argc [in] Number of command line arguments (default: 1)
	 * @param argv [in] Array of command line argument strings (default: NULL)
	 * 
	 * @throws Exits program if parameter file cannot be opened
	 */
	Parameters(std::string filename, int argc = 1, char* argv[] = NULL);
	
	/**
	 * @brief Destructor that closes the parameter file
	 */
	~Parameters();

	/**
	 * @brief Closes the parameter file stream
	 */
	void close();

	/**
	 * @brief Finds a parameter and sets up stream for reading its value
	 * 
	 * @param parameterName [in] Name of parameter to find
	 * @param defaultValue [in] Default value if parameter not found (default: "")
	 * 
	 * @return Reference to this Parameters object for chaining
	 * 
	 * @throws Exits program if parameter not found and no default provided
	 */
	Parameters& findParameter(std::string parameterName, std::string defaultValue = "");

	/**
	 * @brief Gets parameter value and stores it in a typed variable
	 * 
	 * @tparam T Type of the variable to store parameter value
	 * @param parameterName [in] Name of parameter to retrieve
	 * @param variable [out] Variable to store the parameter value
	 * @param mustBeFound [in] If true, exits program when parameter not found (default: false)
	 * 
	 * @return true if parameter was found and parsed successfully, false otherwise
	 * 
	 * @note Supports automatic type conversion for common types (int, double, bool, string, enums)
	 */
	template <typename T>
	bool getParameter(std::string parameterName, T &variable, bool mustBeFound = false);

};


/**
 * @brief Converts string parameter value to typed variable
 * 
 * @tparam T Type to convert the string to
 * @param parameter_value_string [in] String value from parameter file or command line
 * @param variable [out] Variable to store the converted value
 * 
 * @note Template specializations exist for bool, enums, and other special types
 */
template <typename T>
void stringToValue(const std::string& parameter_value_string, T& variable);

/**
 * @brief Input/Output method enumeration
 */
enum class IOMethod {
	ASCII,   /**< ASCII text format */
	Binary,  /**< Binary format */
	Matlab   /**< MATLAB format */
};

/**
 * @brief Matrix inversion method enumeration
 */
enum class InversionMethod {
	ADI,     /**< Alternating Direction Implicit method */
	ADI1,    /**< ADI variant 1 */
	ADI2,    /**< ADI variant 2 */
	MKL,     /**< Intel Math Kernel Library */
	Lapack,  /**< LAPACK library */
	PPFV     /**< positive-preserving finite volume method */
};

/**
 * @brief Density saturation handling enumeration
 */
enum class DensitySaturation {
	Off,              /**< No density saturation */
	WithTimescale,    /**< Saturation with timescale */
	WithoutTimescale, /**< Saturation without timescale */
	LimitSource       /**< Limit source terms */
};

/**
 * @brief Data assimilation data source enumeration
 */
enum class DataAssimilationDataSource {
	DataServer,  /**< Remote data server */
	LocalFiles   /**< Local file system */
};