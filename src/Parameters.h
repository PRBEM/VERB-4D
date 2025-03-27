/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * \file Parameters.h
 * \brief Holds the parameters along with their corresponding file and arguments 
 *
 *
 *  Created on: Nov 17, 2012
 *      Author: dimath
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
@brief Holds the parameters along with their corresponding file and arguments 
*/
class Parameters : public std::stringstream {
private:
	/// store the file contents being read in
	std::ifstream parametersFile;
	/// store the arguments being sent in
	std::vector<std::string> argv;

public:
	Parameters(std::string filename, int argc = 1, char* argv[] = NULL);
	~Parameters();

	void close();

	/// Find the parameters given and return them as a parameter struct
	Parameters& findParameter(std::string parameterName, std::string defaultValue = "");

	/// get the parameter value and store it into variable
	template <typename T>
	bool getParameter(std::string parameterName, T &variable, bool mustBeFound = false);


};


template <typename T>
void stringToValue(const std::string& parameter_value_string, T& variable);

enum class IOMethod {ASCII, Binary, Matlab};
enum class InversionMethod {ADI, ADI1, ADI2, MKL, Lapack};
enum class DensitySaturation {Off, WithTimescale, WithoutTimescale, LimitSource};
enum class DataAssimilationDataSource {DataServer, LocalFiles};