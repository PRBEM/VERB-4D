// SPDX-FileCopyrightText: 2015 UCLA
// SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * \file Parameters.cpp
 * \brief Holds the Parameters along with their corresponding file and arguments 
 *
 *  Created on: Nov 17, 2012
 *      Author: dimath
 */
#include "Parameters.h"
#include "BoundaryConditionType.hpp"

/** Constructor for Parameters: 
* opens file and saves file to Parameters field 'parametersFile',
* stores the char* arguments into a vector 'argv'
*/
Parameters::Parameters(std::string filename, int argc, char* argv[]) {

	// open the parameter file
	this->parametersFile.open(filename.c_str());
			//if (this->parametersFile == NULL) {
	// if file not found - matrix is zero
	if (!this->parametersFile.is_open()) {		
		Logger::error << filename << " not found." << std::endl;
		exit(EXIT_FAILURE);
	}

	// for every argument - store them into argv data field
	for (int i = 1; i<argc; i++) {
		this->argv.push_back(argv[i]);
	}
}

/// Destructor - tells this parameter to close its file
Parameters::~Parameters() {
	this->close();
}

/// close parameter file
void Parameters::close() {
	parametersFile.close();
}

/**
* Finds the parameters with parameterName in parametersFile and saves it to argv
* \param parameterName - name of parameter
*/
Parameters& Parameters::findParameter(std::string parameterName, std::string defaultValue) {

	this->flush();
	this->clear();
	this->seekg( 0, std::ios::beg );

	parametersFile.clear( );
	parametersFile.seekg( 0, std::ios::beg );

	std::string line;

	// Read from command line arguments
	unsigned int i;
	// for every argument
	for (i = 0; i < this->argv.size(); i++) {
		
		// store the arguments one at a time as a string
		line = argv[i];
		// if argument is empty or commented out, ignore this argument
		if (line.size() == 0 || line[0] == '#') {
			continue;
		}
		// set the parameter value if it can find its name in the current argument
		// return the parameter once it is set
		if (line.find(parameterName) != std::string::npos) {
			this->str(line.substr(line.find("=") + 1));
			return *this;
		}
	}

	// Read from parameters file
	// store the lines one at a time as a string
	while (std::getline(parametersFile, line)) {

		// if line is empty or commented out, ignore this line
		if (line.size() == 0 || line[0] == '#') {
			continue;
		}

		// set the parameter value if it can find its name in the current line
		// return the parameter once it is set
		if (line.find(parameterName) != std::string::npos) {
			this->str(line.substr(line.find("=") + 1));
			return *this;
		}
	}

	// If no parameter value can be found log error and return
	if (defaultValue == "") {
		Logger::error << "Parameter " << parameterName << " not found." << std::endl;
		exit(EXIT_FAILURE);
	}
	return *this;

}


/**
* Gets the value for a  parameters with parameterName and saves it into variable
* \param parameterName - name of parameter
* \param variable - variable for storing parameter
* \param mustBeFound - if true and parameter not found, error is logged and function exited, default is false
*/
template <typename T>
bool Parameters::getParameter(std::string parameterName, T &variable, bool mustBeFound) {

	this->flush();
	this->clear();
	this->seekg( 0, std::ios::beg );

	// start reading parameter file from beginning
	parametersFile.clear( );
	parametersFile.seekg( 0, std::ios::beg );

	// variables to store each line and the value of the parameter
	std::string line;
	std::string parameterValue;

	// Read from command line arguments
	unsigned int i;
	for (i = 0; i < this->argv.size(); i++) {
		line = argv[i];

		// if argument is empty or commented out, save the default value
		if (line.size() == 0 || line[0] == '#') {
			//Logger::message << parameterName << " = " << variable << " (default value)" << std::endl;
			continue;
		}

		// set the parameter value if it can find its name in the current argument
		// log the value and store the value into variable
		if (line.find(parameterName) != std::string::npos) {
			parameterValue = line.substr(line.find("=") + 1);
			Logger::message << parameterName << " = " << parameterValue << std::endl;
			stringToValue(parameterValue, variable);
			return true;
		}
	}

	// Read from parameters file - line by line
	while (std::getline(parametersFile, line)) {

		// if line is empty or commented out, go to next line
		if (line.size() == 0 || line[0] == '#') {
			continue;
		}

		// set the parameter value if it can find its name in the current line
		// log the value and store the value into variable
		if (line.find(parameterName) != std::string::npos) {
			parameterValue = line.substr(line.find("=") + 2); // +2 strips of the whitespace
			Logger::message << parameterName << " = " << parameterValue << std::endl;
			stringToValue(parameterValue, variable);
			return true;
		}
	}

	// all variables that arent found, that can't recieve a default value will log an error
	if (mustBeFound) {
		Logger::error << "Parameter " << parameterName << " not found." << std::endl;
		exit(EXIT_FAILURE);
	}

	// return *this;
	return false;

}

//////////////////////////////////////////
// Implementations
//////////////////////////////////////////

template bool Parameters::getParameter(std::string, double&, bool);
template bool Parameters::getParameter(std::string, int&, bool);
template bool Parameters::getParameter(std::string, long&, bool);
template bool Parameters::getParameter(std::string, bool&, bool);
template bool Parameters::getParameter(std::string, std::string&, bool);
template bool Parameters::getParameter(std::string, BoundaryConditionType&, bool);
template bool Parameters::getParameter(std::string, IOMethod&, bool);
template bool Parameters::getParameter(std::string, InversionMethod&, bool);
template bool Parameters::getParameter(std::string, DensitySaturation&, bool);
template bool Parameters::getParameter(std::string, DataAssimilationDataSource&, bool);
template bool Parameters::getParameter(std::string, Logger::DebugLevel&, bool);


// helper functions to insensitive case compare two strings
bool ichar_equals(char a, char b)
{
    return std::tolower(static_cast<unsigned char>(a)) ==
           std::tolower(static_cast<unsigned char>(b));
}
bool iequals(const std::string& a, const std::string& b)
{
    return std::equal(a.begin(), a.end(), b.begin(), b.end(), ichar_equals);
}

// helper functions to convert paramter strings into values
template <typename T>
void stringToValue(const std::string& parameter_value_string, T& variable) {
	std::stringstream tmp(parameter_value_string);
	tmp >> variable;
}

template <>
void stringToValue(const std::string& parameter_value_string, bool& variable) {
	if (iequals(parameter_value_string, "true")) {
		variable = true;
	} else if (iequals(parameter_value_string, "false")) {
		variable = false;
	} else {
		printf("Encountered invalid value for boolean parameter: %s\n", parameter_value_string.c_str());
		exit(EXIT_FAILURE);
	}
}

template <>
void stringToValue(const std::string& parameter_value_string, DensitySaturation& variable) {
	if (iequals(parameter_value_string, "off")) {
		variable = DensitySaturation::Off;
	} else if (iequals(parameter_value_string, "with_timescale")) {
		variable = DensitySaturation::WithTimescale;
	} else if (iequals(parameter_value_string, "without_timescale")) {
		variable = DensitySaturation::WithoutTimescale;
	} else if (iequals(parameter_value_string, "limit_source")) {
		variable = DensitySaturation::LimitSource;
	} else {
		printf("Encountered invalid value for density saturation parameter: %s\n", parameter_value_string.c_str());
		exit(EXIT_FAILURE);
	}
}

template <>
void stringToValue(const std::string& parameter_value_string, IOMethod& variable) {
	if (iequals(parameter_value_string, "matlab")) {
		variable = IOMethod::Matlab;
	} else if (iequals(parameter_value_string, "binary")) {
		variable = IOMethod::Binary;
	} else if (iequals(parameter_value_string, "ascii")) {
		variable = IOMethod::ASCII;
	} else {
		printf("Encountered invalid value for io method parameter: %s\n", parameter_value_string.c_str());
		exit(EXIT_FAILURE);
	}
}

template <>
void stringToValue(const std::string& parameter_value_string, InversionMethod& variable) {
	if (iequals(parameter_value_string, "ADI")) {
		variable = InversionMethod::ADI;
	} else if (iequals(parameter_value_string, "ADI1")) {
		variable = InversionMethod::ADI1;
	} else if (iequals(parameter_value_string, "ADI2")) {
		variable = InversionMethod::ADI2;
	} else if (iequals(parameter_value_string, "Lapack")) {
		variable = InversionMethod::Lapack;
	} else if (iequals(parameter_value_string, "MKL")) {
		variable = InversionMethod::MKL;
	} else {
		printf("Encountered invalid value for inversion method parameter: %s\n", parameter_value_string.c_str());
		exit(EXIT_FAILURE);
	}
}

template <>
void stringToValue(const std::string& parameter_value_string, DataAssimilationDataSource& variable) {
	if (iequals(parameter_value_string, "data_server")) {
		variable = DataAssimilationDataSource::DataServer;
	}
	else if (iequals(parameter_value_string, "local_files")) {
		variable = DataAssimilationDataSource::LocalFiles;
	} else {
		printf("Encountered invalid value for data assimilation data source parameter: %s\n", parameter_value_string.c_str());
		exit(EXIT_FAILURE);
	}
}

template <>
void stringToValue(const std::string& parameter_value_string, Logger::DebugLevel& variable) {
	int enum_int;
	std::stringstream tmp(parameter_value_string);
	tmp >> enum_int;

	if (enum_int < 0 or enum_int > 4) {
		printf("Encountered invalid enum value!");
		exit(EXIT_FAILURE);
	}

	variable = static_cast<Logger::DebugLevel>(enum_int);
}