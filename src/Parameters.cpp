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
void Parameters::getParameter(std::string parameterName, T &variable, bool mustBeFound) {

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
			return;
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
			return;
		}
	}

	// all variables that arent found, that can't recieve a default value will log an error
	if (mustBeFound) {
		Logger::error << "Parameter " << parameterName << " not found." << std::endl;
		exit(EXIT_FAILURE);
	}

	// return *this;
	return;

}

//////////////////////////////////////////
// Implementations
//////////////////////////////////////////

template void Parameters::getParameter(std::string, double&, bool);
template void Parameters::getParameter(std::string, int&, bool);
template void Parameters::getParameter(std::string, long&, bool);
template void Parameters::getParameter(std::string, bool&, bool);
template void Parameters::getParameter(std::string, std::string&, bool);
template void Parameters::getParameter(std::string, BoundaryConditionType&, bool);
template void Parameters::getParameter(std::string, IOMethod&, bool);
template void Parameters::getParameter(std::string, InversionMethod&, bool);


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
		printf("Encountered invalid value for boolean parameter:%s\n", parameter_value_string.c_str());
		exit(EXIT_FAILURE);
	}
}

template <>
void stringToValue(const std::string& parameter_value_string, DensitySaturation& density_saturation) {
	if (iequals(parameter_value_string, "off")) {
		density_saturation = DensitySaturation::Off;
		return;
	}
	if (iequals(parameter_value_string, "with_timescale")) {
		density_saturation = DensitySaturation::WithTimescale;
		return;
	}
	if (iequals(parameter_value_string, "without_timescale")) {
		density_saturation = DensitySaturation::WithoutTimescale;
		return;
	}
	if (iequals(parameter_value_string, "limit_source")) {
		density_saturation = DensitySaturation::LimitSource;
		return;
	}
}

template <>
void stringToValue(const std::string& parameter_value_string, IOMethod& io_method) {
	if (iequals(parameter_value_string, "matlab")) {
		io_method = IOMethod::Matlab;
		return;
	}
	if (iequals(parameter_value_string, "binary")) {
		io_method = IOMethod::Binary;
		return;
	}
	if (iequals(parameter_value_string, "ascii")) {
		io_method = IOMethod::ASCII;
		return;
	}
}

template <>
void stringToValue(const std::string& parameter_value_string, InversionMethod& inversion_method) {
	if (iequals(parameter_value_string, "ADI")) {
		inversion_method = InversionMethod::ADI;
		return;
	}
	if (iequals(parameter_value_string, "ADI1")) {
		inversion_method = InversionMethod::ADI1;
		return;
	}
	if (iequals(parameter_value_string, "ADI2")) {
		inversion_method = InversionMethod::ADI2;
		return;
	}
		if (iequals(parameter_value_string, "Lapack")) {
		inversion_method = InversionMethod::Lapack;
		return;
	}
	if (iequals(parameter_value_string, "MKL")) {
		inversion_method = InversionMethod::MKL;
		return;
	}
}