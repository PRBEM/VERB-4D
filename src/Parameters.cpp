/**
 * \file Parameters.cpp
 * \brief Holds the Parameters along with their corresponding file and arguments 
 *
 *  Created on: Nov 17, 2012
 *      Author: dimath
 */
#include "Parameters.h"

/** Constructor for Parameters: 
* opens file and saves file to Parameters field 'parametersFile',
* stores the char* arguments into a vector 'argv'
*/
Parameters::Parameters(string filename, int argc, char* argv[]) {

	// open the parameter file
	this->parametersFile.open(filename.c_str());
			//if (this->parametersFile == NULL) {
	// if file not found - matrix is zero
	if (!this->parametersFile.is_open()) {		
		Logger::error << filename << " not found." << endl;
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
Parameters& Parameters::findParameter(string parameterName, string defaultValue) {

	this->flush();
	this->clear();
	this->seekg( 0, std::ios::beg );

	parametersFile.clear( );
	parametersFile.seekg( 0, std::ios::beg );

	string line;

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
		Logger::error << "Parameter " << parameterName << " not found." << endl;
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
void Parameters::getParameter(string parameterName, T &variable, bool mustBeFound) {

	this->flush();
	this->clear();
	this->seekg( 0, std::ios::beg );

	// start reading parameter file from beginning
	parametersFile.clear( );
	parametersFile.seekg( 0, std::ios::beg );

	// variables to store each line and the value of the parameter
	string line;
	string parameterValue;

	// Read from command line arguments
	unsigned int i;
	for (i = 0; i < this->argv.size(); i++) {
		line = argv[i];

		// if argument is empty or commented out, save the default value
		if (line.size() == 0 || line[0] == '#') {
			Logger::message << parameterName << " = " << variable << " (default value)" << endl;
			continue;
		}

		// set the parameter value if it can find its name in the current argument
		// log the value and store the value into variable
		if (line.find(parameterName) != std::string::npos) {
			parameterValue = line.substr(line.find("=") + 1);
			Logger::message << parameterName << " = " << parameterValue << endl;
					// stringstream(parameterValue) >> variable;
			stringstream tmp(parameterValue);
			tmp >> variable;
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
			parameterValue = line.substr(line.find("=") + 1);
			Logger::message << parameterName << " = " << parameterValue << endl;
			// stringstream(parameterValue) >> variable;
			stringstream tmp(parameterValue);
			tmp >> variable;
			return;
		}
	}

	// all variables that arent found, that can't recieve a default value will log an error
	if (mustBeFound) {
		Logger::error << "Parameter " << parameterName << " not found." << endl;
		exit(EXIT_FAILURE);
	}

	// return *this;
	return;

}

//////////////////////////////////////////
// Implementations
//////////////////////////////////////////

template void Parameters::getParameter(string, double&, bool);
template void Parameters::getParameter(string, int&, bool);
template void Parameters::getParameter(string, long&, bool);
//template void Parameters::getParameter(string, bool&, bool);
template void Parameters::getParameter(string, string&, bool);
