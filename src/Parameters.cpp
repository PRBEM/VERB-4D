/*
 * Parameters.cpp
 *
 *  Created on: Nov 17, 2012
 *      Author: dimath
 */

#include "Parameters.h"


using namespace std;

Parameters::Parameters(string filename, int argc, char* argv[]) {
	this->parametersFile.open(filename.c_str());
	//if (this->parametersFile == NULL) {
	if (!this->parametersFile.is_open()) {
		// if file not found - matrix is zero
		cout << filename << " not found." << endl;		
		
		exit(EXIT_FAILURE);
	}

	for (int i = 1; i<argc; i++) {
		this->argv.push_back(argv[i]);
	}
}

Parameters::~Parameters() {
	this->close();
}

void Parameters::close() {
	parametersFile.close();
}

Parameters& Parameters::findParameter(string parameterName, string defaultValue) {

	this->flush();
	this->clear();
	this->seekg( 0, std::ios::beg );

	parametersFile.clear( );
	parametersFile.seekg( 0, std::ios::beg );

	string line;

	// Read from command line arguments
	unsigned int i;
	for (i = 0; i < this->argv.size(); i++) {
		line = argv[i];

		if (line.size() == 0 || line[0] == '#') {
			continue;
		}

		if (line.find(parameterName) != -1) {
			this->str(line.substr(line.find("=") + 1));
			return *this;
		}
	}

	// Read from parameters file
	while (std::getline(parametersFile, line)) {

		if (line.size() == 0 || line[0] == '#') {
			continue;
		}

		if (line.find(parameterName) != -1) {
			this->str(line.substr(line.find("=") + 1));
			return *this;
		}
	}

	if (defaultValue == "") {
		cout << "Parameter " << parameterName << " not found." << endl;
		exit(EXIT_FAILURE);
	}
	return *this;

}


template <typename T>
void Parameters::getParameter(string parameterName, T &variable, bool mustBeFound) {

	this->flush();
	this->clear();
	this->seekg( 0, std::ios::beg );

	parametersFile.clear( );
	parametersFile.seekg( 0, std::ios::beg );

	string line;
	string parameterValue;

	// Read from command line arguments
	unsigned int i;
	for (i = 0; i < this->argv.size(); i++) {
		line = argv[i];

		if (line.size() == 0 || line[0] == '#') {
			cout << parameterName << " = " << variable << " (default value)" << endl;
			continue;
		}

		if (line.find(parameterName) != -1) {
			parameterValue = line.substr(line.find("=") + 1);
			cout << parameterName << " = " << parameterValue << endl;
			// stringstream(parameterValue) >> variable;
			stringstream tmp(parameterValue);
			tmp >> variable;
			return;
		}
	}

	// Read from parameters file
	while (std::getline(parametersFile, line)) {

		if (line.size() == 0 || line[0] == '#') {
			continue;
		}

		if (line.find(parameterName) != -1) {
			parameterValue = line.substr(line.find("=") + 1);
			cout << parameterName << " = " << parameterValue << endl;
			// stringstream(parameterValue) >> variable;
			stringstream tmp(parameterValue);
			tmp >> variable;
			return;
		}
	}

	if (mustBeFound) {
		cout << "Parameter " << parameterName << " not found." << endl;
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
