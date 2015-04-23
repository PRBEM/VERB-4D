/*
 * Parameters.h
 *
 *  Created on: Nov 17, 2012
 *      Author: dimath
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include <stdio.h>
#include <stdlib.h>

using namespace std;

class Parameters : public stringstream {
private:
	ifstream parametersFile;
	vector<string> argv;

public:
	Parameters(string filename, int argc = 1, char* argv[] = NULL);
	~Parameters();

	void close();

	Parameters& findParameter(string parameterName, string defaultValue = "");

	template <typename T>
	void getParameter(string parameterName, T &variable, bool mustBeFound = false);


};


#endif /* PARAMETERS_H_ */
