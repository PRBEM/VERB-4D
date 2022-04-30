/**
 * \file Parameters.h
 * \brief Holds the parameters along with their corresponding file and arguments 
 *
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
	void getParameter(std::string parameterName, T &variable, bool mustBeFound = false);


};


#endif /* PARAMETERS_H_ */
