/*
 * ReadInitialData.h
 *
 *  Created on: Oct 5, 2012
 *      Author: dimath
 */

#ifndef READINITIALDATA_H_
#define READINITIALDATA_H_

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>

//#include <strings.h>
#include <string>
using namespace std;

#include "Logger.h"

#include "UpdatableMatrix.h"
#include "Matrix.h"

#include "Parameters.h"




//bool ReadInitialData(string InputFolder,
//		int &P_size, int &R_size, int &I_size, int &K_size, int &L_size,
//		long int &it_total, double &dt, double &output_time, double &total_time, int &output_step);
bool ReadInitialData(string &InputFolder, string &OutputFolder, int argc, char* argv[],
		double &total_time, double &dt, double &output_time, double &time_first, long int &it_first, int &max_threads,
		string &inversion_method,
		Matrix4D<double> &PSD,
		Matrix4D<double> &P, Matrix4D<double> &R, Matrix4D<double> &V, Matrix4D<double> &K, UpdatableMatrix < Matrix4D<double> > &L,
		int &P_size, int &R_size, int &V_size, int &K_size, int &L_size,
		Matrix3D<double> &P_LBC, Matrix3D<double> &P_UBC, UpdatableMatrix< Matrix3D<double> > &R_LBC, UpdatableMatrix< Matrix3D<double> > &R_UBC, UpdatableMatrix< Matrix3D<double> > &V_LBC, UpdatableMatrix< Matrix3D<double> >&V_UBC, UpdatableMatrix< Matrix3D<double> > &K_LBC, UpdatableMatrix< Matrix3D<double> > &K_UBC, UpdatableMatrix< Matrix3D<double> > &L_LBC, UpdatableMatrix< Matrix3D<double> > &L_UBC,
		string &P_LBC_type, string &P_UBC_type, string &R_LBC_type, string &R_UBC_type, string &V_LBC_type, string &V_UBC_type, string &K_LBC_type, string &K_UBC_type, string &L_LBC_type, string &L_UBC_type,
		UpdatableListMatrix < Matrix4D<double> > &DLL, UpdatableListMatrix < Matrix4D<double> > &DVV, UpdatableListMatrix < Matrix4D<double> > &DKK, UpdatableListMatrix < Matrix4D<double> > &DVK,
		UpdatableMatrix < Matrix4D<double> > &VP, UpdatableMatrix < Matrix4D<double> > &VL,
		UpdatableMatrix < Matrix4D<double> > &G_local, UpdatableMatrix < Matrix4D<double> >  &G_radial,
		UpdatableListMatrix < Matrix4D<double> > &Sources, UpdatableListMatrix < Matrix4D<double> > &Losses);


#endif /* READINITIALDATA_H_ */
