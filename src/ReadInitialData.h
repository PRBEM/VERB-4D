/**
 * \file ReadInitialData.h
 * \brief Reads in all the data files created in Matlab in order to initialize the parameters
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

#include "Logger.h"

#include "UpdatableMatrix.h"
#include "Matrix.h"
#include "BoundaryConditionType.hpp"
#include "Parameters.h"




//bool ReadInitialData(std::string InputFolder,
//		int &P_size, int &R_size, int &I_size, int &K_size, int &L_size,
//		long int &it_total, double &dt, double &output_time, double &total_time, int &output_step);
bool ReadInitialData(std::string &InputFolder, std::string &OutputFolder, int argc, char* argv[],
	double &time_total, double &time_step, double &time_output, double &time_first, long int &it_first, int &max_threads,
	InversionMethod &inversion_method, IOMethod &io_method, IOMethod &PSD0_io_method, DensitySaturation &density_saturation,
	bool &include_boundary, bool &Vl_BC_from_convection, bool &Vu_BC_from_convection, bool &run_remapping,
	bool &run_convection, bool &run_radial_diffusion, bool &run_local_diffusion, bool &positive_PSD, bool &PSD_time_to_lst,
	Matrix4D<double> &PSD,
	Matrix4D<double> &P, Matrix4D<double> &R, Matrix4D<double> &V, Matrix4D<double> &K, UpdatableMatrix < Matrix4D<double> > &L,
	int &P_size, int &R_size, int &V_size, int &K_size, int &L_size,
	Matrix3D<double> &PSD_l_P, Matrix3D<double> &PSD_u_P, 
	UpdatableMatrix<Matrix3D<double>> &PSD_l_R, UpdatableMatrix<Matrix3D<double>> &PSD_u_R, 
	UpdatableMatrix<Matrix3D<double>> &PSD_l_V, UpdatableMatrix<Matrix3D<double>> &PSD_u_V, 
	UpdatableMatrix<Matrix3D<double>> &PSD_l_K, UpdatableMatrix<Matrix3D<double>> &PSD_u_K, 
	UpdatableMatrix<Matrix3D<double>> &PSD_l_L, UpdatableMatrix<Matrix3D<double>> &PSD_u_L,
	BoundaryConditionType &Pl_BC_type, BoundaryConditionType &Pu_BC_type, 
	BoundaryConditionType &Rl_BC_type, BoundaryConditionType &Ru_BC_type, 
	BoundaryConditionType &Vl_BC_type, BoundaryConditionType &Vu_BC_type, 
	BoundaryConditionType &Kl_BC_type, BoundaryConditionType &Ku_BC_type, 
	BoundaryConditionType &Ll_BC_type, BoundaryConditionType &Lu_BC_type,
	UpdatableListMatrix<Matrix4D<double>> &DLL, UpdatableListMatrix<Matrix4D<double>> &DVV, 
	UpdatableListMatrix<Matrix4D<double>> &DKK, UpdatableListMatrix<Matrix4D<double>> &DVK,
	UpdatableMatrix<Matrix4D<double>> &VP, UpdatableMatrix<Matrix4D<double>> &VL,
	UpdatableMatrix<Matrix4D<double>> &G_local, UpdatableMatrix<Matrix4D<double>> &G_radial,
	UpdatableListMatrix<Matrix4D<double>> &Sources, UpdatableListMatrix<Matrix4D<double>> &Losses, 
	UpdatableListMatrix<Matrix4D<double>> &Losses_conv, Matrix4D<double> &SaturationDensity, Matrix4D<double> &SaturationTimescale
);


// Implemented function but not ever used
/* void ReadBoundaryCondition(
		ifstream &input,
		Matrix3D<double> psd_slice,
		// UpdatableMatrix2D &BC,
		Matrix3D<double> &BC,
		std::string &BC_type,
		Matrix3D<double> Q_slice1, Matrix3D<double> Q_slice2, Matrix3D<double> Q_slice3)
*/

#endif /* READINITIALDATA_H_ */
