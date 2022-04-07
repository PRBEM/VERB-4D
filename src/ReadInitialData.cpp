/**
 * \file ReadInitialData.cpp
 * \brief Reads in all the data files created in Matlab in order to initialize the parameters
 *
 *  Created on: Oct 5, 2012
 *      Author: dimath
 */

#include "ReadInitialData.h"

using namespace std;

#if defined(_WIN32) || defined(_WIN64)
	#define strncasecmp _strnicmp
	#define strcasecmp _stricmp
#endif

/**
* Allocates the memory for all the matrices using Matrix4D::AllocateMemory() for everything except boundary conditions which use Matrix3D::AllocateMemory()
*
* \param P,R,V,K,L - matrices for Phi, Radial Distance, Energy, Pitch Angle, Distance due to magnetic field, respectively
* \param P,R,V,K,L _size - sizes for the corresponding matrices
* \param UBC/LBC - Boundary condition for each matrix P,R,V,K,L upper and lower
* \param DLL,DVV,DKK,DVK - diffusion matrices
* \param VP - Phi and Energy Matrix
* \param VL - Phi and distance due to magnetic field Matrix
* \param G_local, G_radial - Jacobians for normalizing data
* \param Sources, Losses - matrices for calculating Sources and Losses(loss cone)
*/
void AllocateMemory(Matrix4D<double> &PSD,
		Matrix4D<double> &P, Matrix4D<double> &R, Matrix4D<double> &V,  Matrix4D<double> &K, Matrix4D<double> &L,
		int &P_size, int &R_size, int &V_size, int &K_size,
		Matrix3D<double> &P_LBC, Matrix3D<double> &P_UBC, Matrix3D<double> &R_LBC, Matrix3D<double> &R_UBC, Matrix3D<double> &V_LBC, Matrix3D<double> &V_UBC, Matrix3D<double> &K_LBC, Matrix3D<double> &K_UBC, Matrix3D<double> &L_LBC, Matrix3D<double> &L_UBC,
		Matrix4D<double> &DLL, Matrix4D<double> &DVV, Matrix4D<double> &DVK, Matrix4D<double> &DKK,
		Matrix4D<double> &VP, Matrix4D<double> &VL,
		Matrix4D<double> &G_local, Matrix4D<double> &G_radial, Matrix4D<double> &Sources, Matrix4D<double> &Losses, Matrix4D<double> &Losses_conv) {


	// Arrays initialization (memory)
	P.AllocateMemory(P_size, R_size, V_size, K_size);
	R.AllocateMemory(P_size, R_size, V_size, K_size);
	V.AllocateMemory(P_size, R_size, V_size, K_size);
	K.AllocateMemory(P_size, R_size, V_size, K_size);
	L.AllocateMemory(P_size, R_size, V_size, K_size);

	PSD.AllocateMemory(P_size, R_size, V_size, K_size);

	DLL.AllocateMemory(P_size, R_size, V_size, K_size);
	DVV.AllocateMemory(P_size, R_size, V_size, K_size);
	DKK.AllocateMemory(P_size, R_size, V_size, K_size);
	DVK.AllocateMemory(P_size, R_size, V_size, K_size);

	VL.AllocateMemory(P_size, R_size, V_size, K_size);
	VP.AllocateMemory(P_size, R_size, V_size, K_size);

	Sources.AllocateMemory(P_size, R_size, V_size, K_size);
	Losses.AllocateMemory(P_size, R_size, V_size, K_size);
    Losses_conv.AllocateMemory(P_size, R_size, V_size, K_size);


	//G_local.AllocateMemory(I_size, K_size); // I,K grid is always the same, so Jacobian is the same
	G_local.AllocateMemory(P_size, R_size, V_size, K_size);
	G_radial.AllocateMemory(P_size, R_size, V_size, K_size); // L-star is different for different P, I, K

	P_LBC.AllocateMemory(R_size, V_size, K_size);
	P_UBC.AllocateMemory(R_size, V_size, K_size);
	R_LBC.AllocateMemory(P_size, V_size, K_size);
	R_UBC.AllocateMemory(P_size, V_size, K_size);
	V_LBC.AllocateMemory(P_size, R_size, K_size);
	V_UBC.AllocateMemory(P_size, R_size, K_size);
	K_LBC.AllocateMemory(P_size, R_size, V_size);
	K_UBC.AllocateMemory(P_size, R_size, V_size);
	L_LBC.AllocateMemory(P_size, V_size, K_size);
	L_UBC.AllocateMemory(P_size, V_size, K_size);

}

/// THIS FUNCTION IS NOT USED ANYWHERE
void ReadBoundaryCondition(
		ifstream &input,
		Matrix3D<double> psd_slice,
		// UpdatableMatrix2D &BC,
		Matrix3D<double> &BC,
		string &BC_type,
		Matrix3D<double> Q_slice1, Matrix3D<double> Q_slice2, Matrix3D<double> Q_slice3) {

	string BC_const;
	input >> BC_type >> BC_const;
	if (BC_type == "BCT_CONSTANT_PSD") {
		// (Matrix2D<double>) BC = psd_slice;
		// BC = psd_slice;
		BC = psd_slice;
		BC_type = "BCT_CONSTANT_VALUE";
	} else if (BC_type == "BCT_FILE") {
		BC.readFromFile(BC_const, Q_slice1, Q_slice2, Q_slice3);
	} else if (BC_type == "BCT_FILE_INI") {
		// XXX
		// BC.readFromIniFile(BC_const, Q_slice1, Q_slice2);
	} else {
		double BC_value;
		stringstream(BC_const) >> BC_value;
		BC = BC_value;
	}

}

//bool ReadInitialData(string InputFolder,
//		int &P_size, int &R_size, int &I_size, int &K_size, int &L_size,
//		long int &it_total, double &dt, double &output_time, double &total_time, int &output_step) {




/** Reads in all the data created in Matlab in order to initialize the parameters
*
* Using Parameters::getParameter() for timesteps, number of threads, input and output folders, boundary condition types
*
* Using AllocateMemory() for all the matrices
*
* Using UpdatableMatrix<MatrixND>::update() for G_local and G_radial
*
* Using UpdatableMatrix<MatrixND>::readFromIniFile() for all boundary conditions and data tables/matrices
*
* \param InputFolder - Folder to get input from
* \param OutputFolder - Folder to write output to
* \param argc - number of arguments given
* \param argv - array of string arguments inputted
* \param time_total,step,output,first - time variables for calculations
* \param max_threads - number of threads that can be used
* \param inversion_method - tool for doing calculations, usually Lapack
* \param PSD - Phase Space Density Matrix
* \param P,R,V,K,L - matrices for Phi, Actual Distance, Energy, Pitch Angle, Distance due to magnetic field, respectively
* \param P,R,V,K,L _size - sizes for the corresponding matrices
* \param PSD_l_ P,R,V,K,L - lower boundaries for PSD with respect to P,R,V,K,L
* \param PSD_u_ P,R,V,K,L - upper boundaries for PSD with respect to P,R,V,K,L
* \param BC_type - Type of boundary condition for each matrix - either BCT_CONSTANT_VALUE or BCT_CONSTANT_DERIVATIVE
* \param DLL,DVV,DKK,DVK - diffusion matrices
* \param G_local, G_radial - Jacobians for normalizing data
* \param Sources, Losses - matrices for calculating Sources and Losses(loss cone)
*/
bool ReadInitialData(string &InputFolder, string &OutputFolder, int argc, char* argv[],
		double &time_total, double &time_step, double &time_output, double &time_first, long int &it_first, int &max_threads,
		string &inversion_method, string &include_boundary, string &Vl_BC_from_convection, string &Vu_BC_from_convection, string &io_method,
		string &run_remapping, string &run_convection, string &run_radial_diffusion, string &run_local_diffusion, string &positive_PSD,
		Matrix4D<double> &PSD,
		Matrix4D<double> &P, Matrix4D<double> &R, Matrix4D<double> &V, Matrix4D<double> &K, UpdatableMatrix < Matrix4D<double> > &L,
		int &P_size, int &R_size, int &V_size, int &K_size, int &L_size,
		Matrix3D<double> &PSD_l_P, Matrix3D<double> &PSD_u_P, UpdatableMatrix< Matrix3D<double> > &PSD_l_R, UpdatableMatrix< Matrix3D<double> > &PSD_u_R, UpdatableMatrix< Matrix3D<double> > &PSD_l_V, UpdatableMatrix< Matrix3D<double> >&PSD_u_V, UpdatableMatrix< Matrix3D<double> > &PSD_l_K, UpdatableMatrix< Matrix3D<double> > &PSD_u_K, UpdatableMatrix< Matrix3D<double> > &PSD_l_L, UpdatableMatrix< Matrix3D<double> > &PSD_u_L,
		string &Pl_BC_type, string &Pu_BC_type, string &Rl_BC_type, string &Ru_BC_type, string &Vl_BC_type, string &Vu_BC_type, string &Kl_BC_type, string &Ku_BC_type, string &Ll_BC_type, string &Lu_BC_type,
		UpdatableListMatrix < Matrix4D<double> > &DLL, UpdatableListMatrix < Matrix4D<double> > &DVV, UpdatableListMatrix < Matrix4D<double> > &DKK, UpdatableListMatrix < Matrix4D<double> > &DVK,
		UpdatableMatrix < Matrix4D<double> > &VP, UpdatableMatrix < Matrix4D<double> > &VL,
		UpdatableMatrix < Matrix4D<double> >  &G_local, UpdatableMatrix < Matrix4D<double> >  &G_radial,
		UpdatableListMatrix < Matrix4D<double> > &Sources, UpdatableListMatrix < Matrix4D<double> > &Losses, UpdatableListMatrix < Matrix4D<double> > &Losses_conv) {

	Parameters parameters("parameters.ini", argc, argv);

	parameters.getParameter("time_total", time_total);
	parameters.getParameter("time_step", time_step);
	parameters.getParameter("time_output", time_output);
	parameters.getParameter("max_threads", max_threads);

	parameters.getParameter("time_first", time_first);
	parameters.getParameter("it_first", it_first);

	parameters.getParameter("inversion_method", inversion_method);

	parameters.getParameter("include_boundary", include_boundary);
	parameters.getParameter("Vl_BC_from_convection", Vl_BC_from_convection);
	parameters.getParameter("Vu_BC_from_convection", Vu_BC_from_convection);

    parameters.getParameter("io_method", io_method);

	parameters.getParameter("input_folder", InputFolder);
	parameters.getParameter("output_folder", OutputFolder);

	parameters.getParameter("run_remapping", run_remapping);
	parameters.getParameter("run_convection", run_convection);
	parameters.getParameter("run_radial_diffusion", run_radial_diffusion);
	parameters.getParameter("run_local_diffusion", run_local_diffusion);
    parameters.getParameter("positive_PSD", positive_PSD);

    string initial_PSD = "PSD0";
    parameters.findParameter("initial_PSD", "PSD0") >> initial_PSD;

    // NOTE: this option is left only for backward compatibility of the code and is excessive
    string use_matlab = "false";
    parameters.getParameter("use_matlab", use_matlab);

    if (use_matlab == "true"){
        io_method = "matlab";
    }

#if !(MATLAB_CAPABLE)
{
	if (io_method.compare("matlab") == 0){
        printf("Error! Trying to use matlab files but the executable was compiled without matlab capabilities.");
        exit(EXIT_FAILURE);
	}
}
#endif

    // Opening up grid.plt in order to get P,R,V,K sizes which are stored in the header
    ifstream input;
    input.open((InputFolder + "grid.plt").c_str());
            //if (input == NULL) {
    // if the file is missing sends message to the logger
    if (!input.is_open()) {
        Logger::warning << "Grid file " << (InputFolder + "grid.plt") << " not found." << endl;
        return false;
    }

    // "I"=K_size, "J"=V_size, "K"=R_size, "L"=P_size from grid.plt from create_grid.m generated in matlab
    string inBuf;
    input >> inBuf;
    while (strcasecmp(inBuf.c_str(), "ZONE") != 0 && !input.eof() ) input >> inBuf;
    while(inBuf.substr(inBuf.size() - 1).compare("I") && !input.eof()) getline(input, inBuf, '=');
    input >> K_size;
    while(inBuf.substr(inBuf.size() - 1).compare("J") && !input.eof()) getline(input, inBuf, '=');
    input >> V_size;
    while(inBuf.substr(inBuf.size() - 1).compare("K") && !input.eof()) getline(input, inBuf, '=');
    input >> R_size;
    while(inBuf.substr(inBuf.size() - 1).compare("L") && !input.eof()) getline(input, inBuf, '=');
    input >> P_size;



    // If the file doesn't have data send error message
    if (input.eof()) {
        Logger::warning << "Grid file error." << endl;
        return false;
    }
    input.close();

	// sets the size of L and R to be the same - both measure distance
	L_size = R_size;
	// records the sizes of P,R,V,K to the logger
	Logger::message << "P_size = " << P_size << ", R_size = " << R_size << ", V_size = " << V_size << ", K_size = " << K_size << endl;


	///////
	// Allocate memory
	//////
	AllocateMemory(PSD,
			P, R, V, K, L,
			P_size, R_size, V_size, K_size,
			PSD_l_P, PSD_u_P, PSD_l_R, PSD_u_R, PSD_l_V, PSD_u_V, PSD_l_K, PSD_u_K, PSD_l_L, PSD_u_L,
			DLL, DVV, DVK, DKK,
			VP, VL,
			G_local, G_radial, Sources, Losses, Losses_conv);

    // Read grid from .plt file only
	P.readFromFile(InputFolder + "grid.plt",  1);
	R.readFromFile(InputFolder + "grid.plt",  2);
	V.readFromFile(InputFolder + "grid.plt",  3);
	K.readFromFile(InputFolder + "grid.plt",  4);

    // Read from Lstar file if Lstar.tab is not present
    if (!L.readFromIniFile(InputFolder + "Lstar.tab", P, R, V, K)){
        L.readFromAnyFile(InputFolder + "Lstar", io_method, P, R, V, K);
	}

	L.update(time_first, P, R, V, K); // Load L-star so it'll be available

    PSD.readFromAnyFile(InputFolder + initial_PSD, "ascii");

	// For all of the following load the .tab file if it exists, if not load the corresponding data file

    if (!DLL.readFromIniFile(InputFolder + "DLL.tab", P, L, V, K)){
        DLL.readFromAnyFile(InputFolder + "DLL", io_method, P, L, V, K);
	}

    if (!DVV.readFromIniFile(InputFolder + "DVV.tab", P, R, V, K)){
		DVV.readFromAnyFile(InputFolder + "DVV", io_method, P, R, V, K);
	}

    if (!DKK.readFromIniFile(InputFolder + "DKK.tab", P, R, V, K)){
		DKK.readFromAnyFile(InputFolder + "DKK", io_method, P, R, V, K);
	}

    if (!DVK.readFromIniFile(InputFolder + "DVK.tab", P, R, V, K)){
		DVK.readFromAnyFile(InputFolder + "DVK", io_method, P, R, V, K);
	}

    if (!VP.readFromIniFile(InputFolder + "VP.tab", P, R, V, K)){
		VP.readFromAnyFile(InputFolder + "VP", io_method, P, R, V, K);
	}

    if (!VL.readFromIniFile(InputFolder + "VR.tab", P, R, V, K)){
		VL.readFromAnyFile(InputFolder + "VR", io_method, P, R, V, K);
	}

    if (!G_local.readFromIniFile(InputFolder + "G_local.tab", P, R, V, K)){
		G_local.readFromAnyFile(InputFolder + "G_local", io_method, P, R, V, K);
	}

    if (!G_radial.readFromIniFile(InputFolder + "G_radial.tab", P, R, V, K)){
		G_radial.readFromAnyFile(InputFolder + "G_radial", io_method, P, R, V, K);
	}

	G_local.update(time_first, P, R, V, K);
	G_radial.update(time_first, P, R, V, K);

    if (!Sources.readFromIniFile(InputFolder + "Sources.tab", P, R, V, K)){
		Sources.readFromAnyFile(InputFolder + "Sources", io_method, P, R, V, K);
	}

    if (!Losses.readFromIniFile(InputFolder + "Losses.tab", P, R, V, K)){
        Losses.readFromAnyFile(InputFolder + "Losses_losscone", io_method, P, R, V, K);
	}

    if (!Losses_conv.readFromIniFile(InputFolder + "Losses_conv.tab", P, R, V, K))
    {
        Logger::error << "Cannot read file " << InputFolder << "Losses_conv.tab" << endl;
        exit(EXIT_FAILURE);
    }

	// BC conditions
	// Always use periodic boundary for Phi
	Pl_BC_type = "BCT_PERIODIC";
	Pu_BC_type = "BCT_PERIODIC";

	// default values
	Rl_BC_type = "BCT_CONSTANT_VALUE";
	Ru_BC_type = "BCT_CONSTANT_VALUE";
	Ll_BC_type = Rl_BC_type;
	Lu_BC_type = Ru_BC_type;

	Vl_BC_type = "BCT_CONSTANT_VALUE";
	Vu_BC_type = "BCT_CONSTANT_VALUE";

	Kl_BC_type = "BCT_CONSTANT_VALUE";
	//K_LBC_type = "BCT_CONSTANT_DERIVATIVE"; // - Careful with this - derivative is df/dK|const(V), which at alpha=0 is not the same as const(energy)
	//K_UBC_type = "BCT_CONSTANT_VALUE";
	Ku_BC_type = "BCT_CONSTANT_DERIVATIVE";

	// store all of the boundary conditions parameters from BC.ini
	Parameters BC_parameters(InputFolder + "BC.ini", argc, argv);

	// Read BC type fro BC.ini file
	BC_parameters.getParameter("Rl_BC_type", Rl_BC_type);
	BC_parameters.getParameter("Ru_BC_type", Ru_BC_type);
	Ll_BC_type = Rl_BC_type;
	Lu_BC_type = Ru_BC_type;
	BC_parameters.getParameter("Ll_BC_type", Ll_BC_type);
	BC_parameters.getParameter("Lu_BC_type", Lu_BC_type);

	BC_parameters.getParameter("Vl_BC_type", Vl_BC_type);
	BC_parameters.getParameter("Vu_BC_type", Vu_BC_type);
	BC_parameters.getParameter("Kl_BC_type", Kl_BC_type);
	BC_parameters.getParameter("Ku_BC_type", Ku_BC_type);

	// BC values
	PSD_l_P = -1e99; // Not in use, periodic boundary conditions
	PSD_u_P = -1e99; //

	// Read BC values from other files
	// For all of the following load the .tab file if it exists, if not load the corresponding data file to set the boundaries

    if (!PSD_l_R.readFromIniFile(InputFolder + "Rl_BC.tab",P.xSlice(0), V.xSlice(0), K.xSlice(0))) {
        PSD_l_R.readFromAnyFile(InputFolder + "Rl_BC", io_method, P.xSlice(0), V.xSlice(0), K.xSlice(0));
    }
    if (!PSD_u_R.readFromIniFile(InputFolder + "Ru_BC.tab", P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1))) {
        PSD_u_R.readFromAnyFile(InputFolder + "Ru_BC", io_method, P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1));
    }
    if (!PSD_l_L.readFromIniFile(InputFolder + "Ll_BC.tab", P.xSlice(0), V.xSlice(0), K.xSlice(0))) {
        PSD_l_L.readFromAnyFile(InputFolder + "Ll_BC", io_method, P.xSlice(0), V.xSlice(0), K.xSlice(0));
    }
    if (!PSD_u_L.readFromIniFile(InputFolder + "Lu_BC.tab", P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1))) {
        PSD_u_L.readFromAnyFile(InputFolder + "Lu_BC", io_method, P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1));
    }
    if (!PSD_l_V.readFromIniFile(InputFolder + "Vl_BC.tab", P.ySlice(0), R.ySlice(0), K.ySlice(0))) {
        PSD_l_V.readFromAnyFile(InputFolder + "Vl_BC", io_method, P.ySlice(0), R.ySlice(0), K.ySlice(0));
    }
    if (!PSD_u_V.readFromIniFile(InputFolder + "Vu_BC.tab", P.ySlice(P.size_y-1), R.ySlice(R.size_y-1), K.ySlice(K.size_y-1))) {
        PSD_u_V.readFromAnyFile(InputFolder + "Vu_BC", io_method, P.ySlice(P.size_y-1), R.ySlice(R.size_y-1), K.ySlice(K.size_y-1));
    }
    if (!PSD_l_K.readFromIniFile(InputFolder + "Kl_BC.tab", P.zSlice(0), R.zSlice(0), V.zSlice(0))) {
        PSD_l_K.readFromAnyFile(InputFolder + "Kl_BC", io_method, P.zSlice(0), R.zSlice(0), V.zSlice(0));
    }
    if (!PSD_u_K.readFromIniFile(InputFolder + "Ku_BC.tab", P.zSlice(P.size_z-1), R.zSlice(R.size_z-1), V.zSlice(V.size_z-1))) {
        PSD_u_K.readFromAnyFile(InputFolder + "Ku_BC", io_method, P.zSlice(P.size_z-1), R.zSlice(R.size_z-1), V.zSlice(V.size_z-1));
    }

	// Log the boundary types
	Logger::message << "Ll_BC = " << Ll_BC_type << "; Lu_BC = " << Lu_BC_type << ";" << endl;
	Logger::message << "Vl_BC = " << Vl_BC_type << "; Vu_BC = " << Vu_BC_type << ";" << endl;
	Logger::message << "Kl_BC = " << Kl_BC_type << "; Ku_BC = " << Ku_BC_type << ";" << endl;

	return true;
}
