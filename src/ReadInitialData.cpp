/*
 * ReadInitialData.cpp
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


void AllocateMemory(Matrix4D<double> &PSD,
		Matrix4D<double> &P, Matrix4D<double> &R, Matrix4D<double> &V,  Matrix4D<double> &K, Matrix4D<double> &L,
		int &P_size, int &R_size, int &V_size, int &K_size,
		Matrix3D<double> &P_LBC, Matrix3D<double> &P_UBC, Matrix3D<double> &R_LBC, Matrix3D<double> &R_UBC, Matrix3D<double> &V_LBC, Matrix3D<double> &V_UBC, Matrix3D<double> &K_LBC, Matrix3D<double> &K_UBC, Matrix3D<double> &L_LBC, Matrix3D<double> &L_UBC,
		Matrix4D<double> &DLL, Matrix4D<double> &DVV, Matrix4D<double> &DVK, Matrix4D<double> &DKK,
		Matrix4D<double> &VP, Matrix4D<double> &VL,
		Matrix4D<double> &G_local, Matrix4D<double> &G_radial, Matrix4D<double> &Sources, Matrix4D<double> &Losses) {


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
bool ReadInitialData(string InputFolder, int argc, char* argv[],
		double &time_total, double &time_step, double &time_output, double &time_first, long int &it_first, int &max_threads,
		string &inversion_method,
		Matrix4D<double> &PSD,
		Matrix4D<double> &P, Matrix4D<double> &R, Matrix4D<double> &V, Matrix4D<double> &K, UpdatableMatrix < Matrix4D<double> > &L,
		int &P_size, int &R_size, int &V_size, int &K_size, int &L_size,
		Matrix3D<double> &PSD_l_P, Matrix3D<double> &PSD_u_P, UpdatableMatrix< Matrix3D<double> > &PSD_l_R, UpdatableMatrix< Matrix3D<double> > &PSD_u_R, UpdatableMatrix< Matrix3D<double> > &PSD_l_V, UpdatableMatrix< Matrix3D<double> >&PSD_u_V, UpdatableMatrix< Matrix3D<double> > &PSD_l_K, UpdatableMatrix< Matrix3D<double> > &PSD_u_K, UpdatableMatrix< Matrix3D<double> > &PSD_l_L, UpdatableMatrix< Matrix3D<double> > &PSD_u_L,
		string &Pl_BC_type, string &Pu_BC_type, string &Rl_BC_type, string &Ru_BC_type, string &Vl_BC_type, string &Vu_BC_type, string &Kl_BC_type, string &Ku_BC_type, string &Ll_BC_type, string &Lu_BC_type,
		UpdatableListMatrix < Matrix4D<double> > &DLL, UpdatableListMatrix < Matrix4D<double> > &DVV, UpdatableListMatrix < Matrix4D<double> > &DKK, UpdatableListMatrix < Matrix4D<double> > &DVK,
		UpdatableMatrix < Matrix4D<double> > &VP, UpdatableMatrix < Matrix4D<double> > &VL,
		UpdatableMatrix < Matrix4D<double> >  &G_local, UpdatableMatrix < Matrix4D<double> >  &G_radial,
		UpdatableListMatrix < Matrix4D<double> > &Sources, UpdatableListMatrix < Matrix4D<double> > &Losses) {

	Parameters parameters(InputFolder + "parameters.ini", argc, argv);

	parameters.getParameter("time_total", time_total);
	parameters.getParameter("time_step", time_step);
	parameters.getParameter("time_output", time_output);
	parameters.getParameter("max_threads", max_threads);

	parameters.getParameter("time_first", time_first);
	parameters.getParameter("it_first", it_first);

	parameters.getParameter("inversion_method", inversion_method);

	//parameters.findParameter("total_time", "4") >> total_time;
	//parameters.findParameter("dt", "0.047") >> dt;
	//parameters.findParameter("output_time", "0.047") >> output_time;
	//parameters.findParameter("max_threads", "8") >> max_threads;

	ifstream input;
	input.open((InputFolder + "grid.dat").c_str());
	//if (input == NULL) {
	if (!input.is_open()) {
		cout << "Grid file " << (InputFolder + "grid.dat") << " not found." << endl;
		abort();
	}
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

	if (input.eof()) {
		cout << "Grid file error." << endl;
		abort();
	}

	input.close();

	L_size = R_size;

	cout << "P_size = " << P_size  << ", R_size = " << R_size << ", V_size = " << V_size << ", K_size = " << K_size << endl;


	///////
	// Allocate memory
	//////
	AllocateMemory(PSD,
			P, R, V, K, L,
			P_size, R_size, V_size, K_size,
			PSD_l_P, PSD_u_P, PSD_l_R, PSD_u_R, PSD_l_V, PSD_u_V, PSD_l_K, PSD_u_K, PSD_l_L, PSD_u_L,
			DLL, DVV, DVK, DKK,
			VP, VL,
			G_local, G_radial, Sources, Losses);


	// Read values
	P.readFromFile(InputFolder + "grid.dat", 1);
	R.readFromFile(InputFolder + "grid.dat", 2);
	V.readFromFile(InputFolder + "grid.dat", 3);
	K.readFromFile(InputFolder + "grid.dat", 4);

	if (!L.readFromIniFile(InputFolder + "Lstar.ini", P, R, V, K))
		L.readFromFile(InputFolder + "Lstar.dat", P, R, V, K);
	L.update(time_first, P, R, V, K); // Load L-star so it'll available

	string initial_PSD = "PSD0.dat";
	//parameters.getParameter("initial_PSD", initial_PSD);
	parameters.findParameter("initial_PSD", "PSD0.dat") >> initial_PSD;

	//PSD.readFromFile(InputFolder + "PSD0.dat", P, R, I, K);
	PSD.readFromFile(InputFolder + initial_PSD);

	//DLL.readFromFile(InputFolder + "DLL.dat", P, R, V, K);
	//DVV.readFromFile(InputFolder + "DVV.dat", P, R, V, K);
	//DKK.readFromFile(InputFolder + "DKK.dat", P, R, V, K);
	//DVK.readFromFile(InputFolder + "DVK.dat", P, R, V, K);

	// DLL is defined on PLVK?
	if (!DLL.readFromIniFile(InputFolder + "DLL.ini", P, L, V, K))
		DLL.readFromFile(InputFolder + "DLL.dat", P, L, V, K);

	if (!DVV.readFromIniFile(InputFolder + "DVV.ini", P, R, V, K))
		DVV.readFromFile(InputFolder + "DVV.dat", P, R, V, K);

	if (!DKK.readFromIniFile(InputFolder + "DKK.ini", P, R, V, K))
		DKK.readFromFile(InputFolder + "DKK.dat", P, R, V, K);

	if (!DVK.readFromIniFile(InputFolder + "DVK.ini", P, R, V, K))
		DVK.readFromFile(InputFolder + "DVK.dat", P, R, V, K);

	if (!VP.readFromIniFile(InputFolder + "VP.ini", P, R, V, K))
		VP.readFromFile(InputFolder + "VP.dat", P, R, V, K);

	if (!VL.readFromIniFile(InputFolder + "VR.ini", P, R, V, K))
		VL.readFromFile(InputFolder + "VR.dat", P, R, V, K);

	if (!G_local.readFromIniFile(InputFolder + "G_local.ini", P, R, V, K))
		G_local.readFromFile(InputFolder + "G_local.dat", P, R, V, K);
	if (!G_radial.readFromIniFile(InputFolder + "G_radial.ini", P, R, V, K))
		G_radial.readFromFile(InputFolder + "G_radial.dat", P, R, V, K);

	// XXX: Have to load Jacobians here, cause right now we update them only when L is updated.
	// XXX: L is updated here so it will not be updated lated and Jacobians will not be loaded
	G_local.update(time_first, P, R, V, K);
	G_radial.update(time_first, P, R, V, K);

	// Try to read from ini-file
	if (!Sources.readFromIniFile(InputFolder + "Sources.ini", P, R, V, K))
		// if failed - try to read from just data-file
		Sources.readFromFile(InputFolder + "Sources.dat", P, R, V, K);

	// Try to read from ini-file
	if (!Losses.readFromIniFile(InputFolder + "Losses.ini", P, R, V, K))
		// if failed - try to read from just data-file
		Losses.readFromFile(InputFolder + "Losses.dat", P, R, V, K);

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
	if (!PSD_l_R.readFromIniFile(InputFolder + "Rl_BC.ini", P.xSlice(0), V.xSlice(0), K.xSlice(0))) {
		PSD_l_R.readFromFile(InputFolder + "Rl_BC.dat", P.xSlice(0), V.xSlice(0), K.xSlice(0));
	}
	if (!PSD_u_R.readFromIniFile(InputFolder + "Ru_BC.ini", P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1))) {
		PSD_u_R.readFromFile(InputFolder + "Ru_BC.dat",     P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1));
	}

	if (!PSD_l_L.readFromIniFile(InputFolder + "Ll_BC.ini", P.xSlice(0), V.xSlice(0), K.xSlice(0))) {
		PSD_l_L.readFromFile(InputFolder + "Ll_BC.dat",     P.xSlice(0), V.xSlice(0), K.xSlice(0));
	}
	if (!PSD_u_L.readFromIniFile(InputFolder + "Lu_BC.ini", P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1))) {
		PSD_u_L.readFromFile(InputFolder + "Lu_BC.dat",     P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1));
	}

	if (!PSD_l_V.readFromIniFile(InputFolder + "Vl_BC.ini", P.ySlice(0), R.ySlice(0), K.ySlice(0))) {
		PSD_l_V.readFromFile(InputFolder + "Vl_BC.dat",     P.ySlice(0), R.ySlice(0), K.ySlice(0));
	}
	if (!PSD_u_V.readFromIniFile(InputFolder + "Vu_BC.ini", P.ySlice(P.size_y-1), R.ySlice(R.size_y-1), K.ySlice(K.size_y-1))) {
		PSD_u_V.readFromFile(InputFolder + "Vu_BC.dat",     P.ySlice(P.size_y-1), R.ySlice(R.size_y-1), K.ySlice(K.size_y-1));
	}

	if (!PSD_l_K.readFromIniFile(InputFolder + "Kl_BC.ini", P.zSlice(0), R.zSlice(0), V.zSlice(0))) {
		PSD_l_K.readFromFile(InputFolder + "Kl_BC.dat",     P.zSlice(0), R.zSlice(0), V.zSlice(0));
	}
	if (!PSD_u_K.readFromIniFile(InputFolder + "Ku_BC.ini", P.zSlice(P.size_z-1), R.zSlice(R.size_z-1), V.zSlice(V.size_z-1))) {
		PSD_u_K.readFromFile(InputFolder + "Ku_BC.dat",     P.zSlice(P.size_z-1), R.zSlice(R.size_z-1), V.zSlice(V.size_z-1));
	}

	cout << "Ll_BC = " << Ll_BC_type << "; Lu_BC = " << Lu_BC_type << ";" << endl;
	cout << "Vl_BC = " << Vl_BC_type << "; Vu_BC = " << Vu_BC_type << ";" << endl;
	cout << "Kl_BC = " << Kl_BC_type << "; Ku_BC = " << Ku_BC_type << ";" << endl;

	return true;
}
