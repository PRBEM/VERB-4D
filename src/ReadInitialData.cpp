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
		string &inversion_method, string &use_matlab, string &include_boundary, 
		Matrix4D<double> &PSD,
		Matrix4D<double> &P, Matrix4D<double> &R, Matrix4D<double> &V, Matrix4D<double> &K, UpdatableMatrix < Matrix4D<double> > &L,
		int &P_size, int &R_size, int &V_size, int &K_size, int &L_size,
		Matrix3D<double> &PSD_l_P, Matrix3D<double> &PSD_u_P, UpdatableMatrix< Matrix3D<double> > &PSD_l_R, UpdatableMatrix< Matrix3D<double> > &PSD_u_R, UpdatableMatrix< Matrix3D<double> > &PSD_l_V, UpdatableMatrix< Matrix3D<double> >&PSD_u_V, UpdatableMatrix< Matrix3D<double> > &PSD_l_K, UpdatableMatrix< Matrix3D<double> > &PSD_u_K, UpdatableMatrix< Matrix3D<double> > &PSD_l_L, UpdatableMatrix< Matrix3D<double> > &PSD_u_L,
		string &Pl_BC_type, string &Pu_BC_type, string &Rl_BC_type, string &Ru_BC_type, string &Vl_BC_type, string &Vu_BC_type, string &Kl_BC_type, string &Ku_BC_type, string &Ll_BC_type, string &Lu_BC_type,
		UpdatableListMatrix < Matrix4D<double> > &DLL, UpdatableListMatrix < Matrix4D<double> > &DVV, UpdatableListMatrix < Matrix4D<double> > &DKK, UpdatableListMatrix < Matrix4D<double> > &DVK,
		UpdatableMatrix < Matrix4D<double> > &VP, UpdatableMatrix < Matrix4D<double> > &VL,
		UpdatableMatrix < Matrix4D<double> >  &G_local, UpdatableMatrix < Matrix4D<double> >  &G_radial,
		UpdatableListMatrix < Matrix4D<double> > &Sources, UpdatableListMatrix < Matrix4D<double> > &Losses) {

	Parameters parameters("parameters.ini", argc, argv);

	parameters.getParameter("time_total", time_total);
	parameters.getParameter("time_step", time_step);
	parameters.getParameter("time_output", time_output);
	parameters.getParameter("max_threads", max_threads);

	parameters.getParameter("time_first", time_first);
	parameters.getParameter("it_first", it_first);

	parameters.getParameter("inversion_method", inversion_method);
	
	// ADDED	
	parameters.getParameter("use_matlab", use_matlab);
	parameters.getParameter("include_boundary", include_boundary);
	// END ADDED
	
	
	parameters.getParameter("input_folder", InputFolder);
	parameters.getParameter("output_folder", OutputFolder);

	FILE *file;
	//if (file = fopen(name.c_str(), "r"))
	
	
	if ( (!(file = fopen((InputFolder + "grid.mat").c_str() , "r"))) || use_matlab == "false")
	{
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

	}
	// ADDED
	else // IF MATLAB SELECTED
	{
		MATFile *mfPtr; /* MAT-file pointer */
   		mxArray *aPtr;  /* mxArray pointer */
		mwSize nElements;       /* number of elements in array */
    	mwIndex eIdx;           /* element index */
    	const mxArray *fPtr;    /* field pointer */
		string arr; 			/* variable name*/
		
		
		
		mfPtr = matOpen((InputFolder + "grid.mat").c_str(), "r");
    	if (mfPtr == NULL) {
        	printf("Error opening file grid.mat\n");
    	}
		
    	string field = "size";
		
	 	arr = "P";
    	aPtr = matGetVariable(mfPtr, arr.c_str());
    	if (aPtr == NULL) {
        	printf("mxArray not found: %s\n", arr.c_str());
   		}   
		   
		if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
        	if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
            	printf("Field not found: %s\n", field.c_str());
        	}
        	else {
       		 		nElements = (mwSize)mxGetNumberOfElements(aPtr);
   					 for (eIdx = 0; eIdx < nElements; eIdx++) {
						fPtr = mxGetField(aPtr, eIdx, field.c_str());
						P_size = int(mxGetScalar(fPtr));
					}
   				 }
	        }
    	
		arr = "V";
    	aPtr = matGetVariable(mfPtr, arr.c_str());
    	if (aPtr == NULL) {
        	arr = "pc";
			aPtr = matGetVariable(mfPtr, arr.c_str());
    		if (aPtr == NULL) {
				printf("mxArray not found: %s\n", arr.c_str());
			}
   		}   
		   
		if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
        	if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
            	printf("Field not found: %s\n", field.c_str());
        	}
        	else 
			{
				nElements = (mwSize)mxGetNumberOfElements(aPtr);
   			 	for (eIdx = 0; eIdx < nElements; eIdx++) {
					fPtr = mxGetField(aPtr, eIdx, field.c_str());
       		 		V_size = int(mxGetScalar(fPtr));
				}
			}
    	}
		arr = "R";
    	aPtr = matGetVariable(mfPtr, arr.c_str());
    	if (aPtr == NULL) {
        	arr = "L";
			aPtr = matGetVariable(mfPtr, arr.c_str());
    		if (aPtr == NULL) {
				printf("mxArray not found: %s\n", arr.c_str());
			}
   		}   
		   
		if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
        	if (mxGetFieldNumber(aPtr, field.c_str()) == -1) {
            	printf("Field not found: %s\n", field.c_str());
        	}
        	else {
    				nElements = (mwSize)mxGetNumberOfElements(aPtr);
   					 for (eIdx = 0; eIdx < nElements; eIdx++) {
       		 			fPtr = mxGetField(aPtr, eIdx, field.c_str()); 
						R_size = int(mxGetScalar(fPtr));
					}
       		 		
   				 }
	        }
    	
		arr = "K";
    	aPtr = matGetVariable(mfPtr, arr.c_str());
    	if (aPtr == NULL) {
        	arr = "al";
			aPtr = matGetVariable(mfPtr, arr.c_str());
    		if (aPtr == NULL) {
				printf("mxArray not found: %s\n", arr.c_str());
			}
   		}   
		   
		if (mxGetClassID(aPtr) == mxSTRUCT_CLASS) {
        	if ( mxGetFieldNumber(aPtr, field.c_str()) == -1) {
            	printf("Field not found: %s\n", field.c_str());
        	}
        	else {
       		 		nElements = (mwSize)mxGetNumberOfElements(aPtr);
   				 	for (eIdx = 0; eIdx < nElements; eIdx++) {
						fPtr = mxGetField(aPtr, eIdx, field.c_str());
						K_size = int(mxGetScalar(fPtr));
					}
   				 }
	        }	
			
			mxDestroyArray(aPtr);
    		matClose(mfPtr);  	   
	}
	fclose(file);
	
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
			G_local, G_radial, Sources, Losses);

	// ADDED
	if (use_matlab == "false")
	{
	// Original Section	
	// Read values from grid.plt - last argument is the column being read in
	P.readFromFile(InputFolder + "grid.plt", 1);
	R.readFromFile(InputFolder + "grid.plt", 2);
	V.readFromFile(InputFolder + "grid.plt", 3);
	K.readFromFile(InputFolder + "grid.plt", 4);
	
	// End original section
	}
	else // if Matlab
	{
		if (!(file = fopen((InputFolder + "grid.mat").c_str() , "r"))) // Check to see if .mat file is there
		{
			P.readFromFile(InputFolder + "grid.plt", 1);
			R.readFromFile(InputFolder + "grid.plt", 2);
			V.readFromFile(InputFolder + "grid.plt", 3);
			K.readFromFile(InputFolder + "grid.plt", 4);
		}
		else
		{
			P.readFromMatlabFile(InputFolder + "grid.mat", 1);
			R.readFromMatlabFile(InputFolder + "grid.mat", 2);
			V.readFromMatlabFile(InputFolder + "grid.mat", 3);
			K.readFromMatlabFile(InputFolder + "grid.mat", 4);
		}
		fclose(file);
	}
	

	
	// ADDED
	// If matlab not selected
	if (use_matlab == "false")
	{
		// Read in Lstar.plt if Lstar.tab is not present
		if (!L.readFromIniFile(InputFolder + "Lstar.tab", P, R, V, K))
			L.readFromFile(InputFolder + "Lstar.plt", P, R, V, K);
	}
	else // if matlab
	{
		if (!(file = fopen((InputFolder + "Lstar.mat").c_str(), "r" ))) // Check to see if .mat file is there 
		{
			if (!L.readFromIniFile(InputFolder + "Lstar.tab", P, R, V, K))
			L.readFromFile(InputFolder + "Lstar.plt", P, R, V, K);
		}
		else
		{
			L.readFromMatlabFile(InputFolder + "Lstar.mat", P, R, V, K);
		}
		fclose(file);
	}

		
	
	L.update(time_first, P, R, V, K); // Load L-star so it'll be available
	
	string initial_PSD = "PSD0.plt";	
	
	// if (use_matlab == "false")
	// {
		parameters.findParameter("initial_PSD", "PSD0.plt") >> initial_PSD;
		PSD.readFromFile(InputFolder + initial_PSD);
	// }
	// else
	// {
	// 	parameters.findParameter("initial_PSD", "PSD0.mat") >> initial_PSD;
	// 	PSD.readFromMatlabFile(InputFolder + initial_PSD);
	// }
	

	// For all of the following load the .tab file if it exists, if not load the corresponding .plt file
	// DLL is defined on PLVK?
	
	
	// ADDED
	// If matlab not selected
	if (use_matlab == "false")
	{
		if (!DLL.readFromIniFile(InputFolder + "DLL.tab", P, L, V, K))
		DLL.readFromFile(InputFolder + "DLL.plt", P, L, V, K);
	}
	else // if matlab
	{
		if (!(file = fopen((InputFolder + "DLL.mat").c_str() , "r"))) // Check to see if .mat file is there 
		{
				if (!DLL.readFromIniFile(InputFolder + "DLL.tab", P, L, V, K))
				DLL.readFromFile(InputFolder + "DLL.plt", P, L, V, K);
		}
		else
		{
			DLL.readFromMatlabFile(InputFolder + "DLL.mat", P, L, V, K);
		}
	}
	fclose(file);


	// ADDED
	// If matlab not selected
	if (use_matlab == "false")
	{
		if (!DVV.readFromIniFile(InputFolder + "DVV.tab", P, R, V, K))
		DVV.readFromFile(InputFolder + "DVV.plt", P, R, V, K);
	}
	else // if matlab
	{
		if (!(file = fopen((InputFolder + "DVV.mat").c_str() , "r"))) // Check to see if .mat file is there 
		{
				if (!DVV.readFromIniFile(InputFolder + "DVV.tab", P, R, V, K))
				DVV.readFromFile(InputFolder + "DVV.plt", P, R, V, K);
		}
		else
		{
			DVV.readFromMatlabFile(InputFolder + "DVV.mat", P, R, V, K);
		}
	}
	fclose(file);


	// ADDED
	// If matlab not selected
	if (use_matlab == "false")
	{
		if (!DKK.readFromIniFile(InputFolder + "DKK.tab", P, R, V, K))
		DKK.readFromFile(InputFolder + "DKK.plt", P, R, V, K);
	}
	else // if matlab
	{
		if (!(file = fopen((InputFolder + "DKK.mat").c_str() , "r"))) // Check to see if .mat file is there 
		{
				if (!DKK.readFromIniFile(InputFolder + "DKK.tab", P, R, V, K))
				DKK.readFromFile(InputFolder + "DKK.plt", P, R, V, K);
		}
		else
		{
			DKK.readFromMatlabFile(InputFolder + "DKK.mat", P, R, V, K);
		}
	}
	fclose(file);

	// ADDED
	// If matlab not selected
	if (use_matlab == "false")
	{
		if (!DVK.readFromIniFile(InputFolder + "DVK.tab", P, R, V, K))
		DVK.readFromFile(InputFolder + "DVK.plt", P, R, V, K);
	}
	else // if matlab
	{
		if (!(file = fopen((InputFolder + "DVK.mat").c_str() , "r"))) // Check to see if .mat file is there  
		{
				if (!DVK.readFromIniFile(InputFolder + "DVK.tab", P, R, V, K))
				DVK.readFromFile(InputFolder + "DVK.plt", P, R, V, K);
		}
		else
		{
			DVK.readFromMatlabFile(InputFolder + "DVK.mat", P, R, V, K);
		}
	}
	fclose(file);
	
	// ADDED
	// If matlab not selected
	if (use_matlab == "false")
	{
		if (!VP.readFromIniFile(InputFolder + "VP.tab", P, R, V, K))
		VP.readFromFile(InputFolder + "VP.plt", P, R, V, K);
	}
	else // if matlab
	{
		if (!(file = fopen((InputFolder + "VP.mat").c_str() , "r"))) // Check to see if .mat file is there 
		{
				if (!VP.readFromIniFile(InputFolder + "VP.tab", P, R, V, K))
				VP.readFromFile(InputFolder + "VP.plt", P, R, V, K);
		}
		else
		{
			VP.readFromMatlabFile(InputFolder + "VP.mat", P, R, V, K);
		}
	}
	fclose(file);
		
	// ADDED
	// If matlab not selected
	if (use_matlab == "false")
	{
		if (!VL.readFromIniFile(InputFolder + "VR.tab", P, R, V, K))
		VL.readFromFile(InputFolder + "VR.plt", P, R, V, K);
	}
	else // if matlab
	{
		if (!(file = fopen((InputFolder + "VR.mat").c_str() , "r"))) // Check to see if .mat file is there 
		{
				if (!VL.readFromIniFile(InputFolder + "VR.tab", P, R, V, K))
				VL.readFromFile(InputFolder + "VR.plt", P, R, V, K);
		}
		else
		{
			VL.readFromMatlabFile(InputFolder + "VR.mat", P, R, V, K);
		}
	}
	fclose(file);	
	
	// ADDED
	// If matlab not selected
	if (use_matlab == "false")
	{
		if (!G_local.readFromIniFile(InputFolder + "G_local.tab", P, R, V, K))
		G_local.readFromFile(InputFolder + "G_local.plt", P, R, V, K);
	}
	else // if matlab
	{
		if (!(file = fopen((InputFolder + "G_local.mat").c_str() , "r"))) // Check to see if .mat file is there 
		{
				if (!G_local.readFromIniFile(InputFolder + "G_local.tab", P, R, V, K))
				G_local.readFromFile(InputFolder + "G_local.plt", P, R, V, K);
		}
		else
		{
			G_local.readFromMatlabFile(InputFolder + "G_local.mat", P, R, V, K);
		}
	}
	fclose(file);
	
	// ADDED
	// If matlab not selected
	if (use_matlab == "false")
	{
		if (!G_radial.readFromIniFile(InputFolder + "G_radial.tab", P, R, V, K))
		G_radial.readFromFile(InputFolder + "G_radial.plt", P, R, V, K);
	}
	else // if matlab
	{
		if (!(file = fopen((InputFolder + "G_radial.mat").c_str() , "r"))) // Check to see if .mat file is there 
		{
				if (!G_radial.readFromIniFile(InputFolder + "G_radial.tab", P, R, V, K))
				G_radial.readFromFile(InputFolder + "G_radial.plt", P, R, V, K);
		}
		else
		{
			G_radial.readFromMatlabFile(InputFolder + "G_radial.mat", P, R, V, K);
		}
	}
	fclose(file);	

	// if (!VL.readFromIniFile(InputFolder + "VR.tab", P, R, V, K))
	// 	VL.readFromFile(InputFolder + "VR.plt", P, R, V, K);

	// if (!G_local.readFromIniFile(InputFolder + "G_local.tab", P, R, V, K))
	// 	G_local.readFromFile(InputFolder + "G_local.plt", P, R, V, K);
	// if (!G_radial.readFromIniFile(InputFolder + "G_radial.tab", P, R, V, K))
	// 	G_radial.readFromFile(InputFolder + "G_radial.plt", P, R, V, K);

	// XXX: Have to load Jacobians here, cause right now we update them only when L is updated.
	// XXX: L is updated here so it will not be updated later and Jacobians will not be loaded
	G_local.update(time_first, P, R, V, K);
	G_radial.update(time_first, P, R, V, K);


	// ADDED
	// If matlab not selected
	if (use_matlab == "false")
	{
		if (!Sources.readFromIniFile(InputFolder + "Sources.tab", P, R, V, K))
		Sources.readFromFile(InputFolder + "Sources.plt", P, R, V, K);
	}
	else // if matlab
	{
		if (!(file = fopen((InputFolder + "Sources.mat").c_str() , "r"))) // Check to see if .mat file is there 
		{
				if (!Sources.readFromIniFile(InputFolder + "Sources.tab", P, R, V, K))
				Sources.readFromFile(InputFolder + "Sources.plt", P, R, V, K);
		}
		else
		{
			Sources.readFromMatlabFile(InputFolder + "Sources.mat", P, R, V, K);
		}
	}
	fclose(file);	
	
	// ADDED
	// If matlab not selected
	// if (use_matlab == "false")
	// {
	// 	if (!Losses.readFromIniFile(InputFolder + "Losses.tab", P, R, V, K))
	// 	Losses.readFromFile(InputFolder + "Losses_losscone.plt", P, R, V, K);
	// }
	// else // if matlab
	// {
	// 	inputcheck.open((InputFolder + "Losses_losscone.mat").c_str()); // Check to see if .mat file is there 
	// 	if (!inputcheck.is_open()) 
	// 	{
	// 			if (!Losses.readFromIniFile(InputFolder + "Losses.tab", P, R, V, K))
	// 			Losses.readFromFile(InputFolder + "Losses_losscone.plt", P, R, V, K);
	// 	}
	// 	else
	// 	{
	// 		Losses.readFromMatlabFile(InputFolder + "Losses_losscone.mat", P, R, V, K);
	// 	}
	// }
	// inputcheck.close();	
	
	
	
	
	// ADDED
	// If matlab selected
	if (!Losses.readFromIniFile(InputFolder + "Losses.tab", P, R, V, K))
	{
		if (fopen((InputFolder + "Losses_losscone.mat").c_str() , "r") && use_matlab == "true" )
		{
			Losses.readFromMatlabFile(InputFolder + "Losses_losscone.mat", P, R, V, K);
		}
		else
		{
			Losses.readFromFile(InputFolder + "Losses_losscone.plt", P, R, V, K);
		}
		fclose(file);	
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
	// For all of the following load the .tab file if it exists, if not load the corresponding .plt file to set the boundaries
	// if (!PSD_l_R.readFromIniFile(InputFolder + "Rl_BC.tab", P.xSlice(0), V.xSlice(0), K.xSlice(0))) {
	// 	PSD_l_R.readFromFile(InputFolder + "Rl_BC.plt", P.xSlice(0), V.xSlice(0), K.xSlice(0));
	// }
	
	
	// ADDED	
	// If matlab not selected
	if (use_matlab == "false")
	{
		if (!PSD_l_R.readFromIniFile(InputFolder + "Rl_BC.tab",P.xSlice(0), V.xSlice(0), K.xSlice(0))) {
		PSD_l_R.readFromFile(InputFolder + "Rl_BC.plt", P.xSlice(0), V.xSlice(0), K.xSlice(0));
		}
		if (!PSD_u_R.readFromIniFile(InputFolder + "Ru_BC.tab", P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1))) {
		PSD_u_R.readFromFile(InputFolder + "Ru_BC.plt",     P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1));
		}
		if (!PSD_l_L.readFromIniFile(InputFolder + "Ll_BC.tab", P.xSlice(0), V.xSlice(0), K.xSlice(0))) {
			PSD_l_L.readFromFile(InputFolder + "Ll_BC.plt",     P.xSlice(0), V.xSlice(0), K.xSlice(0));
		}
		if (!PSD_u_L.readFromIniFile(InputFolder + "Lu_BC.tab", P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1))) {
			PSD_u_L.readFromFile(InputFolder + "Lu_BC.plt",     P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1));
		}
		if (!PSD_l_V.readFromIniFile(InputFolder + "Vl_BC.tab", P.ySlice(0), R.ySlice(0), K.ySlice(0))) {
			PSD_l_V.readFromFile(InputFolder + "Vl_BC.plt",     P.ySlice(0), R.ySlice(0), K.ySlice(0));
		}
		if (!PSD_u_V.readFromIniFile(InputFolder + "Vu_BC.tab", P.ySlice(P.size_y-1), R.ySlice(R.size_y-1), K.ySlice(K.size_y-1))) {
			PSD_u_V.readFromFile(InputFolder + "Vu_BC.plt",     P.ySlice(P.size_y-1), R.ySlice(R.size_y-1), K.ySlice(K.size_y-1));
		}
		if (!PSD_l_K.readFromIniFile(InputFolder + "Kl_BC.tab", P.zSlice(0), R.zSlice(0), V.zSlice(0))) {
			PSD_l_K.readFromFile(InputFolder + "Kl_BC.plt",     P.zSlice(0), R.zSlice(0), V.zSlice(0));
		}
		if (!PSD_u_K.readFromIniFile(InputFolder + "Ku_BC.tab", P.zSlice(P.size_z-1), R.zSlice(R.size_z-1), V.zSlice(V.size_z-1))) {
			PSD_u_K.readFromFile(InputFolder + "Ku_BC.plt",     P.zSlice(P.size_z-1), R.zSlice(R.size_z-1), V.zSlice(V.size_z-1));
		}
	}	
	else // if matlab
	{
		
		// R boundary conditions
		
		
		if (!(file = fopen((InputFolder + "Rl_BC.mat").c_str() , "r"))) // Check to see if .mat file is there 
		{
				if (!PSD_l_R.readFromIniFile(InputFolder + "Rl_BC.tab", P.xSlice(0), V.xSlice(0), K.xSlice(0)))
				PSD_l_R.readFromFile(InputFolder + "Rl_BC.plt", P.xSlice(0), V.xSlice(0), K.xSlice(0));
		}
		else
		{
			PSD_l_R.readFromMatlabFile(InputFolder + "Rl_BC.mat", P.xSlice(0), V.xSlice(0), K.xSlice(0));
		}
		fclose(file);	
		if (!(file = fopen((InputFolder + "Ru_BC.mat").c_str() , "r"))) // Check to see if .mat file is there 
		 
		{
				if (!PSD_u_R.readFromIniFile(InputFolder + "Ru_BC.tab", P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1)))
				PSD_u_R.readFromFile(InputFolder + "Ru_BC.plt", P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1));
		}
		else
		{
			PSD_u_R.readFromMatlabFile(InputFolder + "Ru_BC.mat",P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1));
		}
		fclose(file);
		
		// L boundary conditions
		
		if (!(file = fopen((InputFolder + "Ll_BC.mat").c_str() , "r"))) // Check to see if .mat file is there 
		{
				if (!PSD_l_L.readFromIniFile(InputFolder + "Ll_BC.tab", P.xSlice(0), V.xSlice(0), K.xSlice(0)))
				PSD_l_L.readFromFile(InputFolder + "Ll_BC.plt", P.xSlice(0), V.xSlice(0), K.xSlice(0));
		}
		else
		{
			PSD_l_L.readFromMatlabFile(InputFolder + "Ll_BC.mat", P.xSlice(0), V.xSlice(0), K.xSlice(0));
		}
		fclose(file);
		if (!(file = fopen((InputFolder + "Lu_BC.mat").c_str() , "r"))) // Check to see if .mat file is there  
		{
				if (!PSD_u_L.readFromIniFile(InputFolder + "Lu_BC.tab", P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1)))
				PSD_u_L.readFromFile(InputFolder + "Lu_BC.plt", P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1));
		}
		else
		{
			PSD_u_L.readFromMatlabFile(InputFolder + "Lu_BC.mat",P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1));
		}
		fclose(file);
		
		
		
		// V boundary conditions
		
		if (!(file = fopen((InputFolder + "Vl_BC.mat").c_str() , "r"))) // Check to see if .mat file is there 
		{
				if (!PSD_l_V.readFromIniFile(InputFolder + "Vl_BC.tab", P.ySlice(0), R.ySlice(0), K.ySlice(0)))
				PSD_l_V.readFromFile(InputFolder + "Vl_BC.plt", P.ySlice(0), R.ySlice(0), K.ySlice(0));
		}
		else
		{
			PSD_l_V.readFromMatlabFile(InputFolder + "Vl_BC.mat", P.ySlice(0), R.ySlice(0), K.ySlice(0));
		}
		fclose(file);
		if (!(file = fopen((InputFolder + "Vu_BC.mat").c_str() , "r"))) // Check to see if .mat file is there 
		{
				if (!PSD_u_V.readFromIniFile(InputFolder + "Vu_BC.tab", P.ySlice(P.size_y-1), R.ySlice(R.size_y-1), K.ySlice(K.size_y-1)))
				PSD_u_V.readFromFile(InputFolder + "Vu_BC.plt", P.ySlice(P.size_y-1), R.ySlice(R.size_y-1), K.ySlice(K.size_y-1));
		}
		else
		{
			PSD_u_V.readFromMatlabFile(InputFolder + "Vu_BC.mat",P.ySlice(P.size_y-1), R.ySlice(R.size_y-1), K.ySlice(K.size_y-1));
		}
		fclose(file);
		
		
		
		// K boundary conditions
		
		if (!(file = fopen((InputFolder + "Kl_BC.mat").c_str() , "r"))) // Check to see if .mat file is there 
		{
				if (!PSD_l_K.readFromIniFile(InputFolder + "Kl_BC.tab", P.zSlice(0), R.zSlice(0), V.zSlice(0)))
				PSD_l_K.readFromFile(InputFolder + "Kl_BC.plt", P.zSlice(0), R.zSlice(0), V.zSlice(0));
		}
		else
		{
			PSD_l_K.readFromMatlabFile(InputFolder + "Kl_BC.mat", P.zSlice(0), R.zSlice(0), V.zSlice(0));
		}
		fclose(file);
		if (!(file = fopen((InputFolder + "Ku_BC.mat").c_str() , "r"))) // Check to see if .mat file is there 
		{
				if (!PSD_u_K.readFromIniFile(InputFolder + "Ku_BC.tab", P.zSlice(P.size_z-1), R.zSlice(R.size_z-1), V.zSlice(V.size_z-1)))
				PSD_u_K.readFromFile(InputFolder + "Ku_BC.plt", P.zSlice(P.size_z-1), R.zSlice(R.size_z-1), V.zSlice(V.size_z-1));
		}
		else
		{
			PSD_u_K.readFromMatlabFile(InputFolder + "Ku_BC.mat",P.zSlice(P.size_z-1), R.zSlice(R.size_z-1), V.zSlice(V.size_z-1));
		}		
		fclose(file);
	}
	
	// For Testing
			// P.writeToFile("P.plt" , "test"); 
			// R.writeToFile("R.plt" , "test"); 
			// V.writeToFile("V.plt" , "test"); 
			// K.writeToFile("K.plt" , "test"); 	
			// L.writeToFile("L.plt" , "test"); 
			// DLL.writeToFile("DLL.plt" , "test"); 
			// DVV.writeToFile("DVV.plt" , "test"); 
			// DKK.writeToFile("DKK.plt" , "test"); 
			// DVK.writeToFile("DVK.plt" , "test"); 
			// VP.writeToFile("VP.plt" , "test"); 
			// VL.writeToFile("VL.plt" , "test"); 
			// PSD.writeToFile("PSD.plt" , "test"); 
			// G_local.writeToFile("G_local.plt" , "test"); 
			// G_radial.writeToFile("G_radial.plt" , "test"); 
			// Sources.writeToFile("Sources.plt" , "test"); 
			// Losses.writeToFile("Losses.plt" , "test"); 
			// PSD_l_R.writeToFile("PSD_l_R.plt" , "test"); 
			// PSD_u_R.writeToFile("PSD_u_R.plt" , "test"); 
			// PSD_l_L.writeToFile("PSD_l_L.plt" , "test"); 
			// PSD_u_L.writeToFile("PSD_u_L.plt" , "test"); 
			// PSD_l_V.writeToFile("PSD_l_V.plt" , "test"); 
			// PSD_u_V.writeToFile("PSD_u_V.plt" , "test"); 
			// PSD_l_K.writeToFile("PSD_l_K.plt" , "test"); 
			// PSD_u_K.writeToFile("PSD_u_K.plt" , "test");
	
	
	
	// if (!PSD_u_R.readFromIniFile(InputFolder + "Ru_BC.tab", P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1))) {
	// 	PSD_u_R.readFromFile(InputFolder + "Ru_BC.plt",     P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1));
	// }

	// if (!PSD_l_L.readFromIniFile(InputFolder + "Ll_BC.tab", P.xSlice(0), V.xSlice(0), K.xSlice(0))) {
	// 	PSD_l_L.readFromFile(InputFolder + "Ll_BC.plt",     P.xSlice(0), V.xSlice(0), K.xSlice(0));
	// }
	// if (!PSD_u_L.readFromIniFile(InputFolder + "Lu_BC.tab", P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1))) {
	// 	PSD_u_L.readFromFile(InputFolder + "Lu_BC.plt",     P.xSlice(P.size_x-1), V.xSlice(V.size_x-1), K.xSlice(K.size_x-1));
	// }

	// if (!PSD_l_V.readFromIniFile(InputFolder + "Vl_BC.tab", P.ySlice(0), R.ySlice(0), K.ySlice(0))) {
	// 	PSD_l_V.readFromFile(InputFolder + "Vl_BC.plt",     P.ySlice(0), R.ySlice(0), K.ySlice(0));
	// }
	// if (!PSD_u_V.readFromIniFile(InputFolder + "Vu_BC.tab", P.ySlice(P.size_y-1), R.ySlice(R.size_y-1), K.ySlice(K.size_y-1))) {
	// 	PSD_u_V.readFromFile(InputFolder + "Vu_BC.plt",     P.ySlice(P.size_y-1), R.ySlice(R.size_y-1), K.ySlice(K.size_y-1));
	// }

	// if (!PSD_l_K.readFromIniFile(InputFolder + "Kl_BC.tab", P.zSlice(0), R.zSlice(0), V.zSlice(0))) {
	// 	PSD_l_K.readFromFile(InputFolder + "Kl_BC.plt",     P.zSlice(0), R.zSlice(0), V.zSlice(0));
	// }
	// if (!PSD_u_K.readFromIniFile(InputFolder + "Ku_BC.tab", P.zSlice(P.size_z-1), R.zSlice(R.size_z-1), V.zSlice(V.size_z-1))) {
	// 	PSD_u_K.readFromFile(InputFolder + "Ku_BC.plt",     P.zSlice(P.size_z-1), R.zSlice(R.size_z-1), V.zSlice(V.size_z-1));
	// }
	
	// Log the boundary types
	Logger::message << "Ll_BC = " << Ll_BC_type << "; Lu_BC = " << Lu_BC_type << ";" << endl;
	Logger::message << "Vl_BC = " << Vl_BC_type << "; Vu_BC = " << Vu_BC_type << ";" << endl;
	Logger::message << "Kl_BC = " << Kl_BC_type << "; Ku_BC = " << Ku_BC_type << ";" << endl;

	return true;
}
