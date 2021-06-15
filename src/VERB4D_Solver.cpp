/**
 *
 *@mainpage Documentation for the VERB4D_SOLVER
 * Includes classes, namespaces, and files
 *
 * For use in conjunction with Matlab in order to solve for Phase Space Density, Diffusion, Convection etc.
 *
 * Diffusion + convection code, VERB4D.
 *
 * The code solves the Fokker-Planck equation with convection terms:
 * \f[
 * \frac{\partial f}{\partial t} =
 * v_{\Phi} \frac{\partial f}{\partial \Pi} + v_{R} \frac{\partial f}{\partial R} +
 * \frac{1}{G} \frac{\partial}{\partial L} G D_{LL} \frac{\partial f}{\partial L} +
 * \frac{1}{G} \frac{\partial}{\partial V} G (D_{VV} \frac{\partial f}{\partial V} + D_{VK} \frac{\partial f}{\partial K}) +
 * \frac{1}{G} \frac{\partial}{\partial K} G (D_{KV} \frac{\partial f}{\partial V} + D_{KK} \frac{\partial f}{\partial K}) +
 * Losses * f + Sources
 * \f]
 * where f is a Phase Space Density (PSD); Φ is MLT angle, R is radial distance from Earth center, \f$ V≡μ∙(K+0.5)^2 \f$ is a combination of adiabatic
 * invariants μ and \f$ K≡ {J}{\sqrt{8μ}} \f$; L is the third adiabatic invariant; D are bounce-averaged diffusion coefficients; \f$ v_Φ \f$, \f$ v_R \f$ are drift velocities; G is the
 * Jacobian of the transformation from (μ,J,L) to (V,K,L)
 *
 *
 * Books, required to understand the code:
 * - 	Stroustrup C++
 * - 	Schulz and Lanzerotti, 1974
 *
 * VARIABLES for understanding code
 * - 	PSD - Phase Space Density
 * - 	P - Phi: magnetic local time
 * - 	R - Radial distance actual
 * - 	V - Energy
 * - 	K - Pitch andle
 * - 	L - Radial distance distorted
 * - 	Dxx - Diffusion - subscripts denote with respect to what derivative
 *
 * Convection is solved using Convection_1D_ULTIMATE_QUICKEST6.h
 * which implements
 * Leonard BP (1988) Universal Limiter for transient interpolation modeling of the advective transport equations:
 * the ULTIMATE conservative difference scheme, NASA technical Memorandum 100916 ICOMP-88-11
 * http://www.hadian.ir/teaching/CompHydr/3.pdf
 *
 * Diffusion is solved using Diffusion_1D() and Diffusion_2D() for the standard cases.
 *
 * Three other methods are also used for finding diffusion in 2D:
 * these can be found in Diffusion_ADI1.h, Diffusion_ADI2.h, and Diffusion_ADI3.h
 *
 * The parameters.ini file created in matlab will specify the inversion method.
 * If "Lapack" is specified then Diffusion_2D() will be used.
 * Otherwise 1 of the 3 ADI methods will be used to inversion.
 * The current chosen method of inversion is Diffusion_2D() if "Lapack" is specified. It can be changed in VERB4D_SOLVER.cpp
 *
 * From current testing all of the ADI methods produce bad results - namely negative PSD values.
 * All examples should thus be run with Lapack until the ADI methods are fixed.
 * While using Lapack, be sure that you use the thread-safe version of the Lapack library.
 *
 * Four options were introduced to run enable/disable particular split-operator steps: run_remapping, run_convection, run_radial_diffusion, run_local_diffusion.
 * The values can be specified in the paramters.ini file. They are "true" by default.
 *
 * The main solver function can be found in VERB4D_Solver.cpp
 * It consists of a series of PSD calculations done for every time step that is specified.
 *
 * The matrix classes listed above are used to store the data in 1,2,3 or 4 dimension. One example is a 4D matrix containing P,R,V,K
 *
 * There is also the standardized MatrixND for any size matrices which is used in conjunction with UpdatableMatrix which can be updated with new data.
 *
 * Most of the computation such as numerically approximating derivatives are done with MatrixSolver.h using CalculationMatrix.
 * These calculation matrices are made up of DiagMatrix which is a mapping of an int (which diagonal number) to a 1d matrix of values (for that diagonal).
 *
 * There is also a Parameters class which holds paramters and their value as defined in the file they came from.
 * These parameters get the values specified in the parameters.ini text file and set variables to determine inversion method, number of threads, using .mat files, etc.
 *
 * There are three options for reading input and writing output: ascii (.plt), binary (.plt) or matlab (.mat) files.
 * To guarantee backward compatibility of the code, parameter use_matlab is left in the code. If use_matlab == "true", io_method will be overwritten with "matlab".
*
*/

/**
 * \file VERB4D_Solver.cpp
 * \brief Holds main function that solves all of the matrices using the help of the other classes/functions
 *
 * Diffusion + convection code, VERB4D.
 *
 * The code solves the Fokker-Planck equation with convection terms:
 * \f[
 * \frac{\partial f}{\partial t} =
 * v_{\Phi} \frac{\partial f}{\partial \Pi} + v_{R} \frac{\partial f}{\partial R} +
 * \frac{1}{G} \frac{\partial}{\partial L} G D_{LL} \frac{\partial f}{\partial L} +
 * \frac{1}{G} \frac{\partial}{\partial V} G (D_{VV} \frac{\partial f}{\partial V} + D_{VK} \frac{\partial f}{\partial K}) +
 * \frac{1}{G} \frac{\partial}{\partial K} G (D_{KV} \frac{\partial f}{\partial V} + D_{KK} \frac{\partial f}{\partial K}) +
 * Losses * f + Sources
 * \f]
 * where f is a Phase Space Density (PSD); Φ is MLT angle, R is radial distance from Earth center, \f$ V≡μ∙(K+0.5)^2 \f$ is a combination of adiabatic
 * invariants μ and \f$ K≡ {J}{\sqrt{8μ}} \f$; L is the third adiabatic invariant; D are bounce-averaged diffusion coefficients; \f$ v_Φ \f$, \f$ v_R \f$ are drift velocities; G is the
 * Jacobian of the transformation from (μ,J,L) to (V,K,L)
 *
 * The solver first reads in data using ReadInitialData() then procedes to calculate convection and diffusion based on the data it received.
 * Convection_2D() is calculated followed by Diffusion_1D() and Diffusion_2D(). Lastly, the solver calculates Diffusion in 2 dimensions using one of the 3 defined diffusion_ADI functions.
 * The current solver is set-up to use Diffusion_2D_ADI3() currently.
 */

/*
 * Books, required to understand the code:
 * Schulz and Lanzerotti, 1974
 * Stroustrup C++
 *
 * The main function in the standard case calculates PSD by doing radial and then local diffusion at every timestep.
 * The code is set up to have the PSD output at a time step be the initial PSD for the next time step.
 * The PSD output files are usually analyzed/plotted using Matlab function plot_results which is found in every VERB4D example
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <time.h>
#include <math.h>

#include "Logger.h"
#include "Matrix.h"
#include "Diffusion_1D.h"
#include "Diffusion_2D.h"   // Different solution methods for 2D diffusion
#include "Diffusion_ADI1.h" // Straightforward ADI
#include "Diffusion_ADI2.h" // Fuliang Xiao's ADI
#include "Diffusion_ADI3.h" // Jihye Shin and Sungsoo S. Kim (2008)'s ADI - most stable
#include "Convection_2D.h"
#include "ReadInitialData.h"
#include "UpdatableMatrix.h"
#include "Interpolation.h"
#include <omp.h>

#ifdef DATA_ASSIMILATION
#include "DataAssimilation.h"
#endif 

using namespace std;

#ifdef DATA_ASSIMILATION
namespace da = data_assimilation;
#endif 

//#define DEBUG_MODE

#if defined(_WIN32) || defined(_WIN64)
#define strncasecmp _strnicmp
#define strcasecmp  _stricmp
#endif

// Everything below these values will be considered to be zero for computation purpose and will not be calculated
// XXX: It should be a parameter
#define min_PSD 1e-10
#define min_V   1e-10
#define min_Dxx 1e-10

int main(int argc, char* argv[]) {
    Logger::createInstance();
    Logger::message << "Compilation time: " << __DATE__ << ", " << __TIME__ << endl;

    // Variables
    // Grid, 4D:
    // P - local time
    // R - radial distance
    // V - invariant = mu * (K + const) ^ 2, where mu is the first adiabatic invariant
    // K - invariant
    Matrix4D<double> R, P, V, K;

    // Grid for radial diffusion
    // L - L-star, 4D
    UpdatableMatrix<Matrix4D<double> > L;

    // copy of L, for interpolation of PSD to a new L
    Matrix4D<double> L_copy;

    // PSD - phase space density, 4D (P, R, M, I)
    Matrix4D<double> PSD;

    // Diffusion coefficient, 4D
    UpdatableListMatrix<Matrix4D<double> > DLL, DVV, DKK, DVK;

    // Convection (advection) velocities, 4D arrays
    UpdatableMatrix<Matrix4D<double> > VR, VP;

    // Additional sources and losses
    UpdatableListMatrix<Matrix4D<double> > Sources, Losses, Losses_conv;

    // Jacobians, 4D, everything is 4D, it makes matrix operators convenient
    UpdatableMatrix < Matrix4D<double> >  G_local;
    UpdatableMatrix < Matrix4D<double> >  G_radial; // L-star is different for different P, I, K

    // Boundary conditions
    // Pl_BC means P-direction lower boundary (i.e. first index on P-direction)
    // Pu_BC means P-direction upper boundary (i.e. last index on P-direction)
    UpdatableMatrix<Matrix3D<double> > Pl_BC, Pu_BC;
    UpdatableMatrix<Matrix3D<double> > Rl_BC, Ru_BC;
    UpdatableMatrix<Matrix3D<double> > Vl_BC, Vu_BC;
    UpdatableMatrix<Matrix3D<double> > Kl_BC, Ku_BC;
    UpdatableMatrix<Matrix3D<double> > Ll_BC, Lu_BC;
    // Type of boundary condition either BCT_CONSTANT_VALUE or BCT_CONSTANT_DERIVATIVE
    string Pl_BC_type, Pu_BC_type;
    string Rl_BC_type, Ru_BC_type;
    string Vl_BC_type, Vu_BC_type;
    string Kl_BC_type, Ku_BC_type;
    string Ll_BC_type, Lu_BC_type;

    // Parameters input
    int P_size, R_size, V_size, K_size, L_size;

    // time step
    int output_step;
    long int it_total, it_first = 1;

    double dt          = 1.0 / 24;  // 1-hour time interval
    double time_output = 1.0 / 24;  // 1-hour time interval
    double time_total  = 4;         // in days
    double time_first  = 0;
    int max_threads    = omp_get_num_threads();

    // The inversion method can be Lapack or ADI
    // Lapack cannot be used with multiple threads (must be updated to scalapack)
    // Instead ADI should be used for parallelization - It should be changed in the matlab Conv_Dif.m file if using more than 1 thread
    string inputFolder      = "./VERB4D_input/";
    string outputFolder     = "./VERB4D_output/";
    string inversion_method = "Lapack";
    string include_boundary = "true";
    string Vl_BC_from_convection = "false";
    string Vu_BC_from_convection = "false";
    string io_method        = "ascii";
    string run_remapping = "true";
    string run_convection = "true";
    string run_radial_diffusion = "true";
    string run_local_diffusion = "true";
    string positive_PSD = "false";

    bool initialLoad = false; // Check the load of the initial files

    // Read all the inputs - store them into variables
    // These inputs come from the matlab files that are generated when running Conv_Dif.m examples
    initialLoad = ReadInitialData(inputFolder, outputFolder, argc, argv, time_total, dt, time_output, time_first, it_first, max_threads,
            inversion_method, include_boundary,  Vl_BC_from_convection, Vu_BC_from_convection, io_method, run_remapping, run_convection,
            run_radial_diffusion, run_local_diffusion, positive_PSD,
            PSD,P, R, V, K, L,
            P_size, R_size, V_size, K_size, L_size, Pl_BC, Pu_BC, Rl_BC, Ru_BC,
            Vl_BC, Vu_BC, Kl_BC, Ku_BC, Ll_BC, Lu_BC, Pl_BC_type, Pu_BC_type, Rl_BC_type, Ru_BC_type, Vl_BC_type,
            Vu_BC_type, Kl_BC_type, Ku_BC_type, Ll_BC_type, Lu_BC_type, DLL, DVV, DKK, DVK, VP, VR, G_local, G_radial,
            Sources, Losses, Losses_conv);

    // Check that all nesesarry files were loaded
    if (!initialLoad) {
        Logger::error << "Error: ReadInitialData return false. Check the initial files." << endl;
        exit(EXIT_FAILURE);
    }

    #ifdef DATA_ASSIMILATION
    da::DataAssimilationManagerConvection daManagerConvection {
        "parameters_da.ini", time_first, time_first + time_total, V.wxSlice(0,0), K.wxSlice(0,0), P_size, R_size
    };
    #endif

    // Copy L-star so we can later interpolate PSD to a new L-star,
    // to account for adiabatic transport if L-star changes
    L_copy = L;

    // minimum step is 1 hour
    it_total = it_first + round(double(time_total) / dt); // total number of hours given it_first offset
    output_step = round(double(time_output) / dt);
    if (output_step < 1) {
        output_step = 1;
    }

    // logs the timestep and output step information
    Logger::message << "Total time " << time_total << ". Time step " << dt << " (" << it_total << " steps)." << endl;
    Logger::message << "Output each " << output_step << " step. " << endl;

    PSD.writeToFile(outputFolder + "PSD0.plt", P, R, V, K);

    // Output zero step - writing PSD_0 file
    ostringstream PSD_filename, time_string;

    time_string.precision(5);
    time_string.setf(ios::fixed);

    PSD_filename << outputFolder << "PSD_" << setw(5) << setfill('0') << 0;
    Logger::message << "Writing results: " << PSD_filename.str() << endl;
    time_string.str("");
    time_string << time_first;
    PSD.writeToAnyFile(PSD_filename.str(), io_method, time_string.str());

    // When to apply loss term:
    // It's better to apply it during pitch-angle diffusion, unless we don't have any pitch-angle diffusion
    // In that case we apply it during radial diffusion
    int radial_losses, local_losses;
    if (V_size >= 3 && K_size >= 3) {
        radial_losses = 0;
        local_losses = 1;
    } else {
        radial_losses = 1;
        local_losses = 0;
    }

    // If local diffusion is not run, local losses should be applied at the radial diffusion step
    if (run_local_diffusion == "false"){
        local_losses = 0;
        radial_losses = 1;
    }

    // For recording length of time for all calculations to complete
    std::clock_t start_time, tic, tok;
    start_time = clock();
    double wall_timer = omp_get_wtime();

    // Setting up parellization
    omp_set_num_threads(max_threads);
#pragma omp parallel
    {
#pragma omp master
        Logger::message << "Number of threads: " << omp_get_num_threads() << endl;
    }

    // Check time-step for ADI method - the stable time step is completely empirical (i.e. made-up)
    // This isn't really working I think.
    // Only throws error if time step is too large - else has no effect on remaining calculations
    if (inversion_method == "ADI") {
        // This check should be included for other inversion_methods
        if (sqrt(1.0 / V.size_y) <= dt) {
            Logger::error << "Calculating with ADI, time step " << dt << " is too large." << endl;
            exit(EXIT_FAILURE);
        }
        Logger::message << "Calculating with " << "Diffusion_2D_ADI3" << endl;
    } else {
        Logger::message << inversion_method << endl;
    }

    // Indexers to keep track of P,R,V,K,L
    int iP, iR, iV, iK, iL;
    double time;

    // variables to show the progress of calculation
    int progress_count, progress_total;

    // Main loop
    // Start time
    for (long int it = it_first; it < it_total; it++) {
        // update time by dt every iteration
        time = time_first + it * dt;
        Logger::message << endl << std::setprecision(15) << "Time[" << it << "/" << it_total << "]: " << time << " (days)" << endl;

        // Update boundary conditions and diffusion coefficients

        // Update magnetic field (update R)
        if (run_remapping == "true") {
            if (L.update(time, P, R, V, K)) {
                // If we've updated L, we most likely need to update Jacobians
                // XXX: Do we need to update Jacobians if we didn't update L?
                G_local.update(time, P, R, V, K);
                G_radial.update(time, P, R, V, K);

                // If L was updated - interpolate PSD to new L
                progress_count = 0;
                progress_total = P_size * V_size * K_size; // total size of solution matrix
                Logger::message << "Interpolation to new L (adiabatic transport): ";
                cout << "           ";

                Matrix1D<double> old_L_1d(L_size), PSD_L(L_size), new_L_1d(L_size);
                // Aparently it's not thread-safe
                //#pragma omp parallel for private(iP, iR, iV, iK, iL, PSD_L) shared(progress_total, progress_count) schedule(dynamic,1) collapse(3)
                for (iP = 0; iP < P_size; iP++) {
                    for (iV = 0; iV < V_size; iV++) {
                        for (iK = 0; iK < K_size; iK++) {
                            // show progress % if 0 threads
                            if (omp_get_thread_num() == 0) {
                                cout << "\b\b\b\b\b\b\b\b\b" << setw(8)
                                    << (int) ((double) progress_count / progress_total * 100) << "\%" << flush;
                            } else {
                                cout << "thread" << omp_get_thread_num();
                            }

                            // 1d slice to get L from matrix4d (P,L,V,K)
                            new_L_1d = L.wyzSlice(iP, iV, iK);
                            old_L_1d = L_copy.wyzSlice(iP, iV, iK);
                            PSD_L    = PSD.wyzSlice(iP, iV, iK);

                            // 1d interpolationCubic1D
                            PSD_L = Cubic1D(old_L_1d, PSD_L, new_L_1d);

                            // copy results back into PSD adding the 1d list PSD_L for all values of iP,iV,iK
                            for (iL = 0; iL < L_size; iL++) {
                                PSD[iP][iL][iV][iK] = PSD_L[iL];
                            }

                            // Progress output - will update P_size * V_size * K_size times
                            //#pragma omp critical
                            progress_count += 1;
                        }
                    }
                }
                cout << "\b\b\b\b\b\b\b\b\b" << setw(8) << (int) ((double) progress_count / progress_total * 100) << "\%" << endl;

                // Copy the new L into L_copy for future interpolations
                L_copy = L;
            }
        }

        if (positive_PSD == "true") {
            for (iP = 0; iP < P_size; iP++) {
                for (iR = 0; iR < R_size; iR++) {
                    for (iV = 0; iV < V_size; iV++) {
                        for (iK = 0; iK < K_size; iK++) {
                            if (PSD[iP][iR][iV][iK] < 1e-21) PSD[iP][iR][iV][iK] = 1e-21;
                        }
                    }
                }
            }
        }

        // Update convection velocities VP and VR and log the maximum absolute values
        if ((3 < P_size || 3 < R_size) && (run_convection == "true")) {
            VP.update(time, P, R, V, K);
            Logger::message << "max(VP) = " << VP.maxabs() << endl;
            VR.update(time, P, R, V, K);
            Logger::message << "max(VR) = " << VR.maxabs() << endl;
        }

        // Diffusion coefficients
        if ((3 < L_size) && (run_radial_diffusion == "true")) {
            DLL.update(time, P, L, V, K);
            Logger::message << "max(DLL) = " << DLL.maxabs() << endl;
        }

        if ((3 < V_size && 3 < K_size) && (run_local_diffusion == "true")) {
            DVV.update(time, P, R, V, K);
            Logger::message << "max(DVV) = " << DVV.maxabs() << endl;
            DVK.update(time, P, R, V, K);
            Logger::message << "max(DVK) = " << DVK.maxabs() << endl;
            DKK.update(time, P, R, V, K);
            Logger::message << "max(DKK) = " << DKK.maxabs() << endl;
        }

        // Sources and losses
        Sources.update(time, P, R, V, K);
        Losses.update(time, P, R, V, K);
        Losses_conv.update(time, P, R, V, K);

        // Boundary conditions
        // By default - it's constant PSD
        if (R_size > 3) {
            // Put PSD-slice into the boundary condition. It will be updated from the ini-file, if there is something in the ini-file.
            //R_LBC.original_arr =(PSD.xSlice(0));
            Rl_BC.update(time, P.xSlice(0), V.xSlice(0), K.xSlice(0));
            //R_UBC.original_arr =(PSD.xSlice(PSD.size_x - 1));
            Ru_BC.update(time, P.xSlice(R.size_x - 1), V.xSlice(R.size_x - 1), K.xSlice(R.size_x - 1));
        }
        if (L_size > 3) {
            //L_LBC.original_arr =(PSD.xSlice(0));
            Ll_BC.update(time, P.xSlice(0), V.xSlice(0), K.xSlice(0));
            //L_UBC.original_arr =(PSD.xSlice(PSD.size_x - 1));
            Lu_BC.update(time, P.xSlice(R.size_x - 1), V.xSlice(R.size_x - 1), K.xSlice(R.size_x - 1));

        }
        if (V_size > 3) {
            //V_LBC.original_arr =(PSD.ySlice(0));
            Vl_BC.update(time, P.ySlice(0), R.ySlice(0), K.ySlice(0));
            //V_UBC.original_arr =(PSD.ySlice(PSD.size_y - 1));
            Vu_BC.update(time, P.ySlice(V.size_y - 1), R.ySlice(V.size_y - 1), K.ySlice(V.size_y - 1));
        }
        if (K_size > 3) {
            //K_LBC.original_arr =(PSD.zSlice(0));
            Kl_BC.update(time, P.zSlice(0), R.zSlice(0), V.ySlice(0));
            //K_UBC.original_arr =(PSD.zSlice(PSD.size_z - 1));
            Ku_BC.update(time, P.zSlice(K.size_z - 1), R.zSlice(K.size_z - 1), V.zSlice(K.size_z - 1));
        }

        // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Check that DQ2Q3^2 < DQ2Q2 * DQ3Q3
        // The computation will be unstable otherwise

        // checking if DVV^2 - min(DVK^2) is greater than 0, if not log and exit
        double check = (DVV.times(DKK) - DVK.times(DVK)).min();
        if (check < 0) {
            Logger::error << "Fatal error: (DVV*DKK - DVK*DVK) = " << check << " < 0." << endl;
            Logger::error << "The computation would be unstable." << endl;
            exit(EXIT_FAILURE);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Convection, for each V and K (therefore for each mu and J)
        if ((3 < P_size && 3 < R_size) && (run_convection == "true")){
            progress_count = 0;
            progress_total = V_size * K_size;
            Logger::message << "Convection:" << endl;;
            cout << "           ";

            Matrix2D<double> PSD_PR(P_size, R_size);

            // If you have OPENMP 3.0 with collapse functionality you can use this. (Visual Studios does not)
#if defined _OPENMP && _OPENMP >= 200711
#pragma omp parallel for private(iP, iR, iV, iK, iL, PSD_PR) shared(progress_total, progress_count) schedule(dynamic,1) collapse(2)

            // Looping it backward allows to speed-up the multithread simulation
            for (iV = V_size-1; iV >= 0; iV--) {
                // due to the highest energies being the slowest to calculate - calculating highest energy first
                for (iK = 0; iK < K_size; iK++) {
                    // Output current progress percentage when number of threads = 0
                    if (omp_get_thread_num() == 0) {
                        cout << "\b\b\b\b\b\b\b\b\b" << setw(8)
                             << (int) ((double) progress_count / progress_total * 100) << "\%" << flush;
                    }

                    // 2d slice of PSD to pass to the Convection calculation function
                    PSD_PR = PSD.yzSlice(iV, iK);

                    // Inputs to the Convection_2D, Sources and Losses set to 0
                    Convection_2D(PSD_PR, P.yzSlice(iV, iK), R.yzSlice(iV, iK), P_size, R_size, Pl_BC.yzSlice(iV, iK),
                            Pu_BC.yzSlice(iV, iK), // R, I, K
                            Rl_BC.yzSlice(iV, iK), Ru_BC.yzSlice(iV, iK), // P, I, K
                            Pl_BC_type, Pu_BC_type, Rl_BC_type, Ru_BC_type, VP.yzSlice(iV, iK), VR.yzSlice(iV, iK),
                            Sources.yzSlice(iV, iK) * 0, Losses_conv.yzSlice(iV, iK), dt, min_PSD, min_V);

                    // copy results back into PSD adding the 2d list PSD_PR for all values of iV,iK
                    for (iP = 0; iP < P_size; iP++) {
                        for (iR = 0; iR < R_size; iR++) {
                            PSD[iP][iR][iV][iK] = PSD_PR[iP][iR];
                        }
                    }

#pragma omp critical
                    progress_count += 1;
                }
            }

#else // If you dont have the latest OPENMP will manually collapse indeces for you
#pragma omp parallel for private(iP, iR, iV, iK, iL, PSD_PR) shared(progress_total, progress_count) schedule(dynamic,1)

            // Looping it backward allows to speed-up the multithread simulation
            for (int index = (V_size * K_size) - 1; index >= 0; index--) {
                // due to the highest energies being the slowest to calculate - calculating highest energy first

                // Manually collapsing into a single for loop in order to optimize
                iV = index % V_size;
                iK = (K_size - 1) - (index / V_size);

                // Output current progress percentage when number of threads = 0
                if (omp_get_thread_num() == 0) {
                    cout << "\b\b\b\b\b\b\b\b\b" << setw(8)
                        << (int) ((double) progress_count / progress_total * 100) << "\%" << flush;
                }

                // 2d slice of PSD to pass to the Convection calculation function
                PSD_PR = PSD.yzSlice(iV, iK);

                // Inputs to the Convection_2D, Sources and Losses set to 0
                Convection_2D(PSD_PR, P.yzSlice(iV, iK), R.yzSlice(iV,  iK), P_size, R_size, Pl_BC.yzSlice(iV,  iK),
                        Pu_BC.yzSlice(iV, iK), // R, I, K
                        Rl_BC.yzSlice(iV, iK), Ru_BC.yzSlice(iV,  iK), // P, I, K
                        Pl_BC_type, Pu_BC_type, Rl_BC_type, Ru_BC_type, VP.yzSlice(iV, iK), VR.yzSlice(iV, iK),
                        Sources.yzSlice(iV, iK) * 0, Losses_conv.yzSlice(iV,  iK), dt, min_PSD, min_V);

                // copy results back into PSD adding the 2d list PSD_PR for all values of iV,iK
                for (iP = 0; iP < P_size; iP++) {
                    for (iR = 0; iR < R_size; iR++) {
                        PSD[iP][iR][iV][iK] = PSD_PR[iP][iR];
                    }
                }

#pragma omp critical
                // update progress count every iteration through
                progress_count += 1;
            }
#endif
            // Output final progress (it should be 100%)
            cout << "\b\b\b\b\b\b\b\b\b" << setw(8) << (int) ((double) progress_count / progress_total * 100) << "\%" << endl;
#pragma omp master
            {
            if(Vl_BC_from_convection == "true" && (Vl_BC_type == "BCT_CONSTANT_VALUE")){ //rewrite boundary conditions at lower V
                Vl_BC = PSD.ySlice(0);
                cout << "Vl_BC from convection are used: max(Vl_BC) = " << Vl_BC.max() << endl;
            }
            if(Vu_BC_from_convection == "true" && (Vu_BC_type == "BCT_CONSTANT_VALUE")){ //rewrite boundary conditions at lower V
                Vu_BC = PSD.ySlice(V_size-1);
                cout << "Vu_BC from convection are used: max(Vu_BC) = " << Vu_BC.max() << endl;
            }
            }
        }

        if (positive_PSD == "true") {
            for (iP = 0; iP < P_size; iP++) {
                for (iR = 0; iR < R_size; iR++) {
                    for (iV = 0; iV < V_size; iV++) {
                        for (iK = 0; iK < K_size; iK++) {
                            if (PSD[iP][iR][iV][iK] < 1e-21) PSD[iP][iR][iV][iK] = 1e-21;
                        }
                    }
                }
            }
        }

        // RADIAL DIFFUSION STEP

        // ADDED FOR TESTING
        //  PSD.writeToFile(to_string(int(it / output_step)) +  "PSD_before_radial.plt");

        if ((L_size >= 3) && (run_radial_diffusion == "true")) {
            progress_count = 0;
            progress_total = P_size * V_size * K_size; // total size of solution matrix
            Logger::message << "Radial diffusion:" << endl;;
            cout<< "           ";

            Matrix1D<double> PSD_L(L_size);

#if defined _OPENMP && _OPENMP >= 200711
#pragma omp parallel for private(iP, iR, iV, iK, iL, PSD_L) shared(progress_total, progress_count) schedule(dynamic,1) collapse(3)
            for (iP = 0; iP < P_size; iP++) {
                for (iV = 0; iV < V_size; iV++) {
                    for (iK = 0; iK < K_size; iK++) {
                        // print percentage done
                        if (omp_get_thread_num() == 0){
                            cout << "\b\b\b\b\b\b\b\b\b" << setw(8) << (int) ((double) progress_count / progress_total * 100) << "\%" << flush;
                        }

                        // 1d slice
                        PSD_L = PSD.wyzSlice(iP, iV, iK);

                        // 1d diffusion
                        Diffusion_1D(PSD_L, L.wyzSlice(iP, iV, iK), L_size, Ll_BC[iP][iV][iK], Lu_BC[iP][iV][iK],
                                Ll_BC_type, Lu_BC_type, DLL.wyzSlice(iP, iV, iK), G_radial.wyzSlice(iP, iV, iK),
                                Sources.wyzSlice(iP, iV, iK) * radial_losses,
                                Losses.wyzSlice(iP, iV, iK) * radial_losses, dt);


                        // copy results back
                        for (iL = 0; iL < L_size; iL++) {
                            PSD[iP][iL][iV][iK] = PSD_L[iL];
                        }

#pragma omp critical
                        progress_count += 1;
                    }
                }
            }

            // ADDED FOR TESTING
            //  PSD.writeToFile(to_string(int(it / output_step)) +  "PSD_after_radial.plt");
#else
#pragma omp parallel for private(iP, iV, iK, iL, PSD_L) shared(progress_total, progress_count)
            for (int index = 0; index < P_size * V_size * K_size; index++) {
                iP = index / (V_size * K_size);
                iV = (index / K_size) % V_size;
                iK = index % K_size;

                // print percentage done
                if (omp_get_thread_num() == 0){
                    cout << "\b\b\b\b\b\b\b\b\b" << setw(8) << (int) ((double) progress_count / progress_total * 100) << "\%" << flush;
                }

                // 1d slice
                PSD_L = PSD.wyzSlice(iP, iV, iK);

                // 1d diffusion
                Diffusion_1D(PSD_L, L.wyzSlice(iP, iV, iK), L_size, Ll_BC[iP][iV][iK], Lu_BC[iP][iV][iK],
                        Ll_BC_type, Lu_BC_type, DLL.wyzSlice(iP, iV, iK), G_radial.wyzSlice(iP, iV, iK),
                        Sources.wyzSlice(iP, iV, iK) * radial_losses,
                        Losses.wyzSlice(iP, iV, iK) * radial_losses, dt);

                // copy results back
                for (iL = 0; iL < L_size; iL++) {
                    PSD[iP][iL][iV][iK] = PSD_L[iL];
                }

#pragma omp critical
                // Progress output
                progress_count += 1;
            }
#endif
            cout << "\b\b\b\b\b\b\b\b\b" << setw(8) << (int) ((double) progress_count / progress_total * 100) << "\%" << endl;
//#pragma omp master
//            {
//            if((Vl_BC_from_convection == "true") && (Vl_BC_type == "BCT_CONSTANT_VALUE")) { //rewrite boundary conditions at lower V
//                Vl_BC = PSD.ySlice(1); // NOTE: low-V boundary conditions are taken from the next-to-lower-V slice after radial diffusion.
//                                       // This eliminates unrealistic gradient there, but it is probably NOT A CORRECT IMPLEMENTATION.
//                //cout << "Vl_BC after radial diffusion are used: max(Vl_BC) = " << Vl_BC.max() << endl;
//                Logger::message << "Vl_BC after radial diffusion are used: max(Vl_BC) = " << Vl_BC.max() << endl;
//            }
//            }
        }

        if (positive_PSD == "true") {
            for (iP = 0; iP < P_size; iP++) {
                for (iR = 0; iR < R_size; iR++) {
                    for (iV = 0; iV < V_size; iV++) {
                        for (iK = 0; iK < K_size; iK++){
                            if (PSD[iP][iR][iV][iK] < 1e-21) PSD[iP][iR][iV][iK] = 1e-21;
                        }
                    }
                }
            }
        }

        // ADDED FOR TESTING
        //  PSD.writeToFile(to_string(int(it / output_step)) +  "PSD_after_radial.plt");

        // LOCAL DIFFUSION
        if ((V_size >= 3 && K_size >= 3) && (run_local_diffusion == "true")){
            int number_of_skipped_points = 0;
            progress_count = 0;
            progress_total = P_size * R_size;
            Logger::message << "Local diffusion: " << endl;
            cout << "           ";

            Matrix2D<double> PSD_IK(V_size, K_size);

#if defined _OPENMP && _OPENMP >= 200711
#pragma omp parallel for private(iP, iR, iV, iK, iL, PSD_IK) shared(progress_total, progress_count, number_of_skipped_points) schedule(dynamic,1) collapse(2)
            for (iP = 0; iP < P_size; iP++) {
                for (iR = 0; iR < R_size; iR++) {
                    // If we dont want to include the boundary in the L/R axis
                    // we can exclude those points by setting include_boundary
                    // to false in matlab main file or parameter.ini file
                    if ((include_boundary == "false") && (iR == 0 || iR == R_size-1)) {
                        continue;
                    }

                    if (omp_get_thread_num() == 0) {
                        cout << "\b\b\b\b\b\b\b\b\b" << setw(8)
                            << (int) ((double) progress_count / progress_total * 100) << "\%" << flush;
                    }

#pragma omp critical
                    progress_count += 1;

                    // 2d slice
                    PSD_IK = PSD.wxSlice(iP, iR);

                    // Don't calculate anything, if diffusion coefficient are zero or PSD is zero
                    // 					if (PSD_IK.max() < min_PSD || (DKK.wxSlice(iP, iR).max() < min_Dxx && DVV.wxSlice(iP, iR).max() < min_Dxx)) {
                    // #pragma omp critical // Progress count
                    // 						number_of_skipped_points++;
                    // 						continue;
                    // 					}


                    // 2d diffusion
                    // If parameters.ini specify "Lapack" then Lapack will be used
                    if (inversion_method == "Lapack") {
                        Diffusion_2D(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
                                Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR), // P, R, K
                                Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR), // P, R, I
                                Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
                                DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
                                Sources.wxSlice(iP, iR) * local_losses, Losses.wxSlice(iP, iR) * local_losses, dt);
                    }
                    // Currently setup to calculate 2d Diffusion using Diffusion_2D_ADI3
                    // Can change to Diffusion_2D_ADI1 or Diffusion_2D_ADI2 for different methods of inversion
                    else if (inversion_method == "ADI") {
                        Diffusion_2D_ADI3(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
                                Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR), // P, R, K
                                Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR), // P, R, I
                                Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
                                DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
                                Sources.wxSlice(iP, iR) * local_losses, Losses.wxSlice(iP, iR) * local_losses, dt);
                    }
                    else if (inversion_method == "ADI2") {
                        Diffusion_2D_ADI2(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
                                Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR), // P, R, K
                                Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR), // P, R, I
                                Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
                                DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
                                Sources.wxSlice(iP, iR) * local_losses, Losses.wxSlice(iP, iR) * local_losses, dt);
                    }
                    else if (inversion_method == "ADI1") {
                        Diffusion_2D_ADI1(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
                                Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR), // P, R, K
                                Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR), // P, R, I
                                Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
                                DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
                                Sources.wxSlice(iP, iR) * local_losses, Losses.wxSlice(iP, iR) * local_losses, dt);
                    }
                    else {
                        Logger::error << "Error: Unknown inversion method " << inversion_method << endl;
                    }

                    // copy results back
                    for (iV = 0; iV < V_size; iV++) {
                        for (iK = 0; iK < K_size; iK++) {
                            PSD[iP][iR][iV][iK] = PSD_IK[iV][iK];
                        }
                    }
                }
            }
#else
#pragma omp parallel for private(iP, iR, iV, iK, iL, PSD_IK) shared(progress_total, progress_count, number_of_skipped_points) schedule(dynamic,1)
            for (int index = 0; index < P_size * R_size; index++) {
                iP = index / R_size;
                iR = index % R_size;

                // If we dont want to include the boundary in the L/R axis we
                // can exclude those points by setting include_boundary to
                // false in matlab main file or parameter.ini file
                if ((include_boundary == "false") && (iR == 0 || iR == R_size-1)) {
                    continue;
                }

                if (omp_get_thread_num() == 0) {
                    cout << "\b\b\b\b\b\b\b\b\b" << setw(8)
                        << (int) ((double) progress_count / progress_total * 100) << "\%" << flush;
                }
#pragma omp critical
                progress_count += 1;

                // 2d slice
                PSD_IK = PSD.wxSlice(iP, iR);

                // Don't calculate anything, if diffusion coefficient are zero or PSD is zero
                //					if (PSD_IK.max() < min_PSD || (DKK.wxSlice(iP, iR).max() < min_Dxx && DVV.wxSlice(iP, iR).max() < min_Dxx)) {
                //#pragma omp critical // Progress count
                //						number_of_skipped_points++;
                //						continue;
                //					}

                // 2d diffusion
                // If parameters.ini specify "Lapack" then Lapack will be used
                if (inversion_method == "Lapack") {
                    Diffusion_2D(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
                            Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR), // P, R, K
                            Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR), // P, R, I
                            Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
                            DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
                            Sources.wxSlice(iP, iR) * local_losses, Losses.wxSlice(iP, iR) * local_losses, dt);
                }
                // Currently setup to calculate 2d Diffusion using Diffusion_2D_ADI3
                // Can change to Diffusion_2D_ADI1 or Diffusion_2D_ADI2 for different methods of inversion
                else if(inversion_method == "ADI"){
                    Diffusion_2D_ADI3(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
                            Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR), // P, R, K
                            Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR), // P, R, I
                            Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
                            DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
                            Sources.wxSlice(iP, iR) * local_losses, Losses.wxSlice(iP, iR) * local_losses, dt);
                }
                else if (inversion_method == "ADI2") {
                    Diffusion_2D_ADI2(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
                            Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR), // P, R, K
                            Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR), // P, R, I
                            Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
                            DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
                            Sources.wxSlice(iP, iR) * local_losses, Losses.wxSlice(iP, iR) * local_losses, dt);
                }
                else if (inversion_method == "ADI1") {
                    Diffusion_2D_ADI1(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
                            Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR), // P, R, K
                            Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR), // P, R, I
                            Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
                            DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
                            Sources.wxSlice(iP, iR) * local_losses, Losses.wxSlice(iP, iR) * local_losses, dt);
                }
                else {
                    Logger::error << "Error: Unknown inversion method " << inversion_method << endl;
                }

                // copy results back
                for (iV = 0; iV < V_size; iV++) {
                    for (iK = 0; iK < K_size; iK++) {
                        PSD[iP][iR][iV][iK] = PSD_IK[iV][iK];
                    }
                }
            }
#endif
            cout << "\b\b\b\b\b\b\b\b\b" << setw(8) << (int) ((double) progress_count / progress_total * 100) << "\%" << endl;
            cout << "Number of skipped points: " << (int) ((double) number_of_skipped_points/progress_total * 100) << "\%" << endl;
        }

        if (positive_PSD == "true") {
            for (iP = 0; iP < P_size; iP++) {
                for (iR = 0; iR < R_size; iR++) {
                    for (iV = 0; iV < V_size; iV++) {
                        for (iK = 0; iK < K_size; iK++) {
                            if (PSD[iP][iR][iV][iK] < 1e-21) PSD[iP][iR][iV][iK] = 1e-21;
                        }
                    }
                }
            }
        }

        int number_of_negative_points = 0;
        for (iP = 0; iP < P_size; iP++) {
            for (iR = 0; iR < R_size; iR++) {
                for (iV = 0; iV < V_size; iV++) {
                    for (iK = 0; iK < K_size; iK++) {
                        if (PSD[iP][iR][iV][iK] < 0) number_of_negative_points++;
                    }
                }
            }
        }
        //cout << "Number of negative points: " << number_of_negative_points << endl;
        Logger::message << endl << "Number of negative points: " << number_of_negative_points << " of " << P_size*R_size*V_size*K_size << endl;

#ifdef DATA_ASSIMILATION
        daManagerConvection.assimilate(time, PSD, P, R, VP, VR, Losses_conv, dt);
        if (positive_PSD == "true") {
            for (iP = 0; iP < P_size; iP++) {
                for (iR = 0; iR < R_size; iR++) {
                    for (iV = 0; iV < V_size; iV++) {
                        for (iK = 0; iK < K_size; iK++) {
                            if (PSD[iP][iR][iV][iK] < 1e-21) PSD[iP][iR][iV][iK] = 1e-21;
                        }
                    }
                }
            }
        }
#endif

        // Output the PSD data for each timestep into the output folder
        if ((it % output_step) == 0) {
            time_string.str("");
            time_string << time;

            PSD_filename.str("");
            PSD_filename << outputFolder << "PSD_" << setw(5) << setfill('0') << int(it / output_step);
            Logger::message << endl << "Writing results: " << PSD_filename.str() << endl;
            PSD.writeToAnyFile(PSD_filename.str(), io_method, time_string.str());
        }
    }

    // logger records if everything went correctly
    Logger::message << "Program was terminated correctly." << endl;
    Logger::message << "Wall-clock time: " << (omp_get_wtime() - wall_timer) << " sec; ";
    Logger::message << "CPU time: " << (double)((clock() - start_time) / CLOCKS_PER_SEC) << " sec." << endl;

    Logger::deleteInstance();
    return 0;
}
