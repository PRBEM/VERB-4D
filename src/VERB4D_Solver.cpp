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

#include <omp.h>

#include <chrono>
#include <cmath>
#include <ctime>
#include <future>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>

// enable alternative tokens
#ifdef _MSC_VER
    #include<iso646.h>
#endif

#include "Convection_2D.h"
#include "Convection_3D.h"
#include "Diffusion_1D.h"
#include "Diffusion_2D.h"    // Different solution methods for 2D diffusion
#include "Diffusion_2D_MKL.hpp"    // Intel MKL sparse solve 
#include "Diffusion_ADI1.h"  // Straightforward ADI
#include "Diffusion_ADI2.h"  // Fuliang Xiao's ADI
#include "Diffusion_ADI3.h"  // Jihye Shin and Sungsoo S. Kim (2008)'s ADI - most stable
#include "Interpolation.h"
#include "Logger.h"
#include "Matrix.h"
#include "ReadInitialData.h"
#include "UpdatableMatrix.h"

#ifdef DATA_ASSIMILATION
#include "DataAssimilation.h"
#endif

using namespace std;
#include "BoundaryConditionType.hpp"

#ifdef DATA_ASSIMILATION
namespace da = data_assimilation;
#endif 
#ifdef MKL_FOUND
#include "mkl_sparse_qr.h"
#endif
#include <omp.h>

using namespace std;

// #define DEBUG_MODE

#if defined(_WIN32) || defined(_WIN64)
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

// Everything below these values will be considered to be zero for computation purpose and will not be calculated
// XXX: It should be a parameter
#define min_PSD 1e-10
#define min_V 1e-10
#define min_Dxx 1e-10

void write_PSD_output(const IOMethod &io_method, const std::string &PSD_filename, const double &time, const Matrix4D<double> PSD, const bool &PSD_time_to_lst, const std::string &output_folder)
{
    std::ostringstream time_string;
    time_string.precision(5);
    time_string.setf(std::ios::fixed);
    time_string.str("");
    time_string << time;

    PSD.writeToAnyFile(PSD_filename, io_method, time_string.str());
    if (PSD_time_to_lst)
    {        
        PSD.writeToLstFile(PSD_filename, io_method, time_string.str(), output_folder);
    }

};

int main(int argc, char *argv[])
{
    Logger::createInstance();
    std::string compile_type = "Debug";
    std::string fast_convection = "";
#ifdef RELEASE
    compile_type = "Release";
#endif
#ifdef FAST_CONVECTION
    fast_convection = " (with fast convection)";
#endif
    Logger::message << compile_type << fast_convection << " compiled on: " << __DATE__ << ", " << __TIME__ << std::endl;

    // Variables
    // Grid, 4D:
    // P - local time
    // R - radial distance
    // V - invariant = mu * (K + const) ^ 2, where mu is the first adiabatic invariant
    // K - invariant
    Matrix4D<double> R, P, V, K;

    // Grid for radial diffusion
    // L - L-star, 4D
    UpdatableMatrix<Matrix4D<double>> L;

    // copy of L, for interpolation of PSD to a new L
    Matrix4D<double> L_copy;

    // PSD - phase space density, 4D (P, R, M, I)
    Matrix4D<double> PSD;

    // Diffusion coefficient, 4D
    UpdatableListMatrix<Matrix4D<double>> DLL, DVV, DKK, DVK;

    // Convection (advection) velocities, 4D arrays
    UpdatableMatrix<Matrix4D<double>> VR, VP, VV;

    // Additional sources and losses
    UpdatableListMatrix<Matrix4D<double>> Sources, Losses_conv, Losses_local, Losses_radial;

    // Jacobians, 4D, everything is 4D, it makes matrix operators convenient
    UpdatableMatrix<Matrix4D<double>> G_local;
    UpdatableMatrix<Matrix4D<double>> G_radial; // L-star is different for different P, I, K
    UpdatableMatrix<Matrix4D<double>> G_conv;   // solving the convection in conservative form requires a Jacobian 

    // Boundary conditions
    // Pl_BC means P-direction lower boundary (i.e. first index on P-direction)
    // Pu_BC means P-direction upper boundary (i.e. last index on P-direction)
    UpdatableMatrix<Matrix3D<double>> Pl_BC, Pu_BC;
    UpdatableMatrix<Matrix3D<double>> Rl_BC, Ru_BC;
    UpdatableMatrix<Matrix3D<double>> Vl_BC, Vu_BC;
    UpdatableMatrix<Matrix3D<double>> Kl_BC, Ku_BC;
    UpdatableMatrix<Matrix3D<double>> Ll_BC, Lu_BC;
    // Type of boundary condition either BCT_CONSTANT_VALUE or BCT_CONSTANT_DERIVATIVE
    BoundaryConditionType Pl_BC_type, Pu_BC_type;
    BoundaryConditionType Rl_BC_type, Ru_BC_type;
    BoundaryConditionType Vl_BC_type, Vu_BC_type;
    BoundaryConditionType Kl_BC_type, Ku_BC_type;
    BoundaryConditionType Ll_BC_type, Lu_BC_type;

    // Matrices for plasmasphere simulations
    Matrix4D<double> SaturationDensity;
    Matrix4D<double> SaturationTimescale;

    // Parameters input
    int P_size, R_size, V_size, K_size, L_size;

    // time step
    int output_step;
    long int it_total, it_first = 1;

    double dt = 1.0 / 24;          // 1-hour time interval
    double sub_dt_diffusion = 1.0 / 24;
    double time_output = -1;       // 1-hour time interval
    double time_total = 4;         // in days
    double time_first = 0;
    int max_threads = omp_get_num_threads();
    bool minimal_output = false;

    // The inversion method can be Lapack or ADI
    // Lapack cannot be used with multiple threads (must be updated to scalapack)
    // Instead ADI should be used for parallelization - It should be changed in the matlab Conv_Dif.m file if using more than 1 thread
    std::string inputFolder = "./VERB4D_input/";
    std::string outputFolder = "./VERB4D_output/";
    InversionMethod inversion_method = InversionMethod::Lapack;
    IOMethod io_method = IOMethod::ASCII;
    IOMethod PSD0_io_method = IOMethod::ASCII;
    DensitySaturation density_saturation = DensitySaturation::Off;
    bool include_boundary = true;
    bool Vl_BC_from_convection = false;
    bool Vu_BC_from_convection = false;
    bool run_remapping = true;
    bool run_convection = true;
    bool run_coulomb_collision = false;
    bool run_radial_diffusion = true;
    bool run_local_diffusion = true;
    bool positive_PSD = false;
    bool PSD_time_to_lst = false;

    bool initialLoad = false; // Check the load of the initial files

    // Read all the inputs - store them into variables
    // These inputs come from the matlab files that are generated when running Conv_Dif.m examples
    initialLoad = ReadInitialData(
        minimal_output, inputFolder, outputFolder, argc, argv, time_total, dt, sub_dt_diffusion, time_output, time_first, it_first, max_threads,
        inversion_method, io_method, PSD0_io_method, density_saturation, include_boundary, Vl_BC_from_convection, Vu_BC_from_convection, run_remapping, run_convection,
        run_radial_diffusion, run_local_diffusion, run_coulomb_collision, positive_PSD, PSD_time_to_lst,
        PSD, P, R, V, K, L,
        P_size, R_size, V_size, K_size, L_size, Pl_BC, Pu_BC, Rl_BC, Ru_BC,
        Vl_BC, Vu_BC, Kl_BC, Ku_BC, Ll_BC, Lu_BC, Pl_BC_type, Pu_BC_type, Rl_BC_type, Ru_BC_type, Vl_BC_type,
        Vu_BC_type, Kl_BC_type, Ku_BC_type, Ll_BC_type, Lu_BC_type, DLL, DVV, DKK, DVK, VP, VR, VV, G_local, G_radial, G_conv,
        Sources, Losses_local, Losses_radial, Losses_conv, SaturationDensity, SaturationTimescale);

    // Check that all necessary files were loaded
    if (!initialLoad)
    {
        Logger::error << "Error: ReadInitialData return false. Check the initial files." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (minimal_output) {
    	Logger::setDebugLevel(Logger::DebugLevel::DEBUG_LEVEL_MESSAGE);
    }

#ifdef DATA_ASSIMILATION
    da::DataAssimilationManagerConvection daManagerConvection{
        "parameters_da.ini", time_first, time_first + time_total, P.yzSlice(0, 0), R.yzSlice(0, 0), V.wxSlice(0, 0), K.wxSlice(0, 0), outputFolder};
#endif

#ifdef SAVE_PSD_LOST_CONV
    Matrix4D<double> PSD_lost_conv;
    PSD_lost_conv.AllocateMemory(P_size, R_size, V_size, K_size);
#endif

#ifdef LU_CACHING
    CalculationMatrix dummy(V_size, K_size, 1, 1);
    long m_size = dummy[0].size_q1;
	long kl = -dummy.begin()->first; // first diagonal
	long ku = dummy.rbegin()->first;
    std::vector<std::vector<double*>> lu_cache(P_size);
    for(auto& r_slice : lu_cache)
    {
        r_slice = std::vector<double*>(R_size);
        for(auto& pr_slice : r_slice)
        {
            pr_slice = new double[(kl+ku+kl+1)*m_size];
        }
    }
    std::vector<std::vector<long*>> index_cache(P_size);
    for(auto& r_slice : index_cache)
    {
        r_slice = std::vector<long*>(R_size);
        for(auto& pr_slice : r_slice)
        {
            pr_slice = new long[m_size];
        }
    }

    bool recompute_lu_decomp = true;
#endif

    // Copy L-star so we can later interpolate PSD to a new L-star,
    // to account for adiabatic transport if L-star changes
    L_copy = L;

    // minimum step is 1 hour
    it_total = (long)(it_first + round(double(time_total) / dt)); // total number of hours given it_first offset
    output_step = (int)round(double(time_output) / dt);
    if (output_step < 1)
    {
        output_step = 1;
    }

    // logs the timestep and output step information
	Logger::message << std::endl;
	Logger::writeSeparator();
	Logger::message << std::setw(57) << "Time step info" << std::endl;
	Logger::writeSeparator();

    Logger::message << "Total time " << time_total << ". Time step " << dt << " (" << it_total << " steps)." << std::endl;
    Logger::message << "Output each " << output_step << " step. " << std::endl;

    if (PSD0_io_method == IOMethod::ASCII)
        PSD.writeToFile(outputFolder + "PSD0.plt", P, R, V, K);
    else if (PSD0_io_method == IOMethod::Binary)
        PSD.writeToBinaryFile(outputFolder + "PSD0.pltb");

    // Output zero step - writing PSD_0 file
    std::future<void> output_writer;
    std::ostringstream PSD_filename, time_string;

    // time_string.precision(5);
    // time_string.setf(ios::fixed);

    PSD_filename << outputFolder << "PSD_" << std::setw(5) << std::setfill('0') << 0;
    Logger::debug << "Writing results: " << PSD_filename.str() << std::endl;
    output_writer = std::async(std::launch::async, write_PSD_output, io_method, PSD_filename.str(), time_first, PSD, PSD_time_to_lst, outputFolder);

    // When to apply loss term:
    // It's better to apply it during pitch-angle diffusion, unless we don't have any pitch-angle diffusion
    // In that case we apply it during radial diffusion
    if (not Losses_radial.initialized) { // Losses_radial.tab was not found -> using Losses_local either during radial or local diffusion
        Losses_radial = Losses_local;
        
        if (run_local_diffusion and V_size >= 3 and K_size >= 3) {
            Losses_radial = 0;
            Losses_radial.clearMatricesList();
        } else {
            Losses_local = 0;
            Losses_local.clearMatricesList();
        }
    }

    // If local diffusion is not run, sources should be applied at the radial diffusion step
    int radial_sources = 0;
    int local_sources = 1;
    if (not run_local_diffusion or V_size < 3 or K_size < 3) {
        // Apply source during radial diffusion step
        radial_sources = 1;
        local_sources = 0;
    }

    // For recording length of time for all calculations to complete
    std::clock_t start_time;
    start_time = clock();
    double wall_timer = omp_get_wtime();

    // Setting up parellization
    omp_set_num_threads(max_threads);
    int num_threads;
#pragma omp parallel
    {
        num_threads = omp_get_num_threads();
#pragma omp master
        Logger::debug << "Number of threads: " << num_threads << endl;
    }

    // Check time-step for ADI method - the stable time step is completely empirical (i.e. made-up)
    // This isn't really working I think.
    // Only throws error if time step is too large - else has no effect on remaining calculations
    if (inversion_method == InversionMethod::ADI)
    {
        // This check should be included for other inversion_methods
        if (sqrt(1.0 / V.size_y) <= dt)
        {
            Logger::error << "Calculating with ADI, time step " << dt << " is too large." << std::endl;
            exit(EXIT_FAILURE);
        }
        Logger::debug << "Calculating with "
                        << "Diffusion_2D_ADI3" << std::endl;
    }

    // Indexers to keep track of P,R,V,K,L
    double time_current = time_first - dt;

    // variables to show the progress of calculation
    int progress_count, progress_total;

    Matrix3D<double> P_lowerR = P.xSlice(0);
    Matrix3D<double> V_lowerR = V.xSlice(0);
    Matrix3D<double> K_lowerR = K.xSlice(0);
    Matrix3D<double> P_upperR = P.xSlice(R.size_x - 1);
    Matrix3D<double> V_upperR = V.xSlice(R.size_x - 1);
    Matrix3D<double> K_upperR = K.xSlice(R.size_x - 1);

    Matrix3D<double> P_lowerV = P.ySlice(0);
    Matrix3D<double> R_lowerV = R.ySlice(0);
    Matrix3D<double> K_lowerV = K.ySlice(0);
    Matrix3D<double> P_upperV = P.ySlice(V.size_y - 1);
    Matrix3D<double> R_upperV = R.ySlice(V.size_y - 1);
    Matrix3D<double> K_upperV = K.ySlice(V.size_y - 1);

    Matrix3D<double> P_lowerK = P.zSlice(0);
    Matrix3D<double> R_lowerK = R.zSlice(0);
    Matrix3D<double> V_lowerK = V.zSlice(0);
    Matrix3D<double> P_upperK = P.zSlice(K.size_z - 1);
    Matrix3D<double> R_upperK = R.zSlice(K.size_z - 1);
    Matrix3D<double> V_upperK = V.zSlice(K.size_z - 1);

#ifdef MKL_FOUND  
    // initiliaze sparse column indices, row indices and matrix handles
    // for MKL sparse solver for local diffusion
    std::vector<sparse_matrix_t*> sparse_matrix_handles;
    std::vector<std::vector<int>> column_indices;
    std::vector<std::vector<int>> rows_csr;
    if( inversion_method == InversionMethod::MKL && run_local_diffusion)
    {
        std::vector<double> dummy_values;
        initialize_sparse_values(
            V.wxSlice(0,0), K.wxSlice(0,0), Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type,
            DVV.wxSlice(0,0), DVK.wxSlice(0,0), DVK.wxSlice(0,0), DKK.wxSlice(0,0), 
            G_local.wxSlice(0,0), Losses_local.wxSlice(0,0), dt, dummy_values
        );
        matrix_descr descr{
            SPARSE_MATRIX_TYPE_GENERAL, // mkl sparse solve only avaible for general matrices; 
            SPARSE_FILL_MODE_UPPER, // fill mode and unit diagonal have to be set but 
            SPARSE_DIAG_NON_UNIT // are not important for non-triangular matrix types
        };
        for(int j = 0; j < num_threads; j++)
        {
            column_indices.push_back(std::vector<int>());
            rows_csr.push_back(std::vector<int>());
            initialize_sparse_indices(V_size, K_size, column_indices[j], rows_csr[j]);

            sparse_matrix_handles.push_back(new sparse_matrix_t);
            sparse_status_t status = mkl_sparse_d_create_csr(
                sparse_matrix_handles[j], SPARSE_INDEX_BASE_ZERO,
                V_size * K_size, V_size * K_size,
                rows_csr[j].data(), rows_csr[j].data() + 1,
                column_indices[j].data(), dummy_values.data()
            );
            if(status != SPARSE_STATUS_SUCCESS)
            {
                std::cout << "MKL create csr error " << status << '\n';
                exit(EXIT_FAILURE);
            }
            status = mkl_sparse_qr_reorder(*sparse_matrix_handles[j], descr);
            if(status != SPARSE_STATUS_SUCCESS)
            {
                std::cout << "MKL reorder error " << status << '\n';
                exit(EXIT_FAILURE);
            }
        }
    }
#endif
    
    Logger::message << std::endl;
	Logger::writeSeparator();
	Logger::message << std::setw(55) << "Main loop" << std::endl;
	Logger::writeSeparator();

    // Main loop
    // Start time
    for (long int it = it_first; it < it_total; it++)
    {
        // update time by dt every iteration
        time_current = time_first + it * dt;

        if (minimal_output) {
            if (it != it_first) {
                std::cout << "\r";
            }
            std::cout << std::fixed << std::setprecision(6) << std::setw(30) << "Time[" << it << "/" << it_total << "]: " << time_current << " (days)" << std::flush;
        } else {
            Logger::message << std::endl
                            << std::setprecision(15) << "Time[" << it << "/" << it_total << "]: " << time_current << " (days)" << std::endl;
        }

        // Update boundary conditions and diffusion coefficients

        // Update magnetic field (update R)
        if (run_remapping)
        {
            if (L.update(time_current, P, R, V, K))
            {
                // If we've updated L, we most likely need to update Jacobians
                // XXX: Do we need to update Jacobians if we didn't update L?
                G_local.update(time_current, P, R, V, K);
                G_radial.update(time_current, P, R, V, K);
                G_conv.update(time_current, P, R, V, K);

                // If L was updated - interpolate PSD to new L
                progress_count = 0;
                progress_total = P_size * V_size * K_size; // total size of solution matrix
                Logger::debug << "Interpolation to new L (adiabatic transport): ";
                if (not minimal_output) {
                    std::cout << "           ";
                }

                Matrix1D<double> old_L_1d(L_size), PSD_L(L_size), new_L_1d(L_size);
                // Apparently it's not thread-safe
                // #pragma omp parallel for private(iP, iR, iV, iK, iL, PSD_L) shared(progress_total, progress_count) schedule(dynamic,1) collapse(3)
                for (int iP = 0; iP < P_size; iP++) {
                    for (int iV = 0; iV < V_size; iV++) {
                        for (int iK = 0; iK < K_size; iK++) {
                            // show progress % if 0 threads
                            if (not minimal_output) {
                                if (omp_get_thread_num() == 0)
                                {
                                    std::cout << "\b\b\b\b\b\b\b\b\b" << std::setw(8)
                                            << (int)((double)progress_count / progress_total * 100) << "%" << std::flush;
                                } else {
                                    std::cout << "thread" << omp_get_thread_num();
                                }
                            }

                            // 1d slice to get L from matrix4d (P,L,V,K)
                            new_L_1d = L.wyzSlice(iP, iV, iK);
                            old_L_1d = L_copy.wyzSlice(iP, iV, iK);
                            PSD_L = PSD.wyzSlice(iP, iV, iK);

                            // 1d interpolationCubic1D
                            PSD_L = Cubic1D(old_L_1d, PSD_L, new_L_1d);

                            // copy results back into PSD adding the 1d list PSD_L for all values of iP,iV,iK
                            for (int iL = 0; iL < L_size; iL++)
                            {
                                PSD[iP][iL][iV][iK] = PSD_L[iL];
                            }

                            // Progress output - will update P_size * V_size * K_size times
                            // #pragma omp critical
                            progress_count += 1;
                        }
                    }
                }
                if (not minimal_output) {
                    std::cout << "\b\b\b\b\b\b\b\b\b" << std::setw(8) << (int)((double)progress_count / progress_total * 100) << "%" << std::endl;
                }

                // Copy the new L into L_copy for future interpolations
                L_copy = L;
            }
        }

        if (positive_PSD)
        {
            PSD.max_of(1e-21);
        }

#ifdef SAVE_PSD_LOST_CONV
        // Reset loss
        PSD_lost_conv = 0;
#endif

        [[maybe_unused]] bool has_VX_updated = false;
        // Update convection velocities VP and VR and log the maximum absolute values
        if ((3 < P_size || 3 < R_size) && run_convection)
        {
            bool has_VP_updated = VP.update(time_current, P, R, V, K);
            bool has_VR_updated = VR.update(time_current, P, R, V, K);

            if (not minimal_output) {
                Logger::debug << "max(VP) = " << VP.maxabs() << std::endl;
                Logger::debug << "max(VR) = " << VR.maxabs() << std::endl;
            }

            has_VX_updated = has_VP_updated or has_VR_updated;
        }

        // Update convection velocities VV and log the maximum absolute values
        if ((3 < V_size) && run_coulomb_collision)
        {
            VV.update(time_current, P, R, V, K);
            if (not minimal_output) {
                Logger::debug << "max(VV) = " << VV.maxabs() << std::endl;
            }
        }

        // Diffusion coefficients
        if ((3 < L_size) && (run_radial_diffusion))
        {
            DLL.update(time_current, P, L, V, K);
            if (not minimal_output) {
                Logger::debug << "max(DLL) = " << DLL.maxabs() << std::endl;
            }
        }

        [[maybe_unused]] bool updated_dvv, updated_dvk, updated_dkk, updated_local_loss;
        if ((3 < V_size && 3 < K_size) && (run_local_diffusion))
        {
            updated_dvv = DVV.update(time_current, P, R, V, K);
            updated_dvk = DVK.update(time_current, P, R, V, K);
            updated_dkk = DKK.update(time_current, P, R, V, K);

            if (not minimal_output) {
                Logger::debug << "max(DVV) = " << DVV.maxabs() << std::endl;
                Logger::debug << "max(DVK) = " << DVK.maxabs() << std::endl;
                Logger::debug << "max(DKK) = " << DKK.maxabs() << std::endl;
            }
        }

        // Sources and losses
        Sources.update(time_current, P, R, V, K);
        updated_local_loss = Losses_local.update(time_current, P, R, V, K);
        Losses_radial.update(time_current, P, R, V, K);
        Losses_conv.update(time_current, P, R, V, K);

        if (density_saturation == DensitySaturation::WithTimescale) {
            Sources = SaturationDensity - PSD;
            for (int iP = 0; iP < P_size; ++iP)
                for (int iR = 0; iR < R_size; ++iR)
                    for (int iV = 0; iV < V_size; ++iV)
                        for (int iK = 0; iK < K_size; ++iK) {
                            Sources[iP][iR][iV][iK] = Sources[iP][iR][iV][iK] / SaturationTimescale[iP][iR][iV][iK];
                            if (Sources[iP][iR][iV][iK] < 1e-21)
                                Sources[iP][iR][iV][iK] = 1e-21;
                        }
        } else if (density_saturation == DensitySaturation::WithoutTimescale) {
            for (int iP = 0; iP < P_size; ++iP)
                for (int iR = 0; iR < R_size; ++iR)
                    for (int iV = 0; iV < V_size; ++iV)
                        for (int iK = 0; iK < K_size; ++iK)
                            if (PSD[iP][iR][iV][iK] > SaturationDensity[iP][iR][iV][iK])
                                PSD[iP][iR][iV][iK] = SaturationDensity[iP][iR][iV][iK];

        } else if (density_saturation == DensitySaturation::LimitSource) {
            for (int iP = 0; iP < P_size; ++iP)
                for (int iR = 0; iR < R_size; ++iR)
                    for (int iV = 0; iV < V_size; ++iV)
                        for (int iK = 0; iK < K_size; ++iK) {
                            if (PSD[iP][iR][iV][iK] >= SaturationDensity[iP][iR][iV][iK]) {
                                Sources[iP][iR][iV][iK] = 1e-21;
                            }
                        }
        }

        // Boundary conditions
        // By default - it's constant PSD
        if (R_size > 3)
        {
            // Put PSD-slice into the boundary condition. It will be updated from the ini-file, if there is something in the ini-file.
            // R_LBC.original_arr =(PSD.xSlice(0));
            Rl_BC.update(time_current, P_lowerR, V_lowerR, K_lowerR);
            // R_UBC.original_arr =(PSD.xSlice(PSD.size_x - 1));
            Ru_BC.update(time_current, P_upperR, V_upperR, K_upperR);
        }
        if (L_size > 3) {
            // L_LBC.original_arr =(PSD.xSlice(0));
            Ll_BC.update(time_current, P_lowerR, V_lowerR, K_lowerR);
            // L_UBC.original_arr =(PSD.xSlice(PSD.size_x - 1));
            Lu_BC.update(time_current, P_upperR, V_upperR, K_upperR);
        }
        if (V_size > 3) {
            // V_LBC.original_arr =(PSD.ySlice(0));
            Vl_BC.update(time_current, P_lowerV, R_lowerV, K_lowerV);
            // V_UBC.original_arr =(PSD.ySlice(PSD.size_y - 1));
            Vu_BC.update(time_current, P_upperV, R_upperV, K_upperV);
        }
        if (K_size > 3) {
            // K_LBC.original_arr =(PSD.zSlice(0));
            Kl_BC.update(time_current, P_lowerK, R_lowerK, V_lowerK);
            // K_UBC.original_arr =(PSD.zSlice(PSD.size_z - 1));
            Ku_BC.update(time_current, P_upperK, R_upperK, V_upperK);
        }

        // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Check that DQ2Q3^2 < DQ2Q2 * DQ3Q3
        // The computation will be unstable otherwise

        // checking if DVV^2 - min(DVK^2) is greater than 0, if not log and exit
        double check = (DVV.times(DKK) - DVK.times(DVK)).min();
        if (check < 0)
        {
            Logger::error << "Fatal error: (DVV*DKK - DVK*DVK) = " << check << " < 0." << std::endl;
            Logger::error << "The computation would be unstable." << std::endl;
            exit(EXIT_FAILURE);
        }

    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// COLOUMB CONVECTION WITH CONVECTION IN P, R and V
    if ( (3 < P_size && 3 < R_size && 3 < V_size) && run_convection && run_coulomb_collision ) {
            progress_count = 0;
            progress_total = K_size;
            Logger::debug << "3D Convection (P,R,V):" << endl;
    
            if (not minimal_output) {
                cout << "           ";
            }

#pragma omp parallel shared(progress_total, progress_count) 
        {
            // define the slices needed for convection
            Matrix3D<double> PSD_PRV(P_size, R_size, V_size);
            Matrix2D<double> plow_boundary_slice(R_size, V_size);
            Matrix2D<double> pup_boundary_slice(R_size, V_size);
            Matrix2D<double> rlow_boundary_slice(P_size, V_size);
            Matrix2D<double> rup_boundary_slice(P_size, V_size);
            Matrix2D<double> vlow_boundary_slice(P_size, R_size);
            Matrix2D<double> vup_boundary_slice(P_size, R_size);

            Matrix3D<double> pgrid_slice(P_size, R_size, V_size);
            Matrix3D<double> rgrid_slice(P_size, R_size, V_size);
            Matrix3D<double> vgrid_slice(P_size, R_size, V_size);
            Matrix3D<double> vp_slice(P_size, R_size, V_size);
            Matrix3D<double> vr_slice(P_size, R_size, V_size);
            Matrix3D<double> vv_slice(P_size, R_size, V_size);
            Matrix3D<double> lossconv_slice(P_size, R_size, V_size);
            Matrix3D<double> G_conv_slice(P_size, R_size, V_size);

            // sources are set to zero and not updated during the loop
            Matrix3D<double> source_slice = Sources.zSlice(0) * 0.0;

#pragma omp for schedule(dynamic,1)
             for (int iK = 0; iK < K_size; iK++) {

                // Output current progress percentage when number of threads = 0
                if (not minimal_output) {
                    if (omp_get_thread_num() == 0) {
                        cout << "\b\b\b\b\b\b\b\b\b" << setw(8)
                            << (int) ((double) progress_count / progress_total * 100) << "%" << flush;
                    }
                }

                // update all slices for convection
                Pl_BC.zSlice(plow_boundary_slice, iK);
                Pu_BC.zSlice(pup_boundary_slice, iK);
                Rl_BC.zSlice(rlow_boundary_slice, iK);
                Ru_BC.zSlice(rup_boundary_slice, iK);
                Vl_BC.zSlice(vlow_boundary_slice, iK);
                Vu_BC.zSlice(vup_boundary_slice, iK);

                P.zSlice(pgrid_slice, iK);
                R.zSlice(rgrid_slice, iK);
                V.zSlice(vgrid_slice, iK);
                VP.zSlice(vp_slice, iK);
                VR.zSlice(vr_slice, iK);
                VV.zSlice(vv_slice, iK);
                Losses_conv.zSlice(lossconv_slice, iK);
                PSD.zSlice(PSD_PRV, iK);
                G_conv.zSlice(G_conv_slice, iK);

                Convection_3D(
                    PSD_PRV, pgrid_slice, rgrid_slice, vgrid_slice, P_size, R_size, V_size,
                    plow_boundary_slice, pup_boundary_slice,
                    rlow_boundary_slice, rup_boundary_slice,
                    vlow_boundary_slice, vup_boundary_slice,
                    Pl_BC_type, Pu_BC_type, Rl_BC_type, Ru_BC_type, Vl_BC_type, Vu_BC_type,
                    vp_slice, vr_slice, vv_slice, 
                    source_slice, lossconv_slice,
                    G_conv_slice, dt);     

                // copy results back into PSD adding the 2d list PSD_PR for all values of iV,iK
                for (int iP = 0; iP < P_size; iP++) {
                    for (int iR = 0; iR < R_size; iR++) {
                        for (int iV = V_size - 1; iV >= 0; iV--) {
                            PSD[iP][iR][iV][iK] = PSD_PRV[iP][iR][iV];
                        }
                    }
                }
            }

#pragma omp critical
            progress_count += 1;
        }
            // Output final progress (it should be 100%)
            if (not minimal_output) {
                cout << "\b\b\b\b\b\b\b\b\b" << setw(8) << (int) ((double) progress_count / progress_total * 100) << "%" << endl;
            }
#pragma omp master
            {
            }
	}
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Convection, for each V and K (therefore for each mu and J)
        if ((3 < P_size || 3 < R_size) && run_convection && not run_coulomb_collision)
        {
            progress_count = 0;
            progress_total = V_size * K_size;
            Logger::debug << "Convection:" << std::endl;
            if (not minimal_output) {
                std::cout << "           ";
            }

#pragma omp parallel shared(progress_total, progress_count)
            {
                // define the slices needed for convection
                Matrix2D<double> PSD_PR(P_size, R_size);
                Matrix1D<double> plow_boundary_slice(R_size);
                Matrix1D<double> pup_boundary_slice(R_size);
                Matrix1D<double> rlow_boundary_slice(P_size);
                Matrix1D<double> rup_boundary_slice(P_size);
                Matrix2D<double> pgrid_slice(P_size, R_size);
                Matrix2D<double> rgrid_slice(P_size, R_size);
                Matrix2D<double> vp_slice(P_size, R_size);
                Matrix2D<double> vr_slice(P_size, R_size);
                Matrix2D<double> lossconv_slice(P_size, R_size);
                Matrix2D<double> sources_slice(P_size, R_size);
                Matrix2D<double> G_conv_slice(P_size, R_size);

#ifdef SAVE_PSD_LOST_CONV
                Matrix2D<double> PSD_lost_PR(P_size, R_size);
#endif

                // Looping it backward allows to speed-up the multithread simulation
                // due to the highest energies being the slowest to calculate - calculating highest energy first

// If you have OPENMP 3.0 with collapse functionality you can use this. (Visual Studios does not)
#if defined _OPENMP && _OPENMP >= 200711
#pragma omp for schedule(dynamic, 1) collapse(2)
                for (int iV = V_size - 1; iV >= 0; iV--)
                {
                    for (int iK = 0; iK < K_size; iK++)
                    {
#else
#pragma omp for schedule(dynamic, 1)
                for (long index = (V_size * K_size) - 1; index >= 0; index--) {
                    {  // Double brackes are needed for compatability with OPENMP 3.0 code.
                        // Manually collapsing into a single for loop in order to optimize
                        int iV = index % V_size;
                        int iK = (K_size - 1) - (index / V_size);
#endif
                        // Output current progress percentage when number of threads = 0
                        if (not minimal_output) {
                            if (omp_get_thread_num() == 0) {
                                std::cout << "\b\b\b\b\b\b\b\b\b" << std::setw(8)
                                        << (int)((double)progress_count / progress_total * 100) << "%" << std::flush;
                            }
                        }

                        // update all slices for convection
                        Pl_BC.yzSlice(plow_boundary_slice, iV, iK);
                        Pu_BC.yzSlice(pup_boundary_slice, iV, iK);
                        Rl_BC.yzSlice(rlow_boundary_slice, iV, iK);
                        Ru_BC.yzSlice(rup_boundary_slice, iV, iK);

                        P.yzSlice(pgrid_slice, iV, iK);
                        R.yzSlice(rgrid_slice, iV, iK);
                        VP.yzSlice(vp_slice, iV, iK);
                        VR.yzSlice(vr_slice, iV, iK);
                        Losses_conv.yzSlice(lossconv_slice, iV, iK);
                        Sources.yzSlice(sources_slice, iV, iK);
                        PSD.yzSlice(PSD_PR, iV, iK);
                        G_conv.yzSlice(G_conv_slice, iV, iK);

#ifdef SAVE_PSD_LOST_CONV
                        PSD_lost_conv.yzSlice(PSD_lost_PR, iV, iK);

                        Convection_2D(
                            PSD_PR, PSD_lost_PR, pgrid_slice, rgrid_slice, P_size, R_size,
                            plow_boundary_slice, pup_boundary_slice,
                            rlow_boundary_slice, rup_boundary_slice,
                            Pl_BC_type, Pu_BC_type, Rl_BC_type, Ru_BC_type,
                            vp_slice, vr_slice,
                            sources_slice, lossconv_slice, G_conv_slice,
                            dt);
#else
                        Convection_2D(
                            PSD_PR, pgrid_slice, rgrid_slice, P_size, R_size,
                            plow_boundary_slice, pup_boundary_slice,
                            rlow_boundary_slice, rup_boundary_slice,
                            Pl_BC_type, Pu_BC_type, Rl_BC_type, Ru_BC_type,
                            vp_slice, vr_slice,
                            sources_slice, lossconv_slice, G_conv_slice,
                            dt);
#endif

                        // copy results back into PSD adding the 2d list PSD_PR for all values of iV,iK
                        for (int iP = 0; iP < P_size; iP++) {
                            for (int iR = 0; iR < R_size; iR++) {
                                PSD[iP][iR][iV][iK] = PSD_PR[iP][iR];
#ifdef SAVE_PSD_LOST_CONV
                                PSD_lost_conv[iP][iR][iV][iK] = PSD_lost_PR[iP][iR];

#endif
                            }
                        }

#pragma omp critical
                        progress_count += 1;
                    }
                }
            }
            // Output final progress (it should be 100%)
            if (not minimal_output) {
                std::cout << "\b\b\b\b\b\b\b\b\b" << std::setw(8) << (int)((double)progress_count / progress_total * 100) << "%" << std::endl;
            }
#pragma omp master
            {
                if (Vl_BC_from_convection && (Vl_BC_type == BoundaryConditionType::ConstantValue)) {  // rewrite boundary conditions at lower V
                    Vl_BC = PSD.ySlice(0);
                    if (not minimal_output) {
                        Logger::debug << "Vl_BC from convection are used: max(Vl_BC) = " << Vl_BC.max() << std::endl;
                    }
                }
                if (Vu_BC_from_convection && (Vu_BC_type == BoundaryConditionType::ConstantValue)) {  // rewrite boundary conditions at lower V
                    Vu_BC = PSD.ySlice(V_size - 1);
                    if (not minimal_output) {
                        Logger::debug << "Vu_BC from convection are used: max(Vu_BC) = " << Vu_BC.max() << std::endl;
                    }
                }
            }
        }

        if (positive_PSD)
        {
            PSD.max_of(1e-21);
        }

        // RADIAL DIFFUSION STEP

        // ADDED FOR TESTING
        //  PSD.writeToFile(to_string(int(it / output_step)) +  "PSD_before_radial.plt");

        if ((L_size >= 3) && (run_radial_diffusion))
        {
            progress_count = 0;
            progress_total = P_size * V_size * K_size;  // total size of solution matrix
            Logger::debug << "Radial diffusion:" << std::endl;
            if (not minimal_output) {
                std::cout << "           ";
            }

#pragma omp parallel shared(progress_total, progress_count)
            {
                Matrix1D<double> PSD_L(L_size);
                Matrix1D<double> lgrid_slice(L_size);
                Matrix1D<double> dll_slice(L_size);
                Matrix1D<double> source_slice(L_size);
                Matrix1D<double> gradial_slice(L_size);
                Matrix1D<double> loss_slice(L_size);

#if defined _OPENMP && _OPENMP >= 200711
#pragma omp for schedule(dynamic, 1) collapse(3)
                for (int iP = 0; iP < P_size; iP++) {
                    for (int iV = 0; iV < V_size; iV++) {
                        for (int iK = 0; iK < K_size; iK++) {
#else
#pragma omp for schedule(dynamic, 1)
                for (long index = 0; index < P_size * V_size * K_size; index++) {
                    {
                        {
                            int iP = index / (V_size * K_size);
                            int iV = (index / K_size) % V_size;
                            int iK = index % K_size;
#endif
                            // print percentage done
                            if (not minimal_output) {
                                if (omp_get_thread_num() == 0) {
                                    std::cout << "\b\b\b\b\b\b\b\b\b" << std::setw(8) << (int)((double)progress_count / progress_total * 100) << "%" << std::flush;
                                }
                            }

                            // 1d slice
                            PSD.wyzSlice(PSD_L, iP, iV, iK);
                            L.wyzSlice(lgrid_slice, iP, iV, iK);
                            DLL.wyzSlice(dll_slice, iP, iV, iK);
                            Sources.wyzSlice(source_slice, iP, iV, iK);
                            G_radial.wyzSlice(gradial_slice, iP, iV, iK);
                            Losses_radial.wyzSlice(loss_slice, iP, iV, iK);

                            source_slice *= radial_sources;

                            // 1d diffusion
                            Diffusion_1D(
                                PSD_L, lgrid_slice, L_size, Ll_BC[iP][iV][iK], Lu_BC[iP][iV][iK],
                                Ll_BC_type, Lu_BC_type, dll_slice, gradial_slice,
                                source_slice, loss_slice, dt);

                            // copy results back
                            for (int iL = 0; iL < L_size; iL++) {
                                PSD[iP][iL][iV][iK] = PSD_L[iL];
                            }

#pragma omp critical
                            progress_count += 1;
                        }
                    }
                }
            }
            // ADDED FOR TESTING
            //  PSD.writeToFile(to_string(int(it / output_step)) +  "PSD_after_radial.plt");
            if (not minimal_output) {
                std::cout << "\b\b\b\b\b\b\b\b\b" << std::setw(8) << (int)((double)progress_count / progress_total * 100) << "%" << std::endl;
            }
            // #pragma omp master
            //             {
            //             if((Vl_BC_from_convection == "true") && (Vl_BC_type == "BCT_CONSTANT_VALUE")) { //rewrite boundary conditions at lower V
            //                 Vl_BC = PSD.ySlice(1); // NOTE: low-V boundary conditions are taken from the next-to-lower-V slice after radial diffusion.
            //                                        // This eliminates unrealistic gradient there, but it is probably NOT A CORRECT IMPLEMENTATION.
            //                 //std::cout << "Vl_BC after radial diffusion are used: max(Vl_BC) = " << Vl_BC.max() << endl;
            //                 Logger::message << "Vl_BC after radial diffusion are used: max(Vl_BC) = " << Vl_BC.max() << endl;
            //             }
            //             }
        }

        if (positive_PSD)
        {
            PSD.max_of(1e-21);
        }

        // ADDED FOR TESTING
        //  PSD.writeToFile(to_string(int(it / output_step)) +  "PSD_after_radial.plt");

        // LOCAL DIFFUSION
        if ((V_size >= 3 && K_size >= 3) && (run_local_diffusion))
        {
            int number_of_skipped_points = 0;
            progress_count = 0;
            progress_total = P_size * R_size;
            Logger::debug << "Local diffusion: " << std::endl;
            if (not minimal_output) {
                std::cout << "           ";
            }

#ifdef LU_CACHING
            recompute_lu_decomp = updated_local_loss or updated_dvv or updated_dvk or updated_dkk;
            if (recompute_lu_decomp) {
                Logger::debug << "Recomputing LU decomposition!" << std::endl;
            }
#endif

            Matrix2D<double> PSD_IK(V_size, K_size);

#if defined _OPENMP && _OPENMP >= 200711
#pragma omp parallel for private(PSD_IK) shared(progress_total, progress_count, number_of_skipped_points) schedule(dynamic, 1) collapse(2)
            for (int iP = 0; iP < P_size; iP++)
            {
                for (int iR = 0; iR < R_size; iR++)
                {
#else
#pragma omp parallel for private(PSD_IK) shared(progress_total, progress_count, number_of_skipped_points) schedule(dynamic, 1)
            for (int index = 0; index < P_size * R_size; index++)
            {
                {
                    int iP = index / R_size;
                    int iR = index % R_size;
#endif
                    // If we dont want to include the boundary in the L/R axis
                    // we can exclude those points by setting include_boundary
                    // to false in matlab main file or parameter.ini file
                    if ((not include_boundary) && (iR == 0 || iR == R_size - 1))
                    {
                        continue;
                    }

                    if (not minimal_output) {
                        if (omp_get_thread_num() == 0) {
                            std::cout << "\b\b\b\b\b\b\b\b\b" << std::setw(8)
                                    << (int)((double)progress_count / progress_total * 100) << "%" << std::flush;
                        }
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
                    if (inversion_method == InversionMethod::Lapack) {
#ifdef LU_CACHING
                        Diffusion_2D(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
                                    Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR),  // P, R, K
                                    Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR),  // P, R, I
                                    Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
                                    DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
                                    Sources.wxSlice(iP, iR) * local_sources, Losses_local.wxSlice(iP, iR), dt, sub_dt_diffusion,
                                    lu_cache[iP][iR], index_cache[iP][iR], recompute_lu_decomp);
#else
                        Diffusion_2D(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
                                    Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR),  // P, R, K
                                    Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR),  // P, R, I
                                    Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
                                    DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
                                    Sources.wxSlice(iP, iR) * local_sources, Losses_local.wxSlice(iP, iR), dt, sub_dt_diffusion);
#endif
                    }
                    // Currently setup to calculate 2d Diffusion using Diffusion_2D_ADI3
                    // Can change to Diffusion_2D_ADI1 or Diffusion_2D_ADI2 for different methods of inversion
                    else if (inversion_method == InversionMethod::ADI) {
                        Diffusion_2D_ADI3(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
                                          Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR),  // P, R, K
                                          Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR),  // P, R, I
                                          Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
                                          DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
                                          Sources.wxSlice(iP, iR) * local_sources, Losses_local.wxSlice(iP, iR), dt);
                    } else if (inversion_method == InversionMethod::ADI2) {
                        Diffusion_2D_ADI2(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
                                          Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR),  // P, R, K
                                          Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR),  // P, R, I
                                          Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
                                          DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
                                          Sources.wxSlice(iP, iR) * local_sources, Losses_local.wxSlice(iP, iR), dt);
                    } else if (inversion_method == InversionMethod::ADI1) {
                        Diffusion_2D_ADI1(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
                                          Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR),  // P, R, K
                                          Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR),  // P, R, I
                                          Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
                                          DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
                                          Sources.wxSlice(iP, iR) * local_sources, Losses_local.wxSlice(iP, iR), dt);
                    } else {
                        Logger::error << "Error: Unknown inversion method!" << std::endl;
                    }

                    // copy results back
                    for (int iV = 0; iV < V_size; iV++)
                    {
                        for (int iK = 0; iK < K_size; iK++)
                        {
                            PSD[iP][iR][iV][iK] = PSD_IK[iV][iK];
                        }
                    }
                }
            }

            if (not minimal_output) {
                std::cout << "\b\b\b\b\b\b\b\b\b" << std::setw(8) << (int)((double)progress_count / progress_total * 100) << "%" << std::endl;
                std::cout << "Number of skipped points: " << (int)((double)number_of_skipped_points / progress_total * 100) << "%" << std::endl;
            }
        }

        int number_of_negative_points = 0;
        if (positive_PSD)
        {
            PSD.max_of(1e-21);
        }
        else
        {
            for (int iP = 0; iP < P_size; iP++)
            {
                for (int iR = 0; iR < R_size; iR++)
                {
                    for (int iV = 0; iV < V_size; iV++)
                    {
                        for (int iK = 0; iK < K_size; iK++)
                        {
                            if (PSD[iP][iR][iV][iK] < 0)
                                number_of_negative_points++;
                        }
                    }
                }
            }
        }
        // std::cout << "Number of negative points: " << number_of_negative_points << endl;
        Logger::debug << std::endl
                        << "Number of negative points: " << number_of_negative_points << " of " << P_size * R_size * V_size * K_size << std::endl;

#ifdef DATA_ASSIMILATION
        daManagerConvection.assimilate(time_current, PSD, VP, VR, has_VX_updated, Sources, SaturationDensity, dt);
        if (positive_PSD) {
            PSD.max_of(1e-21);
        }
#endif

        // Output the PSD data for each timestep into the output folder
        if ((it % output_step) == 0)
        {
            // Wait until the writing of the last output file is finished
            output_writer.wait();

            // PSD_filename already defined
            // std::ostringstream PSD_filename;
            PSD_filename << outputFolder << "PSD_" << std::setw(5) << std::setfill('0') << int(it / output_step);

            Logger::debug << std::endl
                            << "Writing results: " << PSD_filename.str() << std::endl;

            output_writer = std::async(std::launch::async, write_PSD_output, io_method, PSD_filename.str(), time_current, PSD, PSD_time_to_lst, outputFolder);

#ifdef SAVE_PSD_LOST_CONV
            std::ostringstream PSD_lost_filename, time_string;
            time_string.str("");
            time_string << time_current;

            PSD_lost_filename.str("");
            PSD_lost_filename << outputFolder << "lost_conv_PSD_" << std::setw(5) << std::setfill('0') << int(it / output_step);
            Logger::debug << std::endl
                            << "Writing results: " << PSD_lost_filename.str() << std::endl;
            PSD_lost_conv.writeToAnyFile(PSD_lost_filename.str(), io_method, time_string.str());
#endif
        }

    }

#ifdef MKL_FOUND
    for(auto& csr_handle : sparse_matrix_handles)
    {
        mkl_sparse_destroy(*csr_handle);
        // delete csr_handle; //mkl_sparse_destroy seems to delete the pointer
    }
#endif

#ifdef LU_CACHING
    for(auto& p_slice : lu_cache)
    {
        for(auto& pr_slice : p_slice)
        {
            delete[] pr_slice;
        }
    }
    for(auto& p_slice : index_cache)
    {
        for(auto& pr_slice : p_slice)
        {
            delete[] pr_slice;
        }
    }
#endif

    // logger records if everything went correctly
    Logger::message << "\nProgram was terminated correctly." << std::endl;
    Logger::message << "Wall-clock time: " << std::fixed << std::setprecision(2) << (omp_get_wtime() - wall_timer) << " sec; ";
    Logger::message << "CPU time: " << (float)(clock() - start_time) / CLOCKS_PER_SEC << " sec." << std::endl;

    Logger::deleteInstance();
    return 0;
}
