/*
 * Diffusion + convection code, VERB4D.
 *
 * The code solves the Fokker-Planck equation with convection terms:
 * \f[
 * \frac{\partial f}{\partial t} =
 * v_{\Phi} \frac{\partial f}{\partial \pPi} + v_{R} \frac{\partial f}{\partial R} +
 * \frac{1}{G} \frac{\partial}{\partial L} G D_{LL} \frac{\partial f}{\partial L} +
 * \frac{1}{G} \frac{\partial}{\partial V} G (D_{VV} \frac{\partial f}{\partial V} + D_{VK} \frac{\partial f}{\partial K}) +
 * \frac{1}{G} \frac{\partial}{\partial K} G (D_{KV} \frac{\partial f}{\partial V} + D_{KK} \frac{\partial f}{\partial K}) +
 * Losses * f + Sources
 * \f],
 * where f is a Phase Space Density (PSD); Φ is MLT angle, R is radial distance from Earth center, V≡μ∙(K+0.5)2 is a combination of adiabatic
 * invariants μ and K≡J/√8mμ; L is the third adiabatic invariant; D are bounce-averaged diffusion coefficients; v_Φ, v_R are drift velocities; G is the
 * Jacobin of the transformation from (μ,J,L) to (V,K,L)
 */

/*
 * Books, required to understand the code:
 * Schulz and Lanzerotti, 1974
 * Stroustrup C++
 *
 * Continue reading ONLY if you are familiar yourself with these books.
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <time.h>

#include <math.h>

using namespace std;

#include "Matrix.h"

#include "Diffusion_1D.h"

// Different solution methods for 2D diffusion
#include "Diffusion_2D.h"

// Straightforward ADI
//#include "Diffusion_ADI1.h"
// Fuliang Xiao's ADI
//#include "Diffusion_ADI2.h"
// Jihye Shin and Sungsoo S. Kim (2008)'s ADI - most stable
#include "Diffusion_ADI3.h"

#include "Convection_2D.h"

#include "ReadInitialData.h"

#include "UpdatableMatrix.h"

#include "Interpolation.h"

#include <omp.h>

//#define DEBUG_MODE

#if defined(_WIN32) || defined(_WIN64)
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif


// Everything below these values will be considered to be zero for computation purpose and will not be calculated
// XXX: It should be a parameter
#define min_PSD 1e-10
#define min_V 1e-10
#define min_Dxx 1e-10


int main(int argc, char* argv[]) {

	cout << "Compilation time: " << __DATE__ << ", " << __TIME__ << endl;

	// Variables
	// Grid, 4D:
	// P - local time
	// R - radial distance
	// V - invariant = mu / (K + const) ^ 2, where mu is the first adiabatic invariant
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

	// Convection (advection																																																		) velocities, 4D arrays
	UpdatableMatrix<Matrix4D<double> > VR, VP;

	// Additional sources and losses
	UpdatableListMatrix<Matrix4D<double> > Sources, Losses;

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

	double dt = 1.0 / 24;
	double time_output = 1.0 / 24;
	double time_total = 4;
	double time_first = 0;
	int max_threads = omp_get_num_threads();

	string inputFolder = "./VERB4D_input/";
	string outputFolder = "./VERB4D_output/";
	string inversion_method = "Lapack";

	// Read all the inputs
	// TODO: create structures for parameters, so there would be no need to pass them one by one
	ReadInitialData(inputFolder, argc, argv, time_total, dt, time_output, time_first, it_first, max_threads,
			inversion_method, PSD,
			P, R, V, K, L,
			P_size, R_size, V_size, K_size, L_size, Pl_BC, Pu_BC, Rl_BC, Ru_BC,
			Vl_BC, Vu_BC, Kl_BC, Ku_BC, Ll_BC, Lu_BC, Pl_BC_type, Pu_BC_type, Rl_BC_type, Ru_BC_type, Vl_BC_type,
			Vu_BC_type, Kl_BC_type, Ku_BC_type, Ll_BC_type, Lu_BC_type, DLL, DVV, DKK, DVK, VP, VR, G_local, G_radial,
			Sources, Losses);

	// Copy L-star so we can later interpolate PSD to a new L-star,
	// to account for adiabatic transport if L-star changes
	L_copy = L;

	it_total = it_first + round(double(time_total) / dt);
	output_step = round(double(time_output) / dt);
	if (output_step < 1)
		output_step = 1;

	cout << "Total time " << time_total << ". Time step " << dt << " (" << it_total << " steps)." << endl;
	cout << "Output each " << output_step << " step. " << endl;

	// Save initial conditions
	PSD.writeToFile(outputFolder + "PSD0.dat", P, R, V, K);

	// Output zero step
	ostringstream PSD_filename, time_string;
	time_string.precision(5);
	time_string.setf(ios::fixed);
	PSD_filename << outputFolder << "PSD_" << setw(5) << setfill('0') << 0 << ".dat";
	cout << "Writing results: " << PSD_filename.str() << endl;
	time_string.str("");
	time_string << time_first;
	PSD.writeToFile(PSD_filename.str(), time_string.str());

	int radial_losses, local_losses;

	// When to apply loss term:
	// It's better to apply it during pitch-angle diffusion, unless we don't have any pitch-angle diffusion
	// In that case we apply it during radial diffusion
	if (V_size >= 3 && K_size >= 3) {
		radial_losses = 0;
		local_losses = 1;
	} else {
		radial_losses = 1;
		local_losses = 0;
	}

	std::clock_t start_time, tic, tok;
	start_time = clock();
	double wall_timer = omp_get_wtime();

	omp_set_num_threads(max_threads);
#pragma omp parallel
	{
#pragma omp master
		cout << "Number of threads: " << omp_get_num_threads() << endl;
	}

	// Check time-step for ADI method - the stable time step is completely empirical (i.e. made-up)
	// This isn't really working I think.
	if (inversion_method == "ADI") {
		// if (dt > 0.1/24 + 0.001) {
		if (dt > sqrt(1.0 / V.size_y)) {
			cout << "Calculating with ADI, time step " << dt << " is too large." << endl;
			abort();
		} else {
			cout << "Calculating with " << "Diffusion_2D_ADI3" << endl;
		}
	} else {
		cout << inversion_method << endl;
	}

	int iP, iR, iV, iK, iL;
	double time;

	// variables to show the progress of calculation
	int progress_count, progress_total;

	// Main loop
	// Start time
	for (long int it = it_first; it < it_total; it++) {
		time = time_first + it * dt;
		cout << endl << "Time[" << it << "/" << it_total << "]: " << time << " (days)" << endl;

		// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Update boundary conditions and diffusion coefficients

		// Update magnetic field (update R)
		if (L.update(time, P, R, V, K)) {
			// If we've updated L, we most likely need to update Jacobians
			// XXX: Do we need to update Jacobians if we didn't update L?
			G_local.update(time, P, R, V, K);
			G_radial.update(time, P, R, V, K);

			// If L was updated - interpolate PSD to new L
			progress_count = 0;
			progress_total = P_size * V_size * K_size;
			cout << "Interpolation to new L (adiabatic transport)  " << "         ";

			Matrix1D<double> old_L_1d(L_size), PSD_L(L_size), new_L_1d(L_size);
// Aparently it's not thread-safe
// #pragma omp parallel for private(iP, iR, iV, iK, iL, PSD_L) shared(progress_total, progress_count) schedule(dynamic,1) collapse(3)
			for (iP = 0; iP < P_size; iP++) {
				for (iV = 0; iV < V_size; iV++) {
					for (iK = 0; iK < K_size; iK++) {

						if (omp_get_thread_num() == 0)
							cout << "\b\b\b\b\b\b\b\b\b" << setw(8)
									<< (int) ((double) progress_count / progress_total * 100) << "\%" << flush;

						// 1d slice
						new_L_1d = L.wyzSlice(iP, iV, iK);
						old_L_1d = L_copy.wyzSlice(iP, iV, iK);
						PSD_L = PSD.wyzSlice(iP, iV, iK);

						// 1d interpolationCubic1D
						PSD_L = Cubic1D(old_L_1d, PSD_L, new_L_1d);

						// copy results back
						for (iL = 0; iL < L_size; iL++)
							PSD[iP][iL][iV][iK] = PSD_L[iL];

//#pragma omp critical // Progress output
						progress_count += 1;

					}
				}
			}
			cout << "\b\b\b\b\b\b\b\b\b" << setw(8) << (int) ((double) progress_count / progress_total * 100) << "\%"
					<< endl;

			// Copy the new L into L_copy for future interpolations
			L_copy = L;
		}

		// Update convection velocities
		if (P_size >= 3 || R_size >= 3) {
			VP.update(time, P, R, V, K);
			cout << "max(VP) = " << VP.maxabs() << endl;
			VR.update(time, P, R, V, K);
			cout << "max(VR) = " << VR.maxabs() << endl;
		}

		// Diffusion coefficients
		if (L_size >= 3) {
			DLL.update(time, P, L, V, K);
			cout << "max(DLL) = " << DLL.maxabs() << endl;
		}

		if (V_size >= 3 && K_size >= 3) {
			DVV.update(time, P, R, V, K);
			cout << "max(DVV) = " << DVV.maxabs() << endl;
			DVK.update(time, P, R, V, K);
			cout << "max(DVK) = " << DVK.maxabs() << endl;
			DKK.update(time, P, R, V, K);
			cout << "max(DKK) = " << DKK.maxabs() << endl;
		}

		// Sources and losses
		Sources.update(time, P, R, V, K);
		Losses.update(time, P, R, V, K);

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
		if (true) {
			double check = (DVV.times(DKK) - DVK.times(DVK)).min();
			if (check < 0) {
				cout << "Fatal error: (DVV*DKK - DVK*DVK) = " << check << " < 0." << endl;
				cout << "The computation would be unstable." << endl;
				abort();
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Convection, for each V and K (therefore for each mu and J)
		if (P_size >= 3 && R_size >= 3) {
			progress_count = 0;
			progress_total = V_size * K_size;
			cout << "Convection  " << "         ";

			Matrix2D<double> PSD_PR(P_size, R_size);
			// omp-paraller loop
#pragma omp parallel for private(iP, iR, iV, iK, iL, PSD_PR) shared(progress_total, progress_count) schedule(dynamic,1) collapse(2)
			for (iV = V_size-1; iV >= 0; iV--) { 	// Looping it backward allows to speed-up the multithread simulation
												// due to the highest energies being the slowest to calculate
				for (iK = 0; iK < K_size; iK++) {

					// Output current progress
					if (omp_get_thread_num() == 0)
						cout << "\b\b\b\b\b\b\b\b\b" << setw(8)
								<< (int) ((double) progress_count / progress_total * 100) << "\%" << flush;

					// 2d slice of PSD to pass to the Convection calculation function
					PSD_PR = PSD.yzSlice(iV, iK);

					Convection_2D(PSD_PR, P.yzSlice(iV, iK), R.yzSlice(iV, iK), P_size, R_size, Pl_BC.yzSlice(iV, iK),
							Pu_BC.yzSlice(iV, iK), // R, I, K
							Rl_BC.yzSlice(iV, iK), Ru_BC.yzSlice(iV, iK), // P, I, K
							Pl_BC_type, Pu_BC_type, Rl_BC_type, Ru_BC_type, VP.yzSlice(iV, iK), VR.yzSlice(iV, iK),
							Sources.yzSlice(iV, iK) * 0, Losses.yzSlice(iV, iK) * 0, dt, min_PSD, min_V);

					// copy results back
					for (iP = 0; iP < P_size; iP++)
						for (iR = 0; iR < R_size; iR++)
							PSD[iP][iR][iV][iK] = PSD_PR[iP][iR];

#pragma omp critical // Progress count
					progress_count += 1;
				}
			}

			// Output final progress (it should be 100%)
			cout << "\b\b\b\b\b\b\b\b\b" << setw(8) << (int) ((double) progress_count / progress_total * 100) << "\%" << endl;
		}

		// Radial diffusion
		if (L_size >= 3) {
			progress_count = 0;
			progress_total = P_size * V_size * K_size;
			cout << "Radial diffusion  " << "         ";

			Matrix1D<double> PSD_L(L_size);
#pragma omp parallel for private(iP, iR, iV, iK, iL, PSD_L) shared(progress_total, progress_count) schedule(dynamic,1) collapse(3)
			for (iP = 0; iP < P_size; iP++) {
				for (iV = 0; iV < V_size; iV++) {
					for (iK = 0; iK < K_size; iK++) {

						if (omp_get_thread_num() == 0)
							cout << "\b\b\b\b\b\b\b\b\b" << setw(8)
									<< (int) ((double) progress_count / progress_total * 100) << "\%" << flush;

						// 1d slice
						PSD_L = PSD.wyzSlice(iP, iV, iK);

						// 1d diffusion
						Diffusion_1D(PSD_L, L.wyzSlice(iP, iV, iK), L_size, Ll_BC[iP][iV][iK], Lu_BC[iP][iV][iK],
								Ll_BC_type, Lu_BC_type, DLL.wyzSlice(iP, iV, iK), G_radial.wyzSlice(iP, iV, iK),
								Sources.wyzSlice(iP, iV, iK) * radial_losses,
								Losses.wyzSlice(iP, iV, iK) * radial_losses, dt);

						// copy results back
						for (iL = 0; iL < L_size; iL++)
							PSD[iP][iL][iV][iK] = PSD_L[iL];

#pragma omp critical // Progress output
						progress_count += 1;

					}
				}
			}
			cout << "\b\b\b\b\b\b\b\b\b" << setw(8) << (int) ((double) progress_count / progress_total * 100) << "\%"
					<< endl;
		}

		// Local diffusion
		if (V_size >= 3 && K_size >= 3) {
			int number_of_skipped_points = 0;
			progress_count = 0;
			progress_total = P_size * R_size;
			cout << "Local diffusion  " << "         ";

			Matrix2D<double> PSD_IK(V_size, K_size);
#pragma omp parallel for private(iP, iR, iV, iK, iL, PSD_IK) shared(progress_total, progress_count, number_of_skipped_points) schedule(dynamic,1) collapse(2)
			for (iP = 0; iP < P_size; iP++) {
				for (iR = 0; iR < R_size; iR++) {

					if (omp_get_thread_num() == 0)
						cout << "\b\b\b\b\b\b\b\b\b" << setw(8)
								<< (int) ((double) progress_count / progress_total * 100) << "\%" << flush;

#pragma omp critical // Progress count
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
					if (inversion_method == "Lapack") {
						Diffusion_2D(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
								Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR), // P, R, K
								Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR), // P, R, I
								Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
								DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
								Sources.wxSlice(iP, iR) * local_losses, Losses.wxSlice(iP, iR) * local_losses, dt);
					} else {
						Diffusion_2D_ADI3(PSD_IK, V.wxSlice(iP, iR), K.wxSlice(iP, iR), V_size, K_size,
								Vl_BC.xySlice(iP, iR), Vu_BC.xySlice(iP, iR), // P, R, K
								Kl_BC.xySlice(iP, iR), Ku_BC.xySlice(iP, iR), // P, R, I
								Vl_BC_type, Vu_BC_type, Kl_BC_type, Ku_BC_type, DVV.wxSlice(iP, iR),
								DKK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), DVK.wxSlice(iP, iR), G_local.wxSlice(iP, iR),
								Sources.wxSlice(iP, iR) * local_losses, Losses.wxSlice(iP, iR) * local_losses, dt);
					}

					// copy results back
					for (iV = 0; iV < V_size; iV++)
						for (iK = 0; iK < K_size; iK++)
							PSD[iP][iR][iV][iK] = PSD_IK[iV][iK];

				}
			}
			cout << "\b\b\b\b\b\b\b\b\b" << setw(8) << (int) ((double) progress_count / progress_total * 100) << "\%" << endl;
			cout << "Number of skipped points: " << (int) ((double) number_of_skipped_points/progress_total * 100) << "\%" << endl;

		}

		// Output
		if ((it % output_step) == 0) {
			time_string.str("");
			time_string << time;

			PSD_filename.str("");
			PSD_filename << outputFolder << "PSD_" << setw(5) << setfill('0') << int(it / output_step) << ".dat";
			cout << "Writing results: " << PSD_filename.str() << endl;
			PSD.writeToFile(PSD_filename.str(), time_string.str());
		}
	}

	cout << "Program was terminated correctly." << endl;
	cout << "Wall-clock time: " << (omp_get_wtime() - wall_timer) << " sec; ";
	cout << "CPU time: " << (double) ((clock() - start_time) / CLOCKS_PER_SEC) << " sec." << endl;

	return 0;
}
