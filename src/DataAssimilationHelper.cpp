#include "DataAssimilationHelper.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "CustomDate.h"
#include "DataAssimilation.h"
#include "MatrixOperations.h"
#include "PMF.h"

using std::string;

extern const double FILLVAL{-1e31};
namespace internal = data_assimilation::internal;

Matrix2D<double> internal::getModelMatrixConvection2DNoStabilityCheck(
    const Matrix2D<double>& VP,
    const Matrix2D<double>& VR,
    const Matrix2D<double>& Loss,
    double timeStep,
    double dP,
    double dR) {
    int P_size = VP.size_q1;
    int R_size = VP.size_q2;

    Matrix2D<double> result(P_size * R_size, P_size * R_size);
    result = 0.;

    for (auto iP = 0; iP < P_size; ++iP) {
        for (auto iR = 0; iR < R_size; ++iR) {
            auto iP_c = iP;
            size_t iP_l, iP_r;
            if (iP == 0 || iP == P_size - 1) {
                iP_l = P_size - 2;
                iP_r = 1;
            } else {
                iP_l = iP - 1;
                iP_r = iP + 1;
            }

            auto iR_c = iR;
            auto iR_l = iR - 1;
            auto iR_r = iR + 1;

            auto VP_c = VP[iP_c][iR_c];
            auto VR_c = VR[iP_c][iR_c];

            double m_cc = 1. - std::abs(VP_c) * timeStep / dP - std::abs(VR_c) * timeStep / dR;
            //auto Loss_c = Loss[iP_c][iR_c] * 0;
            // + timeStep * Loss_c;

            double m_lc, m_rc, m_cl, m_cr;
            if (VP_c > 0.) {
                m_lc = VP_c * timeStep / dP;
                m_rc = 0.;
            } else {
                m_lc = 0.;
                m_rc = -VP_c * timeStep / dP;
            }

            if (VR_c > 0.) {
                m_cl = VR_c * timeStep / dR;
                m_cr = 0.;
            } else {
                m_cl = 0.;
                m_cr = -VR_c * timeStep / dR;
            }

            if (m_cc < 0 || m_lc < 0 || m_rc < 0 || m_cl < 0 || m_cr < 0) {
                std::cout << "This model matrix is unstable. Probably due to a misscalculation in the CLF condition." << std::endl;
                exit(EXIT_FAILURE);
            }

            if (iR != 0 && iR != R_size - 1) {
                result[iR_c + iP_c * R_size][iR_c + iP_c * R_size] = m_cc;
                result[iR_c + iP_c * R_size][iR_c + iP_l * R_size] = m_lc;
                result[iR_c + iP_c * R_size][iR_c + iP_r * R_size] = m_rc;
                result[iR_c + iP_c * R_size][iR_l + iP_c * R_size] = m_cl;
                result[iR_c + iP_c * R_size][iR_r + iP_c * R_size] = m_cr;
            } else {
                result[iR_c + iP_c * R_size][iR_c + iP_c * R_size] = 1.;
            }
        }
    }

    return result;
}

std::pair<int, double> internal::splitTimeStepCourantCondition(
    double maximumCourantNumber,
    double timeStep,
    const Matrix2D<double>& VP,
    const Matrix2D<double>& VR,
    double dP,
    double dR) {
    auto compareAbs = [](double a, double b) {
        return std::abs(a) < std::abs(b);
    };

    Matrix2D<double> c_p = VP * timeStep / dP;
    Matrix2D<double> c_r = VR * timeStep / dR;
    Matrix2D<double> c_pr = c_p.abs() + c_r.abs();

    double cour_num = c_pr.maxabs();

    int nt = (cour_num <= maximumCourantNumber)
                 ? 1
                 : std::ceil(cour_num / maximumCourantNumber);

    double timeStep_new = timeStep / static_cast<double>(nt);

    return {nt, timeStep_new};
}

Matrix2D<double> internal::getModelMatrixConvection2D(
    const Matrix2D<double>& VP,
    const Matrix2D<double>& VR,
    const Matrix2D<double>& Loss,
    double timeStep,
    double dP,
    double dR) {
    const double maximumCourantNumber = 1.0;
    auto stepSplit = splitTimeStepCourantCondition(
        maximumCourantNumber, timeStep,
        VP, VR, dP, dR);
    const int& nt = stepSplit.first;
    const double& timeStep_new = stepSplit.second;

    Matrix2D<double> result_1step = getModelMatrixConvection2DNoStabilityCheck(
        VP, VR, Loss, timeStep_new, dP, dR);

    // Matrix2D<double> result = result_1step;
    // for (auto it = 1; it < nt; ++it) {
    //      result = result * result_1step;
    // }
    Matrix2D<double> result = result_1step^nt;

    return result;
}

double internal::conditionalMean(
    const Matrix1D<double>& a,
    const Matrix1D<bool>& flag) {
    if (a.size_q1 != flag.size_q1) {
        std::cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        std::cout << "Matrix sizes are different. ";
        std::cout << "a.size_q1 = " << a.size_q1 << ", ";
        std::cout << "flag.size_q1 = " << flag.size_q1 << std::endl;
        ;
        exit(EXIT_FAILURE);
    }

    double result{0.};
    int counter{0};
    for (auto i = 0; i < a.size_q1; ++i) {
        if (flag[i] && !std::isnan(a[i])) {
            result += a[i];
            ++counter;
        }
    }

    if (!counter) {
        // result = FILLVAL;
        result = NAN;  // std::nan(NULL);
    } else {
        result /= counter;
    }

    return result;
}

Matrix2D<double> internal::bin(
    const data_assimilation::Observations& data,
    const Matrix2D<double>& P,
    const Matrix2D<double>& R,
    const string& type) {
    size_t P_size = P.size_q1;
    size_t R_size = R.size_q2;

    double dP = P[1][0] - P[0][0];
    double dR = R[0][1] - R[0][0];

    std::vector<Matrix1D<bool>> P_flags(P_size);
    for (size_t iP = 1; iP < P_size - 1; ++iP) {
        double P_c = P[iP][0];
        P_flags[iP] = (data.P <= P_c + dP * 0.5) && (data.P >= P_c - dP * 0.5);
    }
    P_flags[0] = (data.P <= dP * 0.5) || (data.P >= P[P_size - 1][0] - dP * 0.5);
    P_flags[P_size - 1] = P_flags[0];

    // maybe? - replace boolean vectors with 1.0 - 0.0 vectors and perform dot product instead of conditional mean
    std::vector<Matrix1D<bool>> R_flags(R_size);
    for (size_t iR = 0; iR < R_size; ++iR) {
        double R_c = R[0][iR];
        R_flags[iR] = (data.R <= R_c + dR * 0.5) && (data.R >= R_c - dR * 0.5);
    }

    Matrix2D<double> result{P_size, R_size};
    if (type == "linear") {
        for (size_t iP = 0; iP < P_size; ++iP) {
            const Matrix1D<bool>& P_flag = P_flags[iP];
            for (size_t iR = 0; iR < R_size; ++iR) {
                const Matrix1D<bool>& R_flag = R_flags[iR];
                result[iP][iR] = conditionalMean(data.PSD, P_flag && R_flag);
            }
        }
    } else if (type == "log10") {
        for (size_t iP = 0; iP < P_size; ++iP) {
            const Matrix1D<bool>& P_flag = P_flags[iP];
            for (size_t iR = 0; iR < R_size; ++iR) {
                const Matrix1D<bool>& R_flag = R_flags[iR];
                result[iP][iR] = std::pow(10.,
                                          conditionalMean(log10(data.PSD), P_flag && R_flag));
            }
        }
    } else {
        std::cout << "Unknown binning type: " << type << std::endl;
        exit(EXIT_FAILURE);
    }
    return result;
}

data_assimilation::ObservationSpace internal::convertToObservationSpace(
    const Matrix2D<double>& data) {
    size_t size{0};
    for (size_t iP = 0; iP < data.size_q1; ++iP) {
        for (size_t iR = 0; iR < data.size_q2; ++iR) {
            if (data[iP][iR] != FILLVAL && !std::isnan(data[iP][iR])) {
                ++size;
            }
        }
    }

    if (size == 0) {
        ObservationSpace result{
            Matrix1D<double>{1},
            Matrix2D<double>{1, data.size_q1 * data.size_q2}};
        result.H = 0.0;
        result.data = 0.0;

        return result;
    }

    ObservationSpace result{
        Matrix1D<double>{size},
        Matrix2D<double>{size, data.size_q1 * data.size_q2}};

    result.H = 0.;

    int counter{0};
    for (size_t iP = 0; iP < data.size_q1; ++iP) {
        for (size_t iR = 0; iR < data.size_q2; ++iR) {
            int H_idx = iR + iP * data.size_q2;
            if (data[iP][iR] != FILLVAL && !std::isnan(data[iP][iR])) {
                result.data[counter] = data[iP][iR];
                result.H[counter][H_idx] = 1.0;
                ++counter;
            }
        }
    }
    return result;
}

bool internal::str2bool(const string& str) {
    return (str == "true" || str == "True" || str == "TRUE");
}

void internal::getParameterBool(
    ParametersIni& parameters,
    const string& name,
    bool& parameter,
    bool required) {
    string parameterString;
    parameters.getParameter(name, parameterString, required);
    parameter = str2bool(parameterString);
}

data_assimilation::ProcessedMatFileData internal::readData(
    double timeStart,
    double timeEnd,
    const pmf::Parameters& parameters) {
    std::cout << "Reading data from " << parameters.satellite << " " << parameters.instrument << "...\n";
    ProcessedMatFileData result;
    std::cout << "\tMLT ";
    result.MLT = readProcessedMatFiles1D("MLT", timeStart, timeEnd, parameters);
    // std::cout << "xGEO ";
    // auto xGEO = readProcessedMatFiles2D("xGEO", timeStart, timeEnd, parameters);
    std::cout << "R0 ";
    result.R = readProcessedMatFiles1D("R0", timeStart, timeEnd, parameters);
    std::cout << "InvMu ";
    auto invMu = readProcessedMatFiles3D("InvMu", timeStart, timeEnd, parameters);
    std::cout << "InvK ";
    result.K = readProcessedMatFiles2D("InvK", timeStart, timeEnd, parameters);
    std::cout << "PSD\n";
    result.PSD = readProcessedMatFiles3D("PSD", timeStart, timeEnd, parameters);

    if (!result.MLT.initialized || !result.R.initialized || !invMu.initialized || !result.K.initialized || !result.PSD.initialized) {
        // return an empty mat file if any necessary file was missing
        return ProcessedMatFileData();
    }

    // result.R.AllocateMemory(xGEO.size_q1);
    // for (auto it = 0; it < xGEO.size_q1; ++it) {
    //     result.R[it] = std::sqrt(
    //         xGEO[it][0] * xGEO[it][0] +
    //         xGEO[it][1] * xGEO[it][1] +
    //         xGEO[it][2] * xGEO[it][2]);
    // }
    result.V.AllocateMemory(invMu.size_q1, invMu.size_q2, invMu.size_q3);

    for (size_t it = 0; it < invMu.size_q1; ++it) {
        for (size_t ie = 0; ie < invMu.size_q2; ++ie) {
            for (size_t ia = 0; ia < invMu.size_q3; ++ia) {
                result.V[it][ie][ia] = invMu[it][ie][ia] * std::pow(result.K[it][ia] + 0.5, 2.);
                // replace all zeros with nan's
                if (result.PSD[it][ie][ia] < 1e-21) result.PSD[it][ie][ia] = std::log10(-1.0);

                // unit conversion
                result.PSD[it][ie][ia] *= 3e7;
            }
        }
    }
    // unit conversion
    return result;
}

data_assimilation::ProcessedMatFileData internal::cat(const std::vector<ProcessedMatFileData>& pmfDataSplit) {
    if (pmfDataSplit.size() == 0) {
        return {};
    }

    ProcessedMatFileData result = pmfDataSplit[0];
    for (size_t it = 1; it < pmfDataSplit.size(); ++it) {
        result.MLT = cat(result.MLT, pmfDataSplit[it].MLT);
        result.R = cat(result.R, pmfDataSplit[it].R);
        result.K = cat(result.K, pmfDataSplit[it].K);
        result.V = cat(result.V, pmfDataSplit[it].V);
        result.PSD = cat(result.PSD, pmfDataSplit[it].PSD);
    }
    return result;
}

double internal::interp2d_four_corners(
    const Matrix2D<double>& V_in, const Matrix2D<double>& K_in,
    const Matrix2D<double>& PSD_in, double V_out, double K_out,
    bool use_badval) {
    /*
    We assume K_in to vary only in the the second dimension. Find the two K values in K_in that neighbour
    the requested K_out value. For each of those two K values find the two V values that neighbour the
    requested V_out value.
    Perform linear interpolation of PSD wrt V for each of the two K values.
    Finally interpolate the V-interpolated PSD wrt K.

    Return Nan if there are less than two non-nan K values, less than two non-nan V values for the
    two K neighbours or if either K_out or V_out is outside of K_in or V_in respectively
    */

    // linear interpolation as a lambda
    // as we're interpolating log values we can sometimes expect -infinity
    auto interpolation_lambda = [](double x1, double y1, double x2, double y2, double x_out) {
        // if (y1 == -1.0 / 0.0 || y2 == -1.0 / 0.0 || x1 == x2) return y1;
        // Explisitly define -infinity, otherwise it results in the compilation error
        if (y1 == -std::numeric_limits<double>::infinity() || y2 == -std::numeric_limits<double>::infinity() || x1 == x2) return y1;                
        return ((y2 - y1) * x_out + (y1 * x2 - y2 * x1)) / (x2 - x1);
    };

    double minK = K_in.min();
    double maxK = K_in.max();
    if (K_out < minK || K_out > maxK) return std::log(-1.0);

    size_t nV = V_in.size_q1;
    size_t nK = K_in.size_q2;

    int lower_k_neighbour = -1;
    int upper_k_neighbour = -1;
    double k_lower = -std::numeric_limits<double>::infinity();
    double k_upper = std::numeric_limits<double>::infinity();

    for (size_t iK = 0; iK < nK; iK++) 
    {
        // we check all V channels altough K is independent of V
        // in case of missing data in the first V channel K_in[0][iK]
        for(size_t iV = 0; iV < nV; iV++)
        {
            if (K_in[iV][iK] <= K_out && K_in[iV][iK] > k_lower) {
                k_lower = K_in[iV][iK];
                lower_k_neighbour = iK;
                break;
            }
        }
        for(size_t iV = 0; iV < nV; iV++)
        {
            if (K_in[iV][iK] >= K_out && K_in[iV][iK] < k_upper) {
                k_upper = K_in[iV][iK];
                upper_k_neighbour = iK;
                break;
            }
        }
    }
    if (lower_k_neighbour < 0 || upper_k_neighbour < 0) return std::log(-1);

    /*
    Now we know the K-index of the lower and upper neighbour of the K_out grid point. For both indices
    check if there are at least two non-nans in the corresponding V and PSD arrays.
    */
    int non_nan_v_for_lower_k = 0;
    int non_nan_v_for_upper_k = 0;
    for (size_t iV = 0; iV < nV; iV++) {
        if (!std::isnan(V_in[iV][lower_k_neighbour]) && !std::isnan(PSD_in[iV][lower_k_neighbour])) non_nan_v_for_lower_k++;
        if (!std::isnan(V_in[iV][upper_k_neighbour]) && !std::isnan(PSD_in[iV][upper_k_neighbour])) non_nan_v_for_upper_k++;
    }
    if (non_nan_v_for_lower_k < 2 || non_nan_v_for_upper_k < 2) return std::log(-1);

    int lower_v_neighbour = -1;
    int upper_v_neighbour = -1;
    double v_lower = -std::numeric_limits<double>::infinity();
    double v_upper = std::numeric_limits<double>::infinity();
  
    for (int iV = 0; iV < nV; iV++) {
        if (V_in[iV][lower_k_neighbour] <= V_out && V_in[iV][lower_k_neighbour] > v_lower) {
            v_lower = V_in[iV][lower_k_neighbour];
            lower_v_neighbour = iV;
        }
        if (V_in[iV][lower_k_neighbour] >= V_out && V_in[iV][lower_k_neighbour] < v_upper) {
            v_upper = V_in[iV][lower_k_neighbour];
            upper_v_neighbour = iV;
        }
    }
    if (lower_v_neighbour < 0 || upper_v_neighbour < 0) return std::log(-1);

    // We know the V neighbours for the lower K value. Interpolate PSD for the lower K wrt V
    double y1 = PSD_in[lower_v_neighbour][lower_k_neighbour];
    double y2 = PSD_in[upper_v_neighbour][lower_k_neighbour];

    double PSD_lower = interpolation_lambda(v_lower, y1, v_upper, y2, V_out);
  
    // Now repeat for the upper K neighbour
    v_lower = -std::numeric_limits<double>::infinity();
    v_upper = std::numeric_limits<double>::infinity();
    lower_v_neighbour = -1;
    upper_v_neighbour = -1;

    for (size_t iV = 0; iV < nV; iV++) {
        if (V_in[iV][upper_k_neighbour] <= V_out && V_in[iV][upper_k_neighbour] > v_lower) {
            v_lower = V_in[iV][upper_k_neighbour];
            lower_v_neighbour = iV;
        }
        if (V_in[iV][upper_k_neighbour] >= V_out && V_in[iV][upper_k_neighbour] < v_upper) {
            v_upper = V_in[iV][upper_k_neighbour];
            upper_v_neighbour = iV;
        }
    }
    if (lower_v_neighbour < 0 || upper_v_neighbour < 0) return std::log(-1);

    y1 = PSD_in[lower_v_neighbour][upper_k_neighbour];
    y2 = PSD_in[upper_v_neighbour][upper_k_neighbour];

    double PSD_upper = interpolation_lambda(v_lower, y1, v_upper, y2, V_out);

    // Knowing PSD at V_out for the lower and upper K neighbours, we can now do one final
    // interpolation between them wrt K

    y1 = PSD_lower;
    y2 = PSD_upper;

    return interpolation_lambda(k_lower, y1, k_upper, y2, K_out);
}

std::vector<std::vector<data_assimilation::Observations>> internal::interpolate(
    const std::vector<ProcessedMatFileData>& data,
    const Matrix2D<double>& V_grid,
    const Matrix2D<double>& K_grid) {
    size_t V_size = V_grid.size_q1;
    size_t K_size = K_grid.size_q2;

    int nT = 0;
    for (auto& instrumentData : data) nT += instrumentData.MLT.size_q1;
    Matrix1D<double> R(nT);
    Matrix1D<double> P(nT);
    std::vector<Matrix3D<double>> K_in;
    K_in.reserve(data.size());
    int counter = 0;
    for (auto& instrumentData : data) {
        K_in.emplace_back(instrumentData.MLT.size_q1, instrumentData.V.size_q2, instrumentData.K.size_q2);
        for (size_t it = 0; it < instrumentData.MLT.size_q1; ++it) {
            R[counter] = instrumentData.R[it];
            P[counter] = fmod(instrumentData.MLT[it] + 12., 24.) * M_PI / 12.;
            counter++;
            for (auto iV = 0; iV < instrumentData.V.size_q2; ++iV) {
                for (auto iK = 0; iK < instrumentData.K.size_q2; ++iK) {
                    K_in.back()[it][iV][iK] = instrumentData.K[it][iK];
                }
            }
        }
    }
    std::vector<std::vector<Observations>> result(V_size, std::vector<Observations>(K_size));
    // int counter = 0;
#pragma omp parallel for schedule(dynamic, 1) collapse(2)
    for (auto iV = 0; iV < V_size; ++iV) {
        for (auto iK = 0; iK < K_size; ++iK) {
            Observations& obs = result[iV][iK];
            obs.R = R;
            obs.P = P;
            obs.PSD.AllocateMemory(nT);
            int counter = 0;
            for (size_t instrumentIndex = 0; instrumentIndex < data.size(); instrumentIndex++) {
                auto& instrumentData = data[instrumentIndex];
                for (auto it = 0; it < instrumentData.MLT.size_q1; ++it) {
                    auto PSD_in = instrumentData.PSD.xSlice(it);
                    obs.PSD[counter++] = pow(10., interp2d_four_corners(
                        log10(instrumentData.V.xSlice(it)),
                        K_in[instrumentIndex].xSlice(it),
                        log10(PSD_in),
                        log10(V_grid[iV][iK]),
                        K_grid[iV][iK], true
                    ));
                }
            }
        }
    }
    return result;
}