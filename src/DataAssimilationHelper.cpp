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
    Matrix2D<double> result = result_1step;
    for (auto it = 1; it < nt; ++it) {
         result = result * result_1step;
    }

    return result;
}

double internal::conditionalMean(
    const Matrix1D<double>& a,
    const Matrix1D<bool>& flag) {
    if (a.size_q1 != flag.size_q1) {
        cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        cout << "Matrix sizes are different. ";
        cout << "a.size_q1 = " << a.size_q1 << ", ";
        cout << "flag.size_q1 = " << flag.size_q1 << endl;
        ;
        exit(EXIT_FAILURE);
    }

    double result{0.};
    int counter{0};
    for (auto i = 0; i < a.size_q1; ++i) {
        if (flag[i] && !isnan(a[i])) {
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
    int P_size = P.size_q1;
    int R_size = R.size_q2;

    auto dP = P[1][0] - P[0][0];
    auto dR = R[0][1] - R[0][0];

    vector<Matrix1D<bool>> P_flags(P_size);
    for (auto iP = 1; iP < P_size - 1; ++iP) {
        auto P_c = P[iP][0];
        P_flags[iP] = (data.P <= P_c + dP * 0.5) && (data.P >= P_c - dP * 0.5);
    }
    P_flags[0] = (data.P <= dP * 0.5) || (data.P >= P[P_size - 1][0] - dP * 0.5);
    P_flags[P_size - 1] = P_flags[0];

    // maybe? - replace boolean vectors with 1.0 - 0.0 vectors and perform dot product instead of conditional mean
    vector<Matrix1D<bool>> R_flags(R_size);
    for (auto iR = 0; iR < R_size; ++iR) {
        auto R_c = R[0][iR];
        R_flags[iR] = (data.R <= R_c + dR * 0.5) && (data.R >= R_c - dR * 0.5);
    }

    Matrix2D<double> result{P_size, R_size};
    if (type == "linear") {
        for (auto iP = 0; iP < P_size; ++iP) {
            const auto& P_flag = P_flags[iP];
            for (auto iR = 0; iR < R_size; ++iR) {
                const auto& R_flag = R_flags[iR];
                result[iP][iR] = conditionalMean(data.PSD, P_flag && R_flag);
            }
        }
    } else if (type == "log10") {
        for (auto iP = 0; iP < P_size; ++iP) {
            const auto& P_flag = P_flags[iP];
            for (auto iR = 0; iR < R_size; ++iR) {
                const auto& R_flag = R_flags[iR];
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
    int size{0};
    for (auto iP = 0; iP < data.size_q1; ++iP) {
        for (auto iR = 0; iR < data.size_q2; ++iR) {
            if (data[iP][iR] != FILLVAL && !isnan(data[iP][iR])) {
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
    for (auto iP = 0; iP < data.size_q1; ++iP) {
        for (auto iR = 0; iR < data.size_q2; ++iR) {
            int H_idx = iR + iP * data.size_q2;
            if (data[iP][iR] != FILLVAL && !isnan(data[iP][iR])) {
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
    std::cout << "xGEO ";
    auto xGEO = readProcessedMatFiles2D("xGEO", timeStart, timeEnd, parameters);
    std::cout << "InvMu ";
    auto invMu = readProcessedMatFiles3D("InvMu", timeStart, timeEnd, parameters);
    std::cout << "InvK ";
    result.K = readProcessedMatFiles2D("InvK", timeStart, timeEnd, parameters);
    std::cout << "PSD\n";
    result.PSD = readProcessedMatFiles3D("PSD", timeStart, timeEnd, parameters);

    if (!result.MLT.initialized || !xGEO.initialized || !invMu.initialized || !result.K.initialized || !result.PSD.initialized) {
        // return an empty mat file if any necessary file was missing
        return ProcessedMatFileData();
    }

    result.R.AllocateMemory(xGEO.size_q1);
    for (auto it = 0; it < xGEO.size_q1; ++it) {
        result.R[it] = std::sqrt(
            xGEO[it][0] * xGEO[it][0] +
            xGEO[it][1] * xGEO[it][1] +
            xGEO[it][2] * xGEO[it][2]);
    }
    result.V.AllocateMemory(invMu.size_q1, invMu.size_q2, invMu.size_q3);

    for (auto it = 0; it < invMu.size_q1; ++it) {
        for (auto ie = 0; ie < invMu.size_q2; ++ie) {
            for (auto ia = 0; ia < invMu.size_q3; ++ia) {
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

data_assimilation::ProcessedMatFileData internal::cat(const vector<ProcessedMatFileData>& pmfDataSplit) {
    if (pmfDataSplit.size() == 0) {
        return {};
    }

    ProcessedMatFileData result = pmfDataSplit[0];
    for (auto it = 1; it < pmfDataSplit.size(); ++it) {
        result.MLT = cat(result.MLT, pmfDataSplit[it].MLT);
        result.R = cat(result.R, pmfDataSplit[it].R);
        result.K = cat(result.K, pmfDataSplit[it].K);
        result.V = cat(result.V, pmfDataSplit[it].V);
        result.PSD = cat(result.PSD, pmfDataSplit[it].PSD);
    }
    return result;
}

Matrix1D<double> internal::interp1d_linear(
    const Matrix1D<double>& x_in,
    const Matrix1D<double>& f_in,
    const Matrix1D<double>& x_out,
    bool use_fillval) {
    if (!x_out.initialized)
        return x_out;

    if (x_in.size_q1 < 2) {
        std::cout << "In " << __func__ << ": at least 2 points are required for interpolation.\n";
        exit(EXIT_FAILURE);
    }

    if (x_in.size_q1 != f_in.size_q1) {
        std::cout << "In " << __func__ << ": input array sizes are not consistent.\n";
        exit(EXIT_FAILURE);
    }

    int sz = x_in.size_q1;
    int sz2 = x_out.size_q1;

    int notnan = 0;
    for (int i = 0; i < sz; ++i) {
        // check for FILLVAL?  && x_in[i] != FILLVAL && f_in[i]!= FILLVAL
        if (!isnan(x_in[i]) && !isnan(f_in[i])) {
            ++notnan;
        }
    }

    if (notnan <= 1) {
        Matrix1D<double> f_out(sz2);

        // f_out = FILLVAL;
        f_out = std::log10(-1.0);
        return f_out;
    } else {
        Matrix1D<double> x_in_new(notnan);
        Matrix1D<double> f_in_new(notnan);

        int ind = 0;
        for (int i = 0; i < sz; ++i) {
            if (!isnan(x_in[i]) && !isnan(f_in[i])) {
                x_in_new[ind] = x_in[i];
                f_in_new[ind] = f_in[i];
                ++ind;
            }
        }

        sz = notnan;
        double dx = x_in_new[1] - x_in_new[0];
        double x_t;
        double x1, x2, y1, y2;

        if (dx > 0) {
            for (int i = 1; i < sz - 1; ++i) {
                if (x_in_new[i + 1] - x_in_new[i] <= 0.) {
                    std::cout << "In " << __func__ << ": input grid is not strictly monotonic.\n";
                    exit(EXIT_FAILURE);
                }
            }

            Matrix1D<double> f_out(sz2);
            for (int j = 0; j < sz2; ++j) {
                x_t = x_out[j];

                if ((x_t > x_in_new[sz - 1]) || (x_t < x_in_new[0])) {
                    if (use_fillval)
                        f_out[j] = log10(-1.0);
                    else {
                        f_out[j] = (x_t > x_in_new[sz - 1]) ? f_in_new[sz - 1] : f_in_new[0];
                    }
                    continue;
                }

                for (int i = sz - 1; i >= 0; --i) {
                    if (x_t >= x_in_new[i]) {
                        x1 = x_in_new[i];
                        y1 = f_in_new[i];
                        break;
                    }
                }

                for (int i = 0; i < sz; ++i) {
                    if (x_t <= x_in_new[i]) {
                        x2 = x_in_new[i];
                        y2 = f_in_new[i];
                        break;
                    }
                }
                if (y1 == -1 / 0.0 || y2 == -1 / 0.0)
                    f_out[j] = y1;
                else if (x1 == x2)
                    f_out[j] = y1;
                else {
                    f_out[j] = ((y2 - y1) * x_t + (y1 * x2 - y2 * x1)) / (x2 - x1);
                }
            }
            return f_out;
        } else if (dx < 0) {
            for (int i = 1; i < sz - 1; ++i)
                if (x_in_new[i + 1] - x_in_new[i] >= 0.) {
                    std::cout << "In " << __func__ << ": input grid is not strictly monotonic.\n";
                    exit(EXIT_FAILURE);
                }

            Matrix1D<double> f_out(sz2);
            for (int j = 0; j < sz2; ++j) {
                x_t = x_out[j];

                if ((x_t < x_in_new[sz - 1]) || (x_t > x_in_new[0])) {
                    if (use_fillval)
                        f_out[j] = log10(-1.0);
                    else {
                        f_out[j] = (x_t < x_in_new[sz - 1]) ? f_in_new[sz - 1] : f_in_new[0];
                    }
                    continue;
                }

                for (int i = 0; i < sz; ++i) {
                    if (x_t >= x_in_new[i]) {
                        x2 = x_in_new[i];
                        y2 = f_in_new[i];
                        break;
                    }
                }

                for (int i = sz - 1; i >= 0; --i) {
                    if (x_t <= x_in_new[i]) {
                        x1 = x_in_new[i];
                        y1 = f_in_new[i];
                        break;
                    }
                }

                if (y1 == -1.0 / 0.0 || y2 == -1.0 / 0.0)
                    f_out[j] = y1;
                else if (x1 == x2)
                    f_out[j] = y1;
                else {
                    f_out[j] = ((y2 - y1) * x_t + (y1 * x2 - y2 * x1)) / (x2 - x1);
                }
            }

            return f_out;
        } else {
            std::cout << "In " << __func__ << ": input grid should not contain repeated elements.\n";
            exit(EXIT_FAILURE);
        }
    }
}

double internal::interp1d_linear(
    const Matrix1D<double>& x_in,
    const Matrix1D<double>& f_in,
    double x_out,
    bool use_badval) {
    Matrix1D<double> x_out_tmp(1), f_out_tmp(1);
    x_out_tmp[0] = x_out;

    f_out_tmp = interp1d_linear(x_in, f_in, x_out_tmp, use_badval);
    return f_out_tmp[0];
}

Matrix2D<double> internal::interp2d_linear_dependent(
    const Matrix2D<double>& V_in,
    const Matrix2D<double>& K_in,
    const Matrix2D<double>& PSD_in,
    const Matrix2D<double>& V_out,
    const Matrix2D<double>& K_out,
    bool use_badval) {
    if ((V_in.size_q1 != K_in.size_q1) || (K_in.size_q1 != PSD_in.size_q1) ||
        (V_in.size_q2 != K_in.size_q2) || (K_in.size_q2 != PSD_in.size_q2)) {
        std::cout << "In " << __func__ << ": input array sizes are not consistent." << std::endl;
        exit(EXIT_FAILURE);
    }

    if ((V_out.size_q1 != K_out.size_q1) || (V_out.size_q2 != K_out.size_q2)) {
        std::cout << "In " << __func__ << ": target array sizes are not consistent." << std::endl;
        exit(EXIT_FAILURE);
    }
    int nv_in = V_in.size_q1;
    int nk_in = K_in.size_q2;
    int nv_out = V_out.size_q1;
    int nk_out = K_out.size_q2;

    Matrix2D<double> PSD_out(nv_out, nk_out);

    Matrix3D<double> V_tmp(nv_in, nv_out, nk_out);
    Matrix3D<double> PSD_tmp(nv_in, nv_out, nk_out);

    for (int iV = 0; iV < nv_in; ++iV) {
        Matrix1D<double> K_1d_tmp = K_in.xSlice(iV);
        Matrix1D<double> V_1d_tmp = V_in.xSlice(iV);
        Matrix1D<double> PSD_1d_tmp = PSD_in.xSlice(iV);

        for (int iV2 = 0; iV2 < nv_out; ++iV2) {
            for (int iK2 = 0; iK2 < nk_out; ++iK2) {
                double K_1d_t_tmp = K_out[iV2][iK2];
                V_tmp[iV][iV2][iK2] = interp1d_linear(K_1d_tmp, V_1d_tmp, K_1d_t_tmp, use_badval);
                PSD_tmp[iV][iV2][iK2] = interp1d_linear(K_1d_tmp, PSD_1d_tmp, K_1d_t_tmp, use_badval);
            }
        }
    }

    for (int iV2 = 0; iV2 < nv_out; ++iV2) {
        for (int iK2 = 0; iK2 < nk_out; ++iK2) {
            Matrix1D<double> V_1d_tmp2 = V_tmp.yzSlice(iV2, iK2);
            Matrix1D<double> PSD_1d_tmp2 = PSD_tmp.yzSlice(iV2, iK2);

            PSD_out[iV2][iK2] = interp1d_linear(V_1d_tmp2, PSD_1d_tmp2, V_out[iV2][iK2], use_badval);
        }
    }
    // for (auto j = 0; j < nv_out; j++){
    //     for( auto k = 0; k < nk_out; k++){
    //             cout << PSD_out[j][k] << " ";
    //     }
    //     cout << endl;
    // }
    return PSD_out;
}

double internal::interp2d_linear_dependent(
    const Matrix2D<double>& V_in,
    const Matrix2D<double>& K_in,
    const Matrix2D<double>& PSD_in,
    double V_out,
    double K_out,
    bool use_badval) {
    Matrix2D<double> V_out_matrix{1, 1};
    Matrix2D<double> K_out_matrix{1, 1};
    V_out_matrix[0][0] = V_out;
    K_out_matrix[0][0] = K_out;
    auto result = interp2d_linear_dependent(V_in, K_in, PSD_in,
                                            V_out_matrix, K_out_matrix, use_badval);

    return result[0][0];
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
        if (y1 == -1.0 / 0.0 || y2 == -1.0 / 0.0 || x1 == x2) return y1;
        return ((y2 - y1) * x_out + (y1 * x2 - y2 * x1)) / (x2 - x1);
    };

    double minK = K_in.min();
    double maxK = K_in.max();
    if (K_out < minK || K_out > maxK) return std::log(-1.0);

    int nV = V_in.size_q1;
    int nK = K_in.size_q2;
    int non_nan_k = 0;

    for (int iK = 0; iK < nK; iK++)
        if (!isnan(K_in[0][iK])) {
            non_nan_k++;
        }
    if (non_nan_k < 2) {
        return std::log(-1);
    }
    Matrix1D<double> K_arr(non_nan_k);
    Matrix2D<double> V_arr(V_in.size_q1, non_nan_k);
    Matrix2D<double> PSD_arr(PSD_in.size_q1, non_nan_k);

    int counter = 0;
    for (int iK = 0; iK < nK; iK++) {
        if (!isnan(K_in[0][iK])) {
            K_arr[counter] = K_in[0][iK];
            for (int iV = 0; iV < nV; iV++) {
                V_arr[iV][counter] = V_in[iV][counter];
                PSD_arr[iV][counter] = PSD_in[iV][counter];
            }
            counter++;
        }
    }
    // double k_max = K_arr.max();
    // double k_min = K_arr.min();
    // if (K_out < k_min || K_out > k_max) return std::log(-1);

    // double dK = K_arr[1] - K_arr[0];
    int lower_k_neighbour = -1;
    int upper_k_neighbour = -1;
    double k_lower = -std::numeric_limits<double>::infinity();
    double k_upper = std::numeric_limits<double>::infinity();

    for (int iK = 0; iK < non_nan_k; iK++) {
        if (K_arr[iK] <= K_out && K_arr[iK] > k_lower) {
            k_lower = K_arr[iK];
            lower_k_neighbour = iK;
        }
        if (K_arr[iK] >= K_out && K_arr[iK] < k_upper) {
            k_upper = K_arr[iK];
            upper_k_neighbour = iK;
        }
    }
    if (lower_k_neighbour < 0 || upper_k_neighbour < 0) return std::log(-1);
    // if(dK > 0){
    //     for(int iK = 1; iK < non_nan_k; iK++){
    //         // if(K_arr[iK + 1] <= K_arr[iK]){
    //         //     exit(EXIT_FAILURE);
    //         // }
    //         if(K_arr[iK] >= K_out){
    //             lower_k_neighbour = iK - 1;
    //             upper_k_neighbour = iK;
    //             break;
    //         }
    //     }
    // }
    // else if(dK < 0){
    //     for(int iK = 1; iK < non_nan_k; iK++){
    //         // if(K_arr[iK + 1] >= K_arr[iK]){
    //         //     std::cout << "K_in must be monotonic!\n";
    //         //     for(int iK = 0; iK < nK - 1; iK++){
    //         //         std::cout << K_arr[iK] << ' ';
    //         //     }
    //         //     std::cout << std::endl;
    //         //     exit(EXIT_FAILURE);
    //         // }
    //         if(K_arr[iK] <= K_out){
    //             upper_k_neighbour = iK - 1;
    //             lower_k_neighbour = iK;
    //             break;
    //         }
    //     }
    // }
    // else{
    //     std::cout << "K_in must be monotonic!\n";
    //     exit(EXIT_FAILURE);
    // }

    /*
    Now we know the K-index of the lower and upper neighbour of the K_out grid point. For both indices
    filter nan's from the corresponding V and PSD arrays.
    */
    int non_nan_v_for_lower_k = 0;
    int non_nan_v_for_upper_k = 0;
    for (int iV = 0; iV < nV; iV++) {
        if (!isnan(V_arr[iV][lower_k_neighbour]) && !isnan(PSD_arr[iV][lower_k_neighbour])) non_nan_v_for_lower_k++;
        if (!isnan(V_arr[iV][upper_k_neighbour]) && !isnan(PSD_arr[iV][upper_k_neighbour])) non_nan_v_for_upper_k++;
    }
    if (non_nan_v_for_lower_k < 2 || non_nan_v_for_upper_k < 2) return std::log(-1);

    Matrix1D<double> V_for_lower_K(non_nan_v_for_lower_k);
    Matrix1D<double> V_for_upper_K(non_nan_v_for_upper_k);
    Matrix1D<double> PSD_for_lower_K(non_nan_v_for_lower_k);
    Matrix1D<double> PSD_for_upper_K(non_nan_v_for_upper_k);

    int counter_lower_k = 0;
    int counter_upper_k = 0;
    for (int iV = 0; iV < nV; iV++) {
        if (!isnan(V_arr[iV][lower_k_neighbour]) && !isnan(PSD_arr[iV][lower_k_neighbour])) {
            V_for_lower_K[counter_lower_k] = V_arr[iV][lower_k_neighbour];
            PSD_for_lower_K[counter_lower_k] = PSD_arr[iV][lower_k_neighbour];
            counter_lower_k++;
        }
        if (!isnan(V_arr[iV][upper_k_neighbour]) && !isnan(PSD_arr[iV][upper_k_neighbour])) {
            V_for_upper_K[counter_upper_k] = V_arr[iV][upper_k_neighbour];
            PSD_for_upper_K[counter_upper_k] = PSD_arr[iV][upper_k_neighbour];
            counter_upper_k++;
        }
    }

    int lower_v_neighbour = -1;
    int upper_v_neighbour = -1;
    double v_lower = -std::numeric_limits<double>::infinity();
    double v_upper = std::numeric_limits<double>::infinity();
    // double dV_lower = V_for_lower_K[1] - V_for_lower_K[0];
    // double dV_upper = V_for_upper_K[1] - V_for_upper_K[0];

    // double v_min = V_for_lower_K.min();
    // double v_max = V_for_lower_K.max();
    // if( V_out < v_min || V_out > v_max) return std::log(-1.0);
    // if(dV_lower > 0){
    //     for(int iV = 1; iV < non_nan_v_for_lower_k; iV++){
    //         // if(V_for_lower_K[iV] < V_for_lower_K[iV-1]){
    //         //     std::cout << "V_in must be monotonic!\n";
    //         //     exit(EXIT_FAILURE);
    //         // }
    //         if(V_for_lower_K[iV] >= V_out){
    //             lower_v_neighbour = iV - 1;
    //             upper_v_neighbour = iV;
    //             break;
    //         }
    //     }
    // }
    // else if(dV_lower < 0){
    //     for(int iV = 1; iV < non_nan_v_for_lower_k; iV++){
    //         // if(V_for_lower_K[iV] > V_for_lower_K[iV-1]){
    //         //     std::cout << "V_in must be monotonic!\n";
    //         //     exit(EXIT_FAILURE);
    //         // }
    //         if(V_for_lower_K[iV] <= V_out){
    //             upper_v_neighbour = iV - 1;
    //             lower_v_neighbour = iV;
    //             break;
    //         }
    //     }
    // }
    // else{
    //     std::cout << "V_in must be monotonic!\n";
    //     exit(EXIT_FAILURE);
    // }

    for (int iV = 0; iV < non_nan_v_for_lower_k; iV++) {
        if (V_for_lower_K[iV] <= V_out && V_for_lower_K[iV] > v_lower) {
            v_lower = V_for_lower_K[iV];
            lower_v_neighbour = iV;
        }
        if (V_for_lower_K[iV] >= V_out && V_for_lower_K[iV] < v_upper) {
            v_upper = V_for_lower_K[iV];
            upper_v_neighbour = iV;
        }
    }
    if (lower_v_neighbour < 0 || upper_v_neighbour < 0) return std::log(-1);

    // We know the V neighbours for the lower K value. Interpolate PSD for the lower K wrt V
    double x1 = V_for_lower_K[lower_v_neighbour];
    double x2 = V_for_lower_K[upper_v_neighbour];

    double y1 = PSD_for_lower_K[lower_v_neighbour];
    double y2 = PSD_for_lower_K[upper_v_neighbour];

    double PSD_lower = interpolation_lambda(x1, y1, x2, y2, V_out);
    // Now repeat for the upper K neighbour
    // v_min = V_for_upper_K.min();
    // v_max = V_for_upper_K.max();
    // if( V_out < v_min || V_out > v_max) return std::log10(-1.0);
    // if(dV_upper > 0){
    //     for(int iV = 1; iV < non_nan_v_for_upper_k; iV++){
    //         // if(V_for_upper_K[iV] <= V_for_upper_K[iV-1]){
    //         //     std::cout << "V_in must be monotonic!\n";
    //         //     exit(EXIT_FAILURE);
    //         // }
    //         if(V_for_upper_K[iV] >= V_out){
    //             lower_v_neighbour = iV - 1;
    //             upper_v_neighbour = iV;
    //             break;
    //         }
    //     }
    // }
    // else if(dV_upper < 0){
    //     for(int iV = 1; iV < non_nan_v_for_upper_k; iV++){
    //         // if(V_for_upper_K[iV] >= V_for_upper_K[iV-1]){
    //         //     std::cout << "V_in must be monotonic!\n";
    //         //     exit(EXIT_FAILURE);
    //         // }
    //         if(V_for_upper_K[iV] <= V_out){
    //             upper_v_neighbour = iV;
    //             lower_v_neighbour = iV-1;
    //             break;
    //         }
    //     }
    // }
    // else{
    //     std::cout << "V_in must be monotonic!\n";
    //     exit(EXIT_FAILURE);
    // }
    v_lower = -std::numeric_limits<double>::infinity();
    v_upper = std::numeric_limits<double>::infinity();
    lower_v_neighbour = -1;
    upper_v_neighbour = -1;

    for (int iV = 0; iV < non_nan_v_for_upper_k; iV++) {
        if (V_for_upper_K[iV] <= V_out && V_for_upper_K[iV] > v_lower) {
            v_lower = V_for_upper_K[iV];
            lower_v_neighbour = iV;
        }
        if (V_for_upper_K[iV] >= V_out && V_for_upper_K[iV] < v_upper) {
            v_upper = V_for_upper_K[iV];
            upper_v_neighbour = iV;
        }
    }
    if (lower_v_neighbour < 0 || upper_v_neighbour < 0) return std::log(-1);

    x1 = V_for_upper_K[lower_v_neighbour];
    x2 = V_for_upper_K[upper_v_neighbour];
    y1 = PSD_for_upper_K[lower_v_neighbour];
    y2 = PSD_for_upper_K[upper_v_neighbour];

    double PSD_upper = interpolation_lambda(x1, y1, x2, y2, V_out);

    // Knowing PSD at V_out for the lower and upper K neighbours, we can now do one final
    // interpolation between them wrt K

    x1 = K_arr[lower_k_neighbour];
    x2 = K_arr[upper_k_neighbour];
    y1 = PSD_lower;
    y2 = PSD_upper;

    return interpolation_lambda(x1, y1, x2, y2, K_out);
}
vector<vector<data_assimilation::Observations>> internal::interpolate_old(
    const ProcessedMatFileData& data,
    const Matrix2D<double>& V_grid,
    const Matrix2D<double>& K_grid) {
    int nT = data.MLT.size_q1;

    Matrix1D<double> P{nT};
    for (auto it = 0; it < nT; ++it) {
        P[it] = fmod(data.MLT[it] + 12., 24.) * M_PI / 12.;
    }

    int V_size = V_grid.size_q1;
    int K_size = K_grid.size_q2;

    Matrix3D<double> K_in(nT, data.V.size_q2, data.K.size_q2);
    for (auto it = 0; it < nT; ++it) {
        for (auto iV = 0; iV < K_in.size_q2; ++iV) {
            for (auto iK = 0; iK < K_in.size_q3; ++iK) {
                K_in[it][iV][iK] = data.K[it][iK];
            }
        }
    }

    vector<vector<Observations>> result(V_size, vector<Observations>(K_size));
    // int counter = 0;
#pragma omp parallel for schedule(dynamic, 1) collapse(2)
    for (int iV = 0; iV < V_size; ++iV) {
        for (int iK = 0; iK < K_size; ++iK) {
            Observations& obs = result[iV][iK];
            obs.R = data.R;
            obs.P = P;
            obs.PSD.AllocateMemory(nT);
            for (auto it = 0; it < nT; ++it) {
                auto PSD_in = data.PSD.xSlice(it);
                // Matrix2D<double> K_in (data.V.size_q2, data.K.size_q2);
                // for (auto iV2 = 0; iV2 < K_in.size_q1; ++iV2) {
                //     for (auto iK2 = 0; iK2 < K_in.size_q2; ++iK2){
                //         K_in[iV2][iK2] = data.K[it][iK2];
                //     }
                // }
                // obs.PSD[it] = pow(10., interp2d_linear_dependent(
                //     log10(data.V.xSlice(it)), K_in, log10(PSD_in),
                //     V_grid[iV][iK], K_grid[iV][iK], false));

                obs.PSD[it] = pow(10., interp2d_linear_dependent(
                                           log10(data.V.xSlice(it)), K_in.xSlice(it), log10(PSD_in),
                                           log10(V_grid[iV][iK]), K_grid[iV][iK], true));
                // obs.PSD[it] = interp2d_linear_dependent(
                //     data.V.xSlice(it), K_in.xSlice(it), PSD_in,
                //     V_grid[iV][iK], K_grid[iV][iK], false);
            }
        }
    }
    return result;
}

vector<vector<data_assimilation::Observations>> internal::interpolate(
    const std::vector<ProcessedMatFileData>& data,
    const Matrix2D<double>& V_grid,
    const Matrix2D<double>& K_grid) {
    int V_size = V_grid.size_q1;
    int K_size = K_grid.size_q2;

    int nT = 0;
    for (auto& instrumentData : data) nT += instrumentData.MLT.size_q1;
    Matrix1D<double> R(nT);
    Matrix1D<double> P(nT);
    vector<Matrix3D<double>> K_in;
    K_in.reserve(data.size());
    int counter = 0;
    for (auto& instrumentData : data) {
        K_in.emplace_back(instrumentData.MLT.size_q1, instrumentData.V.size_q2, instrumentData.K.size_q2);
        for (auto it = 0; it < instrumentData.MLT.size_q1; ++it) {
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
    vector<vector<Observations>> result(V_size, vector<Observations>(K_size));
    // int counter = 0;
#pragma omp parallel for schedule(dynamic, 1) collapse(2)
    for (int iV = 0; iV < V_size; ++iV) {
        for (int iK = 0; iK < K_size; ++iK) {
            Observations& obs = result[iV][iK];
            obs.R = R;
            obs.P = P;
            obs.PSD.AllocateMemory(nT);
            int counter = 0;
            for (int instrumentIndex = 0; instrumentIndex < data.size(); instrumentIndex++) {
                auto& instrumentData = data[instrumentIndex];
                for (auto it = 0; it < instrumentData.MLT.size_q1; ++it) {
                    auto PSD_in = instrumentData.PSD.xSlice(it);
                    // obs.PSD[counter++] = pow(10., interp2d_linear_dependent(
                    //     log10(instrumentData.V.xSlice(it)),
                    //     K_in[instrumentIndex].xSlice(it),
                    //     log10(PSD_in),
                    //     log10(V_grid[iV][iK]),
                    //     K_grid[iV][iK], true));
                    obs.PSD[counter++] = pow(10., interp2d_four_corners(
                                                      log10(instrumentData.V.xSlice(it)),
                                                      K_in[instrumentIndex].xSlice(it),
                                                      log10(PSD_in),
                                                      log10(V_grid[iV][iK]),
                                                      K_grid[iV][iK], true));
                    // obs.PSD[it] = interp2d_linear_dependent(
                    //     data.V.xSlice(it), K_in.xSlice(it), PSD_in,
                    //     V_grid[iV][iK], K_grid[iV][iK], false);
                }
            }
        }
    }
    return result;
}