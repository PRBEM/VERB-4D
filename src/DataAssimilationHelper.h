#pragma once
#include "DataAssimilation.h"
#include "Matrix.h"
#include "PMF.h"
#include "Parameters.h"
#include <string>
#include <vector>

extern const double FILLVAL;
using ParametersIni=::Parameters;
namespace data_assimilation::internal {
    Matrix2D<double> getModelMatrixConvection2D(
        const Matrix2D<double>& VP,
        const Matrix2D<double>& VR,
        const Matrix2D<double>& Loss,
        double timeStep,
        double dP,
        double dR
    );
    Matrix2D<double> getModelMatrixConvection2DNoStabilityCheck(
        const Matrix2D<double>& VP,
        const Matrix2D<double>& VR,
        const Matrix2D<double>& Loss,
        double timeStep,
        double dP,
        double dR
    );
    std::pair<int, double> splitTimeStepCourantCondition(
        double maximumCourantNumber, double timeStep, 
        const Matrix2D<double>& VP, const Matrix2D<double>& VR, 
        double dP, double dR
    );
    Matrix2D<double> getModelMatrixConvection2D(
        const Matrix2D<double>& VP,
        const Matrix2D<double>& VR,
        const Matrix2D<double>& Loss,
        double timeStep,
        double dP,
        double dR
    );

    double conditionalMean(const Matrix1D<double>& a, const Matrix1D<bool>& flag);
    
    Matrix2D<double> bin(
        const data_assimilation::Observations& data, 
        const Matrix2D<double>& P, const Matrix2D<double>& R,
        const std::string& type = "linear"
    );

    data_assimilation::ObservationSpace convertToObservationSpace(const Matrix2D<double>& data);

    void getParameterBool(ParametersIni& parameters, const std::string& name, bool& parameter, bool required);

    data_assimilation::ProcessedMatFileData readData(double timeStart, double timeEnd, const pmf::Parameters& parameters);

    data_assimilation::ProcessedMatFileData cat(const std::vector<ProcessedMatFileData>& pmfDataSplit);
    
    Matrix1D<double> interp1d_linear(
        const Matrix1D<double>& x_in, 
        const Matrix1D<double>& f_in,
        const Matrix1D<double>& x_out, bool use_fillval
    );
    double interp1d_linear(
        const Matrix1D<double>& x_in, 
        const Matrix1D<double>& f_in,
        double x_out, bool use_badval
    );
    Matrix2D<double> interp2d_linear_dependent(
        const Matrix2D<double>& V_in, const Matrix2D<double>& K_in, 
        const Matrix2D<double>& PSD_in,const Matrix2D<double>& V_out,
        const Matrix2D<double>& K_out, bool use_badval
    );
    double interp2d_linear_dependent(
        const Matrix2D<double>& V_in, const Matrix2D<double>& K_in,
        const Matrix2D<double>& PSD_in, double V_out, double K_out,
        bool use_badval
    );
   std::vector<std::vector<data_assimilation::Observations>> interpolate_old(
        const data_assimilation::ProcessedMatFileData& data,
        const Matrix2D<double>& V_grid,
        const Matrix2D<double>& K_grid 
    );

   std::vector<std::vector<data_assimilation::Observations>> interpolate(
        const std::vector<data_assimilation::ProcessedMatFileData>& data,
        const Matrix2D<double>& V_grid,
        const Matrix2D<double>& K_grid 
    );
    bool str2bool(const std::string& str);
}