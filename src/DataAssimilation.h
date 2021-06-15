#pragma once

#include "Matrix.h"
#include "PMF.h"
#include "Parameters.h"
#include <string>
#include <vector>
using ParametersIni=::Parameters;
namespace data_assimilation {  
    struct Parameters {
        bool runDataAssimilation;
        bool useLog;
        double timeStep;
        double modelError;
        double observationError;
    };
    struct Observations {
        Matrix1D<double> P;
        Matrix1D<double> R;
        Matrix1D<double> PSD;
    };
    struct ObservationSpace {
        Matrix1D<double> data;
        Matrix2D<double> H;
    };
    struct ProcessedMatFileData {
        Matrix1D<double> MLT;
        Matrix1D<double> R;
        Matrix3D<double> V;
        Matrix2D<double> K;
        Matrix3D<double> PSD; 
    };
    class Convection2DAnalysisCovariances {
        public:
            Convection2DAnalysisCovariances(size_t P_size, size_t R_size, size_t V_size, size_t K_size);
            const std::vector<Matrix2D<double>>& operator[](size_t iV) const;
            std::vector<Matrix2D<double>>& operator[](size_t iV);
        private:
            std::vector<std::vector<Matrix2D<double>>> _data;
    };
    class DataAssimilationManagerConvection {
        public:
            DataAssimilationManagerConvection(const std::string& parametersFile, 
                double time_start, double time_end, const Matrix2D<double>& V, 
                const Matrix2D<double>& K, int P_size, int R_size);

            void assimilate(double time, Matrix4D<double>& PSD, 
                const Matrix4D<double>& P, const Matrix4D<double>& R,
                const Matrix4D<double>& VP, const Matrix4D<double>& VR,
                const Matrix4D<double>& Loss, 
                double dt);

        private:
            double _timeStart;
            double _timeEnd;
            Matrix2D<double> _V;
            Matrix2D<double> _K;
            Convection2DAnalysisCovariances _analysisCovarianceConvection;

            bool _runDataAssimilation;
            double _timePrev;
            double _timeNext;
            Parameters _assimilationParameters;
            std::vector<pmf::Parameters> _dataParameters;
    };
    
    Parameters readParameters(const std::string& filename);

    std::vector<std::vector<Observations>> getObservations_old(
        double timeStart,
        double timeEnd,
        const Matrix2D<double>& V,
        const Matrix2D<double>& K,
        const std::vector<pmf::Parameters>& parameters
    );

    std::vector<std::vector<Observations>> getObservations(
        double timeStart,
        double timeEnd,
        const Matrix2D<double>& V,
        const Matrix2D<double>& K,
        const std::vector<pmf::Parameters>& parameters
    );


    void runKalmanFilterConvection2D (
        Matrix2D<double>& forecast,
        Matrix2D<double>& analysisErrorCovariance,
        const Matrix2D<double>& P,
        const Matrix2D<double>& R,
        const Matrix2D<double>& VP,
        const Matrix2D<double>& VR,
        const Matrix2D<double>& Loss,
        const Observations& observations,
        double timeStep,
        const Parameters& parameters
    );

    /**
     * @brief Applies the Kalman filter updating forecast and Pa
     * @param[in, out]  forecast    the current forecast [n x 1]
     * @param[in, out]  Pa          analysis error covariance matrix [n x n]
     * @param[in]       M           model matrix operator [n x n]
     * @param[in]       Q           model error corvariance matrix [n x n]
     * @param[in]       obs         observations [m x 1]
     * @param[in]       H           model operator, transform forecast to observation space [m x n]
     * @param[in]       R           observation error covariance matrix [m x m]
     */
    void runKalmanFilter(
        Matrix1D<double>& forecast,  
        Matrix2D<double>& Pa,
        const Matrix2D<double>& M, const Matrix2D<double>& Q, 
        const Matrix1D<double>& obs, const Matrix2D<double>& H, 
        const Matrix2D<double>& R
    );
}
