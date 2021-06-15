#pragma once

#include "../src/PMF.h"
#include <string>
#include <vector>

namespace pmf::internal {
    class ProcessedMatFile;
    
    class MatArray;
    
    std::string toupper(std::string str);
    
    std::string generateProcessedMatFileMissionFolder (
            const std::string& mission);
    
    std::string generateProcessedMatFileSatelliteFolder(const std::string& satellite);
    
    std::string generateProcessedMatFileSatelliteName(const std::string& satellite);
    
    std::string generateProcessedMatFileInstrumentName(const std::string& instrument);
    
    std::string generateProcessedMatFileVariableSuffix(const std::string& name);
    
    std::string generateProcessedMatFileMagneticFieldModelName(
        pmf::MagneticFieldModel mfm);
    
    std::string generateProcessedMatFileVersionName(pmf::DataVersion version);

    std::string generateProcessedMatFileName(
        CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters,
            const std::string& name);

    std::string generateProcessedMatFileVariableName(const std::string& name);

    pmf::MagneticFieldModel getMagneticFieldModel(const std::string& name);

    pmf::DataVersion getPmfVersion(const std::string& name);

    pmf::Parameters readParameterOneSatellite(const std::string& filename);
    
    template<typename MatrixType>
    MatrixType limitWithTime(
        const MatrixType& array, const Matrix1D<double>& time,
        CustomDate timeStart, CustomDate timeEnd);              

    template<typename MatrixType>
    MatrixType readOneProcessedMatFile(
        const std::string& name, 
        CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters
    );

    template<typename MatrixType>
    MatrixType readProcessedMatFiles(
        const std::string& name, 
        CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters
    );
}
