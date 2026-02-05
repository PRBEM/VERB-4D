#pragma once

#include "CustomDate.h"
#include "Matrix.h"

#include <string>
#include <vector>

namespace pmf {

enum class MagneticFieldModel {NoExtB, T89, T96, T01, T01S, T04S, TS07Dmid15};

enum class DataVersion {v4};

struct Parameters {
    std::string dataServer;
    std::string mission;
    std::string satellite;
    std::string instrument;
    MagneticFieldModel magneticFieldModel;
    DataVersion version;
};

std::vector<Parameters> readParameters(const std::string& filename);

Matrix1D<double> readProcessedMatFiles1D(
    const std::string& name, 
    CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters
); 

Matrix2D<double> readProcessedMatFiles2D(
    const std::string& name, 
    CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters
); 

Matrix3D<double> readProcessedMatFiles3D(
    const std::string& name, 
    CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters
); 

} 