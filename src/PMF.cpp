#include <iostream>
#include <string>
#include <vector>
#include "mat.h"
#include "matrix.h"

#include "MatrixOperations.h"
#include "Parameters.h"
#include "PMF.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

using ParametersIni=::Parameters;

namespace pmf::internal {};

using namespace pmf::internal;

namespace pmf::internal {
string toupper(string str) {
    for (size_t i = 0; i<str.size(); ++i)
        str[i] = std::toupper(str[i]);

    return str;
}
}

namespace pmf::internal {

class ProcessedMatFile {
    MATFile* file = NULL;

    ProcessedMatFile();
    ProcessedMatFile operator=(const ProcessedMatFile&);
    ProcessedMatFile(const ProcessedMatFile& );
public:
    ProcessedMatFile(std::string path2file){
        file = matOpen(path2file.c_str(), "r");
        if (file == NULL){
            std::cout << "In ProcessedMatFile: cannot open file \"" << 
                path2file << "\"" << "\n";
            exit(EXIT_FAILURE);
        }
    }
    ~ProcessedMatFile() { if (file != NULL) matClose(file); }

    operator MATFile*() const { return file; }
};

class MatArray {
    std::vector<mwSize> dims;
    mxArray *mx_arr = NULL;
    double* elem = NULL;
    std::string name;

    MatArray();
    MatArray operator=(const MatArray&);
    MatArray (const MatArray&);
public:
    MatArray(const ProcessedMatFile& pmf, const std::string& arr_name);
    ~MatArray() { if (mx_arr != NULL) mxDestroyArray(mx_arr); }

    double operator[] (int i) const {return *(elem + i);};
    std::vector<int> size() const;
    int length() const;

    Matrix1D<double> to_Matrix1D() const;
    Matrix2D<double> to_Matrix2D() const;
    Matrix3D<double> to_Matrix3D() const;
};

MatArray::MatArray(const ProcessedMatFile& pmf, const std::string& arr_name)
{
    mx_arr = matGetVariable(pmf, arr_name.c_str());
    if (mx_arr == NULL) {
        std::cout << "In MatArray: field " << arr_name << " does not exist in processed mat file\n";
        exit(EXIT_FAILURE);
    } else {
        name = arr_name;
        mwSize ndim = mxGetNumberOfDimensions(mx_arr);
        dims = std::vector<mwSize>(ndim);
        for (mwSize i = 0; i<ndim; ++i){
            dims[i] = *(mxGetDimensions(mx_arr) + i);
            //std::cout << "mxGETDIMENSIONS, " << i << ": " << dims[i] << std::endl;
        }
        elem = mxGetPr(mx_arr);
    }
}

std::vector<int> MatArray::size() const {
    std::vector<int> size_arr(dims.size());
    for (size_t i = 0; i<size_arr.size(); ++i) {
        size_arr[i] = static_cast<int>(dims[i]);
    }

    return size_arr;
}

int MatArray::length() const {
    int len = 1;
    for (auto i = 0; i < dims.size(); ++i)
        len *= static_cast<int>(dims[i]);

    return len;
}

Matrix1D<double> MatArray::to_Matrix1D() const {
    if (dims.size() != 1){
        for (int i = 1; i < dims.size(); ++i){
            if (dims[i] != 1) {
                std::cout << this->name << std::endl;
                std::cout << "Cannot convert " << dims.size() << "D Matlab array to Matrix1D\n";
                exit(EXIT_FAILURE);
            }
        }
    }

    Matrix1D<double> matrix;

    matrix.AllocateMemory(this->length());

    for (int i=0; i < matrix.size_q1; ++i){
        matrix[i] = *(elem + i);
    }

    return matrix;

}

Matrix2D<double> MatArray::to_Matrix2D() const {
    if (dims.size() != 2) {
        std::cout << "Cannot convert " << dims.size() << "D Matlab array to Matrix2D\n";
        exit(EXIT_FAILURE);
    } else {

        Matrix2D<double> matrix;

        matrix.AllocateMemory(dims[0], dims[1]);

        int lin_idx;
        for (int ix = 0; ix < dims[0]; ++ix)
            for (int iy = 0; iy < dims[1]; ++iy) {
                lin_idx = ix + iy * dims[0];
                matrix[ix][iy] = *(elem + lin_idx);
            }

        return matrix;
    }
}

Matrix3D<double> MatArray::to_Matrix3D() const {
    if (dims.size() != 3) {
        std::cout << "Cannot convert " << dims.size() << "D Matlab array to Matrix2D\n";
        exit(EXIT_FAILURE);
    } else {
        Matrix3D<double> matrix;

        matrix.AllocateMemory(dims[0], dims[1], dims[2]);

        int lin_idx;
        for (int ix = 0; ix < dims[0]; ++ix)
            for (int iy = 0; iy < dims[1]; ++iy)
                for (int iz = 0; iz < dims[2]; ++iz){
                    lin_idx = ix + iy * dims[0] + iz * dims[0] * dims[1];
                    matrix[ix][iy][iz] = *(elem + lin_idx);
                }

        return matrix;
    }
}

template<typename MatrixType> 
MatrixType cast_to(const MatArray&);

template<>
Matrix1D<double> cast_to<Matrix1D<double>>(const MatArray& array) {
    return array.to_Matrix1D();
}

template<>
Matrix2D<double> cast_to<Matrix2D<double>>(const MatArray& array) {
    return array.to_Matrix2D();
}
template<>
Matrix3D<double> cast_to<Matrix3D<double>>(const MatArray& array) {
    return array.to_Matrix3D();
}
}

namespace pmf::internal {

string generateProcessedMatFileMissionFolder (const string& mission) {
    string missionUpper = toupper(mission);
    if (missionUpper == "RBSP" || missionUpper == "ARASE") {
        return missionUpper;
    } else {
        std::cout << "Mission " << mission << " is not supported\n";
        exit(EXIT_FAILURE);
    }
}

string generateProcessedMatFileSatelliteFolder(const string& satellite) {
    string satellite_upper = toupper(satellite);

    if (satellite_upper == "RBSPA") {
        return "rbspa";
    } else if (satellite_upper == "RBSPB") {
        return "rbspb";
    } else if (satellite_upper == "ARASE") {
        return "arase";
    } else if (satellite_upper == "") {
        std::cout << "Empty satellite name!\n";
        exit(EXIT_FAILURE);
    } else {
        std::cout << "Unknown satellite name: " << satellite << "\n";
        exit(EXIT_FAILURE);
    }
}

string generateProcessedMatFileSatelliteName(const string& satellite) {
    string satellite_upper = toupper(satellite);

    if (satellite_upper == "RBSPA") {
        return "rbspa";
    } else if (satellite_upper == "RBSPB") {
        return "rbspb";
    } else if (satellite_upper == "ARASE") {
        return "arase";
    }  else if (satellite_upper == "") {
        std::cout << "Empty satellite name!\n";
        exit(EXIT_FAILURE);
    } else {
        std::cout << "Unknown satellite name: " << satellite << "\n";
        exit(EXIT_FAILURE);
    }
}

string generateProcessedMatFileInstrumentName(const string& instrument) {

    string instrument_upper = toupper(instrument);

    if (instrument_upper == "HOPE") {
        return "hope";
    } else if (instrument_upper == "MAGEIS") {
        return "mageis";
    } else if (instrument_upper == "REPT") {
        return "rept";
    } else if (instrument_upper == "MEPE-L3") {
        return "mepe-l3";
    } else if (instrument_upper == "HEP-L3L") {
        return "hep-l3l";
    } else if (instrument_upper == "") {
        cout << "Empty instrument name!\n";
        exit(EXIT_FAILURE);
    } else {
        cout << "Unknown instrument " << instrument << endl;
        exit(EXIT_FAILURE);
    }
}

string generateProcessedMatFileVariableSuffix(const string& name) {
    string name_upper = toupper(name);
    if (name_upper == "TIME") {
        return "xGEO";
    } else if (name_upper == "MLT") {
        return "mlt";
    } else if (name_upper == "R0") {
        return "R0";
    } else if (name_upper == "XGEO") {
        return "xGEO";
    } else if (name_upper == "INVMU") {            
        return "invmu_and_invk";
    } else if (name_upper == "INVK") {
        return "invmu_and_invk";
    } else if (name_upper == "PSD") {
        return "psd";
    } else {
        cout << "Error! Uknown processed mat file varible name: " << 
            name << endl;
        exit(EXIT_FAILURE);
    }
}

string generateProcessedMatFileMagneticFieldModelName(pmf::MagneticFieldModel mfm) {

    switch (mfm) {
        case pmf::MagneticFieldModel::NoExtB:
            return "n4_4_NoExtB";
        case pmf::MagneticFieldModel::T89:
            return "n4_4_T89";
        case pmf::MagneticFieldModel::T96:
            return "n4_4_T96";
        case pmf::MagneticFieldModel::T01:
            return "n4_4_T01";
        case pmf::MagneticFieldModel::T01S:
            return "n4_4_T01s";
        case pmf::MagneticFieldModel::T04S:
            return "n4_4_T04s";
        case pmf::MagneticFieldModel::TS07Dmid15:
            return "n4_4_TS07Dmid15";
    default:
        cout << "Error in " << __FILE__ << ", line " << __LINE__ << ": " <<
            "Magnetic field model is not supported." << endl;
        exit(EXIT_FAILURE);
    }            

}

string generateProcessedMatFileVersionName(pmf::DataVersion version) {
    switch (version) {
    case pmf::DataVersion::v4:
        return "ver4";    
    default:
        cout << "Error in " << __FILE__ << ", line " << __LINE__ << ": " <<
            "Magnetic field model is not supported." << endl;
        exit(EXIT_FAILURE);
    }
}

string generateProcessedMatFileName(
    CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters,
        const string& name) {
            
    char filesep {'/'};

    string path2server = parameters.dataServer;
    string missionFolder = generateProcessedMatFileMissionFolder(
        parameters.mission
    );
    string satelliteFolder = generateProcessedMatFileSatelliteFolder(
        parameters.satellite
    );
    string processedMatFilesFolder = "Processed_Mat_Files";
    string satelliteName = generateProcessedMatFileSatelliteName(
        parameters.satellite
    );
    string instrumentName = generateProcessedMatFileInstrumentName(
        parameters.instrument
    );
    string dateStringBegin = timeStart.to_date_string();
    string dateStringEnd = timeEnd.to_date_string();
    string variableSuffix = generateProcessedMatFileVariableSuffix(
        name
    );
    string mfmName = generateProcessedMatFileMagneticFieldModelName(
        parameters.magneticFieldModel
    );
    string versionName = generateProcessedMatFileVersionName(
        parameters.version
    );

    string generalPrefix = path2server + filesep + missionFolder + filesep + 
        satelliteFolder + filesep + processedMatFilesFolder + filesep + 
        satelliteName + '_' + instrumentName + '_' + 
        dateStringBegin + "to" + dateStringEnd + '_' + variableSuffix + '_';

    if (variableSuffix == "psd" || variableSuffix == "xGEO") {
        return generalPrefix + versionName + ".mat";
    } else {
        return generalPrefix + mfmName + '_' + versionName + ".mat";
    }
}

string generateProcessedMatFileVariableName(const string& name) {
    string name_upper = toupper(name);
    if (name_upper == "TIME") {
        return "time";
    } else if (name_upper == "MLT") {
        return "MLT";
    } else if (name_upper == "R0") {
        return "R0";
    } else if (name_upper == "XGEO") {
        return "xGEO";
    } else if (name_upper == "INVMU") {            
        return "InvMu";
    } else if (name_upper == "INVK") {
        return "InvK";
    } else if (name_upper == "PSD") {
        return "PSD";
    } else {
        cout << "Error! In function " << __FUNCTION__ << " file " << __FILE__ 
                << " line " << __LINE__ 
                << " Uknown processed mat file varible name: " << 
            name << endl;
        exit(EXIT_FAILURE);
    }
}

pmf::MagneticFieldModel getMagneticFieldModel(const string& name) {
    if (name == "n4_4_NoExtB") {
        return pmf::MagneticFieldModel::NoExtB;
    } else if (name == "n4_4_T89") {
        return pmf::MagneticFieldModel::T89;
    } else if (name == "n4_4_T96") {
        return pmf::MagneticFieldModel::T96;
    } else if (name == "n4_4_T01") {
        return pmf::MagneticFieldModel::T01;
    } else if (name == "n4_4_T01S") {
        return pmf::MagneticFieldModel::T01S;
    } else if (name == "n4_4_T04S") {
        return pmf::MagneticFieldModel::T04S;
    } else if (name == "n4_4_TS07Dmid15") {
        return pmf::MagneticFieldModel::TS07Dmid15;
    } else {
        cout << "Error! In function " << __FUNCTION__ << " file " << __FILE__ 
                << " line " << __LINE__ 
                << " Unknown magnetic field model: " << name << endl;
        exit(EXIT_FAILURE);
    }
}

pmf::DataVersion getPmfVersion(const string& name) {
    if (name == "ver4") {
        return pmf::DataVersion::v4;
    } else {
        cout << "Error! In function " << __FUNCTION__ << " file " << __FILE__ 
                << " line " << __LINE__ 
                << " Unknown version of processed mat file: " << name << endl;
        exit(EXIT_FAILURE);
    }
}

pmf::Parameters readParameterOneSatellite(const string& filename) {

    pmf::Parameters result;

    ParametersIni parameters(filename);

	parameters.getParameter("data_server", result.dataServer, true);
	parameters.getParameter("satellite", result.satellite, true);
	parameters.getParameter("mission", result.mission, true);
	parameters.getParameter("instrument", result.instrument, true);

    string magneticFieldString;
	parameters.getParameter("magnetic_field", magneticFieldString, true);
    result.magneticFieldModel = getMagneticFieldModel(magneticFieldString);

    string versionString;
	parameters.getParameter("version", versionString, true);
    result.version = getPmfVersion(versionString);

    return result;
}

}

std::vector<pmf::Parameters> pmf::readParameters(const std::string& filename) {
    std::fstream filestream(filename);    
    std::vector<pmf::Parameters> result; 
    for(std::string filename; std::getline(filestream, filename); ) {
        result.push_back(readParameterOneSatellite(filename));
    }
    return result;
}

namespace pmf::internal {

template<typename MatrixType>
MatrixType limitWithTime(
    const MatrixType& array, const Matrix1D<double>& time,
    CustomDate timeStart, CustomDate timeEnd, uint32_t memory_seconds = 0);              

template<>
Matrix1D<double> limitWithTime<Matrix1D<double>>(
    const Matrix1D<double>& array, const Matrix1D<double>& time,
    CustomDate timeStart, CustomDate timeEnd, uint32_t memory_seconds) {               
        
    if (array.size_q1 != time.size_q1) {
        cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        cout << "Sizes mismatch:";
        cout << "array.size_q1 = " << array.size_q1 << ", ";
        cout << "time.size_q1 = " << time.size_q1 << endl;
        exit(EXIT_FAILURE);
    }

    int size = {0};
    for (auto it = 0; it < time.size_q1; ++it) {
        CustomDate cd {time[it]};
        if (cd >= timeStart && cd <= timeEnd) {
            ++size;
        }
    }

    Matrix1D<double> result(size);
    int counter {0};
    for (auto it = 0; it < time.size_q1; ++it) {
        CustomDate cd {time[it]};
        if (cd >= timeStart && cd <= timeEnd) {
           result[counter] = array[it]; 
           ++counter;
        }
    }

    return result;
}

template<>
Matrix2D<double> limitWithTime<Matrix2D<double>>(
    const Matrix2D<double>& array, const Matrix1D<double>& time,
    CustomDate timeStart, CustomDate timeEnd, uint32_t memory_seconds) {               
        
    if (array.size_q1 != time.size_q1) {
        cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        cout << "Sizes mismatch:";
        cout << "array.size_q1 = " << array.size_q1 << ", ";
        cout << "time.size_q1 = " << time.size_q1 << endl;
        exit(EXIT_FAILURE);
    }

    int size = {0};
    for (auto it = 0; it < time.size_q1; ++it) {
        CustomDate cd {time[it]};
        if (cd >= timeStart && cd <= timeEnd) {
            ++size;
        }
    }

    Matrix2D<double> result(size, array.size_q2);
    int counter {0};
    for (auto it = 0; it < time.size_q1; ++it) {
        CustomDate cd {time[it]};
        if (cd >= timeStart && cd <= timeEnd) {
           for (auto i = 0; i < result.size_q2; ++i) {
            result[counter][i] = array[it][i]; 
           }
           ++counter;
        }
    }

    return result;
}

template<>
Matrix3D<double> limitWithTime<Matrix3D<double>>(
        const Matrix3D<double>& array, const Matrix1D<double>& time,
        CustomDate timeStart, CustomDate timeEnd, uint32_t memory_seconds) {               
    timeStart = timeStart + memory_seconds;
    timeEnd   = timeEnd   - memory_seconds;
    if (array.size_q1 != time.size_q1) {
        cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        cout << "Sizes mismatch:";
        cout << "array.size_q1 = " << array.size_q1 << ", ";
        cout << "time.size_q1 = " << time.size_q1 << endl;
        exit(EXIT_FAILURE);
    }

    int size = {0};
    for (auto it = 0; it < time.size_q1; ++it) {
        CustomDate cd {time[it]};
        if (cd >= timeStart && cd <= timeEnd) {
            ++size;
        }
    }

    Matrix3D<double> result(size, array.size_q2, array.size_q3);
    int counter {0};
    for (auto it = 0; it < time.size_q1; ++it) {
        CustomDate cd {time[it]};
        if (cd >= timeStart && cd <= timeEnd) {
           for (int i = 0; i < result.size_q2; ++i) {
               for (int j = 0; j < result.size_q3; ++j) {
                result[counter][i][j] = array[it][i][j]; 
               }
           }
           ++counter;
        }
    }
    return result;
}

template<typename MatrixType>
MatrixType readOneProcessedMatFile(
    const string& name, 
    CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters
) {
    std::string filename = generateProcessedMatFileName(
        timeStart, timeEnd, parameters, name);

    string variableName = generateProcessedMatFileVariableName(name);
    ProcessedMatFile pmf(filename);    
    MatArray array(pmf, variableName);   
    MatrixType result = cast_to<MatrixType>(array);
    return result;
}

template<typename MatrixType>
MatrixType readProcessedMatFiles(
    const string& name, 
    CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters
) {
    CustomDate sd {timeStart.to_bom()};
    CustomDate ed {timeEnd.to_eom()};

    MatrixType result;
    while(sd < ed) {
        CustomDate end_current_month = sd.to_eom();
        MatrixType matrix_1m = readOneProcessedMatFile<MatrixType>(
            name, sd, end_current_month, parameters);

        if (sd < timeStart || end_current_month > timeEnd) {
            Matrix1D<double> time = readOneProcessedMatFile<Matrix1D<double>>(                   
                "time", sd, end_current_month, parameters);
            matrix_1m = limitWithTime<MatrixType>(matrix_1m, time, 
                timeStart, timeEnd);
        }

        result = cat(result, matrix_1m);
        sd = sd.add_month();
    }
    return result;
}

}

Matrix1D<double> pmf::readProcessedMatFiles1D(
    const string& name, 
    CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters
) {
    Matrix1D<double> result = readProcessedMatFiles<Matrix1D<double>>(
            name, timeStart, timeEnd, parameters);
    return result;
}

Matrix2D<double> pmf::readProcessedMatFiles2D(
    const string& name, 
    CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters
) {
    auto result = readProcessedMatFiles<Matrix2D<double>>(
            name, timeStart, timeEnd, parameters);
    return result;
}

Matrix3D<double> pmf::readProcessedMatFiles3D(
    const string& name, 
    CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters
) {
    auto result = readProcessedMatFiles<Matrix3D<double>>(
            name, timeStart, timeEnd, parameters);
    return result;
}
