
#include "PMFInternal.h"
#include "test_runner.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

using std::string;

using namespace pmf::internal;

void TestPMFInternalAuxiliary() {
    std::string str = "AbCdEf";
    std::string expected = "ABCDEF";
    ASSERT_EQUAL(expected, toupper(str));
}

void TestPMFMissionFolderGeneration() {
    string expected = "RBSP";
    ASSERT_EQUAL(expected, generateProcessedMatFileMissionFolder("rbsp"));
    ASSERT_EQUAL(expected, generateProcessedMatFileMissionFolder("RBSP"));
    ASSERT_EQUAL(expected, generateProcessedMatFileMissionFolder("rBsP"));
}

void TestPMFSatelliteFolderGeneration() {
    string expected = "rbspa";
    ASSERT_EQUAL(expected, generateProcessedMatFileSatelliteFolder("rbspa"));
    ASSERT_EQUAL(expected, generateProcessedMatFileSatelliteFolder("RBSPA"));
    ASSERT_EQUAL(expected, generateProcessedMatFileSatelliteFolder("rBsPa"));
    
    expected = "rbspb";
    ASSERT_EQUAL(expected, generateProcessedMatFileSatelliteFolder("rbspb"));
    ASSERT_EQUAL(expected, generateProcessedMatFileSatelliteFolder("RBSPB"));
    ASSERT_EQUAL(expected, generateProcessedMatFileSatelliteFolder("rBsPb"));
}

void TestPMFSatelliteNameGeneration() {
    string expected = "rbspa";
    ASSERT_EQUAL(expected, generateProcessedMatFileSatelliteName("rbspa"));
    ASSERT_EQUAL(expected, generateProcessedMatFileSatelliteName("RBSPA"));
    ASSERT_EQUAL(expected, generateProcessedMatFileSatelliteName("rBsPa"));
    
    expected = "rbspb";
    ASSERT_EQUAL(expected, generateProcessedMatFileSatelliteName("rbspb"));
    ASSERT_EQUAL(expected, generateProcessedMatFileSatelliteName("RBSPB"));
    ASSERT_EQUAL(expected, generateProcessedMatFileSatelliteName("rBsPb"));
}

void TestPMFInstrumentNameGeneration() {
    string expected = "hope";
    ASSERT_EQUAL(expected, generateProcessedMatFileInstrumentName("hope"));
    ASSERT_EQUAL(expected, generateProcessedMatFileInstrumentName("HOPE"));
    ASSERT_EQUAL(expected, generateProcessedMatFileInstrumentName("hOpE"));
    
    expected = "mageis";
    ASSERT_EQUAL(expected, generateProcessedMatFileInstrumentName("mageis"));
    ASSERT_EQUAL(expected, generateProcessedMatFileInstrumentName("MAGEIS"));
    ASSERT_EQUAL(expected, generateProcessedMatFileInstrumentName("MaGeIs"));
    
    expected = "rept";
    ASSERT_EQUAL(expected, generateProcessedMatFileInstrumentName("rept"));
    ASSERT_EQUAL(expected, generateProcessedMatFileInstrumentName("REPT"));
    ASSERT_EQUAL(expected, generateProcessedMatFileInstrumentName("RePt"));
}

void TestPMFVariableSuffixGeneration() {
    ASSERT_EQUAL("xGEO", generateProcessedMatFileVariableSuffix("time"));
    ASSERT_EQUAL("xGEO", generateProcessedMatFileVariableSuffix("TIME"));
    ASSERT_EQUAL("xGEO", generateProcessedMatFileVariableSuffix("TimE"));

    ASSERT_EQUAL("mlt", generateProcessedMatFileVariableSuffix("mlt"));
    ASSERT_EQUAL("mlt", generateProcessedMatFileVariableSuffix("MLT"));
    ASSERT_EQUAL("mlt", generateProcessedMatFileVariableSuffix("mLt"));

    ASSERT_EQUAL("xGEO", generateProcessedMatFileVariableSuffix("xgeo"));
    ASSERT_EQUAL("xGEO", generateProcessedMatFileVariableSuffix("XGEO"));
    ASSERT_EQUAL("xGEO", generateProcessedMatFileVariableSuffix("xGeO"));

    ASSERT_EQUAL("invmu_and_invk", 
        generateProcessedMatFileVariableSuffix("invmu"));
    ASSERT_EQUAL("invmu_and_invk", 
        generateProcessedMatFileVariableSuffix("INVMU"));
    ASSERT_EQUAL("invmu_and_invk", 
        generateProcessedMatFileVariableSuffix("InVmU"));

    ASSERT_EQUAL("invmu_and_invk", 
        generateProcessedMatFileVariableSuffix("invk"));
    ASSERT_EQUAL("invmu_and_invk", 
        generateProcessedMatFileVariableSuffix("INVK"));
    ASSERT_EQUAL("invmu_and_invk", 
        generateProcessedMatFileVariableSuffix("InVk"));

    ASSERT_EQUAL("psd", generateProcessedMatFileVariableSuffix("psd"));
    ASSERT_EQUAL("psd", generateProcessedMatFileVariableSuffix("PSD"));
    ASSERT_EQUAL("psd", generateProcessedMatFileVariableSuffix("pSd"));
}

void TestPMFMagneticFieldModelNameGeneration() {
    ASSERT_EQUAL("n4_4_NoExtB", generateProcessedMatFileMagneticFieldModelName(
                                            pmf::MagneticFieldModel::NoExtB));
    ASSERT_EQUAL("n4_4_T89", generateProcessedMatFileMagneticFieldModelName(
                                            pmf::MagneticFieldModel::T89));
    ASSERT_EQUAL("n4_4_T96", generateProcessedMatFileMagneticFieldModelName(
                                            pmf::MagneticFieldModel::T96));
    ASSERT_EQUAL("n4_4_T01", generateProcessedMatFileMagneticFieldModelName(
                                            pmf::MagneticFieldModel::T01));
    ASSERT_EQUAL("n4_4_T01s", generateProcessedMatFileMagneticFieldModelName(
                                            pmf::MagneticFieldModel::T01S));
    ASSERT_EQUAL("n4_4_T04s", generateProcessedMatFileMagneticFieldModelName(
                                            pmf::MagneticFieldModel::T04S));
    ASSERT_EQUAL("n4_4_TS07Dmid15", generateProcessedMatFileMagneticFieldModelName(
                                            pmf::MagneticFieldModel::TS07Dmid15));
}

void TestPMFVersionNameGeneration() {
    ASSERT_EQUAL("ver4", 
            generateProcessedMatFileVersionName(pmf::DataVersion::v4));
}

void TestPMFFullNameGeneration() {
    const pmf::Parameters parameters {
        .dataServer = "/home/wutzig/mnt/data",
        .mission = "rbsp",
        .satellite = "rbspa",
        .instrument = "mageis",
        .magneticFieldModel = pmf::MagneticFieldModel::T89,
        .version = pmf::DataVersion::v4
    };

    string variableName = "invk";
    const CustomDate timeStart {2013, 3, 1};
    const CustomDate timeEnd {2013, 3, 31};

    string expected = "/home/wutzig/mnt/data/RBSP/rbspa/Processed_Mat_Files/"
                  "rbspa_mageis_20130301to20130331_invmu_and_invk_n4_4_T89_ver4.mat";
    ASSERT_EQUAL(expected, 
        generateProcessedMatFileName(
            timeStart, timeEnd, parameters, variableName));

    variableName = "psd";
    expected = "/home/wutzig/mnt/data/RBSP/rbspa/Processed_Mat_Files/"
                  "rbspa_mageis_20130301to20130331_psd_ver4.mat";
    ASSERT_EQUAL(expected, 
        generateProcessedMatFileName(
            timeStart, timeEnd, parameters, variableName));
    
}

void TestPMFVariableNameGeneration() {
    ASSERT_EQUAL("time", generateProcessedMatFileVariableName("time"));
    ASSERT_EQUAL("time", generateProcessedMatFileVariableName("TIME"));
    ASSERT_EQUAL("time", generateProcessedMatFileVariableName("TimE"));

    ASSERT_EQUAL("MLT", generateProcessedMatFileVariableName("mlt"));
    ASSERT_EQUAL("MLT", generateProcessedMatFileVariableName("MLT"));
    ASSERT_EQUAL("MLT", generateProcessedMatFileVariableName("mLt"));

    ASSERT_EQUAL("xGEO", generateProcessedMatFileVariableName("xgeo"));
    ASSERT_EQUAL("xGEO", generateProcessedMatFileVariableName("XGEO"));
    ASSERT_EQUAL("xGEO", generateProcessedMatFileVariableName("xGeO"));

    ASSERT_EQUAL("InvMu", 
        generateProcessedMatFileVariableName("invmu"));
    ASSERT_EQUAL("InvMu", 
        generateProcessedMatFileVariableName("INVMU"));
    ASSERT_EQUAL("InvMu", 
        generateProcessedMatFileVariableName("InVmU"));

    ASSERT_EQUAL("InvK", 
        generateProcessedMatFileVariableName("invk"));
    ASSERT_EQUAL("InvK", 
        generateProcessedMatFileVariableName("INVK"));
    ASSERT_EQUAL("InvK", 
        generateProcessedMatFileVariableName("InVk"));

    ASSERT_EQUAL("PSD", generateProcessedMatFileVariableName("psd"));
    ASSERT_EQUAL("PSD", generateProcessedMatFileVariableName("PSD"));
    ASSERT_EQUAL("PSD", generateProcessedMatFileVariableName("pSd"));
}

void TestPMFGetMagneticFieldModel() {
    ASSERT(pmf::MagneticFieldModel::NoExtB == getMagneticFieldModel("n4_4_NoExtB"));
    ASSERT(pmf::MagneticFieldModel::T89 == getMagneticFieldModel("n4_4_T89"));
    ASSERT(pmf::MagneticFieldModel::T96 == getMagneticFieldModel("n4_4_T96"));
    ASSERT(pmf::MagneticFieldModel::T01 == getMagneticFieldModel("n4_4_T01"));
    ASSERT(pmf::MagneticFieldModel::T01S == getMagneticFieldModel("n4_4_T01S"));
    ASSERT(pmf::MagneticFieldModel::T04S == getMagneticFieldModel("n4_4_T04S"));
    ASSERT(pmf::MagneticFieldModel::TS07Dmid15 == getMagneticFieldModel("n4_4_TS07Dmid15"));
}

void TestPMFGetVersion() {
    ASSERT(pmf::DataVersion::v4 == getPmfVersion("ver4"));
}

void TestParameterReading() {
    pmf::Parameters expected {
        .dataServer = "/home/wutzig/mnt/data",
        .mission = "rbsp",
        .satellite = "rbspa",
        .instrument = "hope",
        .magneticFieldModel = pmf::MagneticFieldModel::T89,
        .version = pmf::DataVersion::v4
    };

    string filename {"ParametersTMP.ini"};
    std::ofstream file (filename);
    if (file.is_open()) {
        file << "data_server = /home/wutzig/mnt/data" << std::endl; 
        file << "satellite = rbspa" << std::endl;
        file << "mission = rbsp" << std::endl;
        file << "instrument = hope" << std::endl;
        file << "magnetic_field = n4_4_T89" << std::endl;
        file << "version = ver4" << std::endl;

        pmf::Parameters result = readParameterOneSatellite(filename);
        std::remove(filename.c_str());

        ASSERT_EQUAL(expected.dataServer, result.dataServer);
        ASSERT_EQUAL(expected.mission, result.mission);
        ASSERT_EQUAL(expected.satellite, result.satellite);
        ASSERT_EQUAL(expected.instrument, result.instrument);
        ASSERT(expected.magneticFieldModel == result.magneticFieldModel);
        ASSERT(expected.version == result.version);
    } else {
        std::cerr << "Cannot open file for writing: " << filename << std::endl;
    }

}

void TestPMFReading() {
    pmf::Parameters parameters {
        .dataServer = "/home/wutzig/mnt/data",
        .mission = "rbsp",
        .satellite = "rbspa",
        .instrument = "hope",
        .magneticFieldModel = pmf::MagneticFieldModel::T89,
        .version = pmf::DataVersion::v4
    };

    CustomDate timeStart {2013, 3, 1};
    CustomDate timeEnd {2013, 3, 31};

    Matrix3D<double> mu = readOneProcessedMatFile<Matrix3D<double>>(
         "InvMu", timeStart, timeEnd, parameters); 
    Matrix1D<double> time = readOneProcessedMatFile<Matrix1D<double>>(
         "time", timeStart, timeEnd, parameters); 
    int a;
}

void TestPMFInternal() {
    TestRunner tr;
    RUN_TEST(tr, TestPMFInternalAuxiliary);
    RUN_TEST(tr, TestPMFMissionFolderGeneration);
    RUN_TEST(tr, TestPMFSatelliteFolderGeneration);
    RUN_TEST(tr, TestPMFSatelliteNameGeneration);
    RUN_TEST(tr, TestPMFInstrumentNameGeneration);
    RUN_TEST(tr, TestPMFVariableSuffixGeneration);
    RUN_TEST(tr, TestPMFMagneticFieldModelNameGeneration);
    RUN_TEST(tr, TestPMFVersionNameGeneration);
    RUN_TEST(tr, TestPMFFullNameGeneration);
    RUN_TEST(tr, TestPMFVariableNameGeneration);
    RUN_TEST(tr, TestPMFGetMagneticFieldModel);
    RUN_TEST(tr, TestPMFGetVersion);
    RUN_TEST(tr, TestParameterReading);
    RUN_TEST(tr, TestPMFReading);
}
