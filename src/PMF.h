/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * @file PMF.h
 * @brief Processed MAT File (PMF) reader for satellite data
 * 
 * This file provides functionality to read processed MAT files containing
 * satellite data from various missions (RBSP, ARASE) and instruments.
 * Supports different magnetic field models and data versions.
 */

#pragma once

#include "CustomDate.h"
#include "Matrix.h"

#include <string>
#include <vector>

/**
 * @namespace pmf
 * @brief Namespace for Processed MAT File operations
 */
namespace pmf {

/**
 * @enum MagneticFieldModel
 * @brief Supported magnetic field models for data processing
 */
enum class MagneticFieldModel {
    NoExtB,      ///< No external magnetic field
    T89,         ///< Tsyganenko 1989 model
    T96,         ///< Tsyganenko 1996 model
    T01,         ///< Tsyganenko 2001 model
    T01S,        ///< Tsyganenko 2001 storm model
    T04S,        ///< Tsyganenko 2004 storm model
    TS07Dmid15   ///< Tsyganenko-Sitnov 2007 model with Dmid=15
};

/**
 * @enum DataVersion
 * @brief Supported data versions for processed MAT files
 */
enum class DataVersion {
    v4  ///< Version 4 of processed data format
};

/**
 * @struct Parameters
 * @brief Configuration parameters for reading processed MAT files
 */
struct Parameters {
    std::string dataServer;              ///< Path to data server containing MAT files
    std::string mission;                 ///< Mission name (e.g., "RBSP", "ARASE")
    std::string satellite;               ///< Satellite identifier (e.g., "RBSPA", "RBSPB", "ARASE")
    std::string instrument;              ///< Instrument name (e.g., "HOPE", "MAGEIS", "REPT")
    MagneticFieldModel magneticFieldModel; ///< Magnetic field model used for processing
    DataVersion version;                 ///< Data format version
};

/**
 * @brief Read PMF parameters from configuration file
 * 
 * @param filename Path to configuration file containing parameter definitions
 * @return Vector of Parameters structures for multiple satellites/instruments
 */
std::vector<Parameters> readParameters(const std::string& filename);

/**
 * @brief Read 1D data array from processed MAT files
 * 
 * @param name Variable name to read (e.g., "MLT", "R0", "time")
 * @param timeStart Start time for data extraction
 * @param timeEnd End time for data extraction
 * @param parameters PMF parameters specifying data source and processing options
 * @return 1D matrix containing the requested data over the time range
 */
Matrix1D<double> readProcessedMatFiles1D(
    const std::string& name, 
    CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters
); 

/**
 * @brief Read 2D data array from processed MAT files
 * 
 * @param name Variable name to read (e.g., "InvK" for K parameter data)
 * @param timeStart Start time for data extraction
 * @param timeEnd End time for data extraction
 * @param parameters PMF parameters specifying data source and processing options
 * @return 2D matrix containing the requested data over the time range
 */
Matrix2D<double> readProcessedMatFiles2D(
    const std::string& name, 
    CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters
); 

/**
 * @brief Read 3D data array from processed MAT files
 * 
 * @param name Variable name to read (e.g., "PSD", "InvMu" for phase space density data)
 * @param timeStart Start time for data extraction
 * @param timeEnd End time for data extraction
 * @param parameters PMF parameters specifying data source and processing options
 * @return 3D matrix containing the requested data over the time range
 */
Matrix3D<double> readProcessedMatFiles3D(
    const std::string& name, 
    CustomDate timeStart, CustomDate timeEnd, const pmf::Parameters& parameters
); 

// Overloaded versions that accept double (MATLAB datenum) for time parameters

/**
 * @brief Read 1D data array from processed MAT files (double time version)
 * 
 * @param name Variable name to read (e.g., "MLT", "R0", "time")
 * @param timeStart Start time as MATLAB datenum
 * @param timeEnd End time as MATLAB datenum
 * @param parameters PMF parameters specifying data source and processing options
 * @return 1D matrix containing the requested data over the time range
 */
inline Matrix1D<double> readProcessedMatFiles1D(
    const std::string& name, 
    double timeStart, double timeEnd, const pmf::Parameters& parameters
) {
    return readProcessedMatFiles1D(name, CustomDate(timeStart), CustomDate(timeEnd), parameters);
}

/**
 * @brief Read 2D data array from processed MAT files (double time version)
 * 
 * @param name Variable name to read (e.g., "InvK" for K parameter data)
 * @param timeStart Start time as MATLAB datenum
 * @param timeEnd End time as MATLAB datenum
 * @param parameters PMF parameters specifying data source and processing options
 * @return 2D matrix containing the requested data over the time range
 */
inline Matrix2D<double> readProcessedMatFiles2D(
    const std::string& name, 
    double timeStart, double timeEnd, const pmf::Parameters& parameters
) {
    return readProcessedMatFiles2D(name, CustomDate(timeStart), CustomDate(timeEnd), parameters);
}

/**
 * @brief Read 3D data array from processed MAT files (double time version)
 * 
 * @param name Variable name to read (e.g., "PSD", "InvMu" for phase space density data)
 * @param timeStart Start time as MATLAB datenum
 * @param timeEnd End time as MATLAB datenum
 * @param parameters PMF parameters specifying data source and processing options
 * @return 3D matrix containing the requested data over the time range
 */
inline Matrix3D<double> readProcessedMatFiles3D(
    const std::string& name, 
    double timeStart, double timeEnd, const pmf::Parameters& parameters
) {
    return readProcessedMatFiles3D(name, CustomDate(timeStart), CustomDate(timeEnd), parameters);
}

} 