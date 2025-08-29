/*
 * SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * @file DataAssimilationHelper.h
 * @brief Helper functions for data assimilation operations
 * 
 * This file contains utility functions for data assimilation including model matrix
 * generation, data processing, interpolation, and observation space conversions.
 */

#pragma once
#include "DataAssimilation.h"
#include "Matrix.h"
#include "PMF.h"
#include "Parameters.h"
#include <string>
#include <vector>

extern const double FILLVAL;
using ParametersIni=::Parameters;

/**
 * @namespace data_assimilation::internal
 * @brief Internal helper functions for data assimilation operations
 */
namespace data_assimilation::internal {
    /**
     * @brief Generate 2D convection model matrix with stability checking
     * 
     * @param VP Velocity field in P direction
     * @param VR Velocity field in R direction  
     * @param Loss Loss term matrix
     * @param saturation_map Boolean map indicating saturation regions
     * @param timeStep Time step for integration
     * @param dP Grid spacing in P direction
     * @param dR Grid spacing in R direction
     * @param use_log Whether to use logarithmic formulation
     * @param correlation_time Correlation time for boundary conditions
     * @return Model matrix for 2D convection with automatic time step splitting
     */
    Matrix2D<double> getModelMatrixConvection2D(
        const Matrix2D<double> &VP,
        const Matrix2D<double> &VR,
        const Matrix2D<double> &Loss,
        const Matrix2D<bool> &saturation_map,
        double timeStep,
        double dP,
        double dR,
        bool use_log,
        double correlation_time);
    /**
     * @brief Generate 2D convection model matrix without stability checking
     * 
     * @param VP Velocity field in P direction
     * @param VR Velocity field in R direction
     * @param Loss Loss term matrix
     * @param saturation_map Boolean map indicating saturation regions
     * @param timeStep Time step for integration
     * @param dP Grid spacing in P direction
     * @param dR Grid spacing in R direction
     * @param use_log Whether to use logarithmic formulation
     * @param correlation_factor Correlation factor for boundary conditions
     * @return Model matrix for single time step without stability checks
     */
    Matrix2D<double> getModelMatrixConvection2DNoStabilityCheck(
        const Matrix2D<double> &VP,
        const Matrix2D<double> &VR,
        const Matrix2D<double> &Loss,
        const Matrix2D<bool> &saturation_map,
        double timeStep,
        double dP,
        double dR,
        bool use_log,
        double correlation_factor);
    /**
     * @brief Split time step based on Courant-Friedrichs-Lewy (CFL) condition
     * 
     * @param maximumCourantNumber Maximum allowed Courant number for stability
     * @param timeStep Original time step
     * @param VP Velocity field in P direction
     * @param VR Velocity field in R direction
     * @param dP Grid spacing in P direction
     * @param dR Grid spacing in R direction
     * @return Pair containing number of sub-steps and adjusted time step size
     */
    std::pair<int, double> splitTimeStepCourantCondition(
        double maximumCourantNumber, double timeStep, 
        const Matrix2D<double>& VP, const Matrix2D<double>& VR, 
        double dP, double dR
    );
    /**
     * @brief Generate 2D convection model matrix (simplified version)
     * 
     * @param VP Velocity field in P direction
     * @param VR Velocity field in R direction
     * @param Loss Loss term matrix
     * @param timeStep Time step for integration
     * @param dP Grid spacing in P direction
     * @param dR Grid spacing in R direction
     * @return Model matrix for 2D convection without saturation mapping
     */
    Matrix2D<double> getModelMatrixConvection2D(
        const Matrix2D<double>& VP,
        const Matrix2D<double>& VR,
        const Matrix2D<double>& Loss,
        double timeStep,
        double dP,
        double dR
    );

    /**
     * @brief Calculate conditional mean of values where flag is true
     * 
     * @param a Array of values to average
     * @param flag Boolean array indicating which values to include
     * @return Mean of values where flag is true, NAN if no valid values
     */
    double conditionalMean(const Matrix1D<double>& a, const Matrix1D<bool>& flag);
    
    /**
     * @brief Bin observations into P-R grid cells
     * 
     * @param observations 2D vector of observation data
     * @param P Grid coordinates in P direction
     * @param R Grid coordinates in R direction
     * @param type Binning type ("linear" or "log10")
     * @return 4D matrix with binned PSD values [P][R][V][K]
     */
    Matrix4D<double> bin(
        const std::vector<std::vector<data_assimilation::Observations>> &observations,
        const Matrix2D<double> &P, const Matrix2D<double> &R,
        const std::string &type = "linear");

    /**
     * @brief Convert 2D data matrix to observation space format
     * 
     * @param data 2D matrix of observation data
     * @return ObservationSpace structure with data vector and observation operator H
     */
    data_assimilation::ObservationSpace convertToObservationSpace(const Matrix2D<double>& data);

    /**
     * @brief Read and process satellite data from MAT files
     * 
     * @param timeStart Start time for data reading
     * @param timeEnd End time for data reading
     * @param parameters PMF parameters specifying satellite and instrument
     * @return Processed data structure with MLT, R, K, V, and PSD arrays
     */
    data_assimilation::ProcessedMatFileData readData(double timeStart, double timeEnd, const pmf::Parameters& parameters);

    /**
     * @brief Concatenate multiple processed MAT file data structures
     * 
     * @param pmfDataSplit Vector of processed data structures to concatenate
     * @return Single concatenated data structure
     */
    data_assimilation::ProcessedMatFileData cat(const std::vector<ProcessedMatFileData>& pmfDataSplit);

    /**
     * @brief Perform 2D bilinear interpolation using four corner points
     * 
     * @param V_in Input V coordinate grid
     * @param K_in Input K coordinate grid
     * @param PSD_in Input PSD values to interpolate
     * @param V_out Target V coordinate for interpolation
     * @param K_out Target K coordinate for interpolation
     * @return Interpolated PSD value, NAN if interpolation not possible
     */
    double interp2d_four_corners(
        const Matrix2D<double>& V_in, const Matrix2D<double>& K_in,
        const Matrix2D<double>& PSD_in, double V_out, double K_out
    );

    /**
     * @brief Interpolate processed data onto specified V-K grid
     * 
     * @param data Vector of processed MAT file data from multiple instruments
     * @param V_grid Target velocity grid for interpolation
     * @param K_grid Target K parameter grid for interpolation
     * @return 2D vector of observations interpolated onto the target grid
     */
    std::vector<std::vector<data_assimilation::Observations>> interpolate(
        const std::vector<data_assimilation::ProcessedMatFileData>& data,
        const Matrix2D<double>& V_grid,
        const Matrix2D<double>& K_grid 
    );
}
