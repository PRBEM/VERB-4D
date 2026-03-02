# SPDX-FileCopyrightText: 2015 UCLA
# SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
#
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required(VERSION 3.14)
project(VERB4D LANGUAGES C CXX)

# --- Project Settings ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# --- Options ---
option(DATA_ASSIMILATION "Enable data assimilation" OFF)
option(DATA_ASSIMILATION_DEBUG "Enable debug output for data assimilation" OFF)
option(USE_PPFV "Compile positive preserving finite volume method" OFF)
option(FAST_CONVECTION "Use larger time steps for convection (less accurate)" OFF)
option(PYTHON_BINDINGS "Build Python bindings using pybind11" OFF)
option(SAVE_PSD_LOST_CONV "Save lost PSD during convection step" OFF)
option(LU_CACHING "Cache LU decompositions during Lapack calls" OFF)

# --- Dependencies ---
find_package(OpenMP)
include(GetBlasType) # this defines BLAS_LIBRARIES_TO_LINK and BLAS_INCLUDE_DIRS

if(USE_PPFV)
    message(STATUS "PPFV enabled: Finding or downloading Eigen3 and xtensor...")
    include(FindOrDownloadEigen)
    include(FindOrDownloadXTensor)
endif()

# --- Source Files ---
file(GLOB_RECURSE VERB4D_HEADERS src/*.h)

set(VERB4D_SOURCES
    src/VERB4D_Solver.cpp
    src/Convection_1D_ULTIMATE_QUICKEST6.cpp
    src/Convection_1D_2ndORDER_NONUNIFORM_GRID.cpp
    src/Convection_2D.cpp
    src/Convection_3D.cpp
    src/Diffusion_1D.cpp
    src/Diffusion_2D.cpp
    src/Diffusion_2D_MKL.cpp
    src/Diffusion_ADI1.cpp
    src/Diffusion_ADI2.cpp
    src/Diffusion_ADI3.cpp
    src/Interpolation.cpp
    src/Logger.cpp
    src/Matrix.cpp
    src/MatrixSolver.cpp
    src/MonotCubicInterpolator.cpp
    src/Parameters.cpp
    src/ReadInitialData.cpp
    src/UpdatableMatrix.cpp
)

if(DATA_ASSIMILATION)
    message(STATUS "Data Assimilation enabled")
    list(APPEND VERB4D_SOURCES 
        src/MatrixOperations.cpp
        src/CustomDate.cpp
        src/DataAssimilationHelper.cpp
        src/DataAssimilation.cpp
        src/DataAssimilationConvection.cpp
        src/DataAssimilationDiffusion.cpp
    )
endif()

if(USE_PPFV)
    list(APPEND VERB4D_SOURCES src/Diffusion_2D_PPFV.cpp)
endif()

# --- Main Executable Target ---
add_executable(VERB4D_Solver ${VERB4D_SOURCES} ${VERB4D_HEADERS})

# --- Compiler Definitions and Options ---
target_compile_definitions(VERB4D_Solver PRIVATE 
    _USE_MATH_DEFINES
    $<$<BOOL:${DATA_ASSIMILATION}>:DATA_ASSIMILATION>
    $<$<BOOL:${DATA_ASSIMILATION_DEBUG}>:DATA_ASSIMILATION_DEBUG>
    $<$<BOOL:${USE_PPFV}>:USE_PPFV>
    $<$<BOOL:${FAST_CONVECTION}>:FAST_CONVECTION>
    $<$<BOOL:${SAVE_PSD_LOST_CONV}>:SAVE_PSD_LOST_CONV>
    $<$<BOOL:${LU_CACHING}>:LU_CACHING>
    $<$<CONFIG:Release>:RELEASE>
)

if(MSVC)
    target_compile_options(VERB4D_Solver PRIVATE /W4)
else()
    target_compile_options(VERB4D_Solver PRIVATE -Wall -Wextra -Wpedantic)
endif()

# --- Include Directories ---
target_include_directories(VERB4D_Solver PUBLIC src)
if(BLAS_INCLUDE_DIRS)
    target_include_directories(VERB4D_Solver PRIVATE ${BLAS_INCLUDE_DIRS})
endif()

# --- Link Libraries ---
target_link_libraries(VERB4D_Solver PRIVATE ${BLAS_LIBRARIES_TO_LINK})

if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found! Parallelization enabled.")
    target_link_libraries(VERB4D_Solver PRIVATE OpenMP::OpenMP_CXX)
else()
    message(WARNING "OpenMP not found. VERB4D will be compiled in SERIAL mode.")
endif()

# --- PPFV ---
if(USE_PPFV)
    set(PPFV_SOURCES 
        src/ppfv/Mesh.cpp
        src/ppfv/Solver.cpp
        src/ppfv/Cases/VERB_ppfv2d.cpp
    )
    add_library(ppfv SHARED ${PPFV_SOURCES} ${VERB4D_HEADERS})
    
    target_link_libraries(ppfv PUBLIC Eigen3::Eigen)
    target_link_libraries(ppfv PUBLIC xtensor)

    # Link the internal library to the main solver
    target_link_libraries(VERB4D_Solver PRIVATE ppfv)
endif()

# --- Optimization (IPO/LTO) ---
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
if(ipo_supported)
    set_property(TARGET VERB4D_Solver PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    message(STATUS "IPO/LTO enabled for VERB4D_Solver")
endif()

# --- Python Bindings ---
if(PYTHON_BINDINGS)
    message(STATUS "Building Python bindings with pybind11...")
    add_subdirectory(extern/pybind11)

    pybind11_add_module(verb4d_solver src/python_bindings/Define_python_modules.cpp ${VERB4D_SOURCES})
    target_link_libraries(verb4d_solver PUBLIC OpenMP::OpenMP_CXX ${BLAS_LIBRARIES_TO_LINK})
    include_directories(${PROJECT_SOURCE_DIR}/src/python_bindings ${PROJECT_SOURCE_DIR}/src)

    if(USE_PPFV)
        target_compile_definitions(verb4d_solver PRIVATE USE_PPFV)
        target_link_libraries(verb4d_solver PUBLIC ppfv)
    endif()
endif()