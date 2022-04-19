# CMakeLists.txt
# Buildsystem for VERB4D - See README.md for build instructions
#
# Author: Anthony Miyaguchi (acmiyaguchi@gmail.com)

cmake_minimum_required(VERSION 2.8)
project(VERB4D C CXX)
add_definitions(-std=c++17)

# Enable/disable tests
include(CTest)
set(ENABLE_TESTING ON)

# Enable/disable data assimilation
set(DATA_ASSIMILATION ON CACHE BOOL "Enable data assimilation")

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Find OpenMP
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

# Find Matlab, Lapack, and BLAS
if(DATA_ASSIMILATION)
    find_package(MATLAB REQUIRED)
else(DATA_ASSIMILATION)
    find_package(MATLAB)
endif(DATA_ASSIMILATION)

if(MATLAB_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MATLAB_CXXFLAGS}")
    include_directories(${MATLAB_INCLUDE_DIR})
    list(APPEND DEPLIBS
        ${MATLAB_MAT_LIBRARY}
        ${MATLAB_MX_LIBRARY}
        ${MATLAB_ENG_LIBRARY})
    # Set the definition for MATLAB_CAPABLE
    add_definitions(-DMATLAB_CAPABLE=true)
endif(MATLAB_FOUND)

if(WIN32)
    # The libraries for lapack and BLAS for windows are in the lib folder
    set(LAPACK_LIBRARIES ${CMAKE_SOURCE_DIR}/lib/clapack.lib)
    set(BLAS_LIBRARIES   ${CMAKE_SOURCE_DIR}/lib/blas.lib)
    set(F2C_LIBRARIES    ${CMAKE_SOURCE_DIR}/lib/libf2c.lib)
else(WIN32)
    # Installed by distribution's package manager
    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)
endif(WIN32)

# Libraries to link against
list(APPEND DEPLIBS
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

# Fortran to C library for windows, needed by lapack
if(WIN32)
    list(APPEND DEPLIBS ${F2C_LIBRARIES})
endif(WIN32)

# Add directories for include files
include_directories(${PROJECT_SOURCE_DIR}/src)

# Recursively search for header files to include for project files.  This isn't
# necessary to compile, but its nice to have in Visual Studio solutions.
file(GLOB_RECURSE VERB4D_HEADERS *.h)

# Add new source files here and re-run cmake
set(VERB4D_SOURCES
    src/VERB4D_Solver.cpp)

set(VERB4D_SOURCES
    ${VERB4D_SOURCES}
    src/Convection_1D_ULTIMATE_QUICKEST6.cpp
    src/Convection_2D.cpp
    src/Diffusion_1D.cpp
    src/Diffusion_2D.cpp
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
    src/UpdatableMatrix.cpp)

# Add data assimilation source files and macros
if(DATA_ASSIMILATION)
    add_definitions(-DDATA_ASSIMILATION)

    set(VERB4D_SOURCES 
        ${VERB4D_SOURCES} 
        src/MatrixOperations.cpp
        src/PMF.cpp
        src/CustomDate.cpp
        src/DataAssimilationHelper.cpp
        src/DataAssimilation.cpp)
endif(DATA_ASSIMILATION)

# Build and link 
add_executable(VERB4D_Solver ${VERB4D_SOURCES} ${VERB4D_HEADERS})
target_link_libraries(VERB4D_Solver ${DEPLIBS})

# Add tests
if(ENABLE_TESTING)
    add_subdirectory(test)
endif()
if(FAST_CONVECTION)
    add_definitions(-DFAST_CONVECTION)
endif(FAST_CONVECTION)
