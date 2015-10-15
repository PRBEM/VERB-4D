cmake_minimum_required(VERSION 2.8)
project(VERB4D C CXX)
add_definitions(-std=c++11)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Find OpenMP
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

# Find Matlab, Lapack, and BLAS
find_package(MATLAB)
if(MATLAB_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MATLAB_CXXFLAGS}")
    include_directories(${MATLAB_INCLUDE_DIR})
    list(APPEND DEPLIBS
        ${MATLAB_MAT_LIBRARY}
        ${MATLAB_MX_LIBRARY}
        ${MATLAB_ENG_LIBRARY})
    add_definitions(-DMATLAB_CAPABLE=true)
endif(MATLAB_FOUND)

if(WIN32)
    set(LAPACK_LIBRARIES ${CMAKE_SOURCE_DIR}/lib/clapack.lib)
    set(BLAS_LIBRARIES ${CMAKE_SOURCE_DIR}/lib/BLAS.lib)
else(WIN32)
    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)
endif(WIN32)

# Libraries to link against
list(APPEND DEPLIBS
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

# Add directories for include files
include_directories(${PROJECT_SOURCE_DIR}/src)

# Set the source files to compile
set(VERB4D_SOURCES
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
    src/UpdatableMatrix.cpp
    src/VERB4D_Solver.cpp)

# Build and link
add_executable(VERB4D_Solver ${VERB4D_SOURCES})
target_link_libraries(VERB4D_Solver ${DEPLIBS})
