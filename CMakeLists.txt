# CMakeLists.txt
# Buildsystem for VERB4D - See README.md for build instructions
#
# Author: Anthony Miyaguchi (acmiyaguchi@gmail.com)

cmake_minimum_required(VERSION 3.9)
project(VERB4D C CXX)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Check if interprocedural optimization (link time optimization) is supported by the compiler
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)

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
    # Set the definition for MATLAB_CAPABLE
    add_definitions(-DMATLAB_CAPABLE=true)
endif(MATLAB_FOUND)

# Add directories for include files
include_directories(${PROJECT_SOURCE_DIR}/src)
# include_directories($ENV{HOME}/intel/oneapi/mkl/latest/include)

# Recursively search for header files to include for project files.  This isn't
# necessary to compile, but its nice to have in Visual Studio solutions.
file(GLOB_RECURSE VERB4D_HEADERS *.h)

# Add new source files here and re-run cmake
set(VERB4D_SOURCES
    src/Convection_1D_ULTIMATE_QUICKEST6.cpp
    src/Convection_2D.cpp
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
    src/VERB4D_Solver.cpp)

# the target executable
add_executable(VERB4D_Solver ${VERB4D_SOURCES} ${VERB4D_HEADERS})

# use interprocidural optimization if supported
if( ipo_supported )
    message(STATUS "IPO enabled")
    set_property(TARGET VERB4D_Solver PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO not supported: <${ipo_error}>")
endif()

if(WIN32)
    # The libraries for lapack and BLAS for windows are in the lib folder
    set(LAPACK_LIBRARIES ${CMAKE_SOURCE_DIR}/lib/clapack.lib)
    set(BLAS_LIBRARIES   ${CMAKE_SOURCE_DIR}/lib/blas.lib)
    set(F2C_LIBRARIES    ${CMAKE_SOURCE_DIR}/lib/libf2c.lib)
    # Libraries to link against
    # including Fortran to C library for windows, needed by lapack
    list(APPEND DEPLIBS
        ${LAPACK_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${F2C_LIBRARIES})
else(WIN32)
    # Set path to openblas library
    add_library(openblas SHARED IMPORTED)
    set_target_properties(openblas PROPERTIES IMPORTED_LOCATION $ENV{HOME}/.local/lib/libopenblas.so)
    set_target_properties(VERB4D_Solver PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

    set(MKL_THREADING sequential)
    set(MKL_INTERFACE lp64)
    find_package(MKL CONFIG REQUIRED)
    target_include_directories(VERB4D_Solver PUBLIC ${MKL_INCLUDE})
    
    list(APPEND DEPLIBS openblas ${MKL_IMPORTED_TARGETS})
    
    set_target_properties(VERB4D_Solver PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
    install(TARGETS VERB4D_Solver DESTINATION bin)
endif(WIN32)

target_link_libraries(VERB4D_Solver ${DEPLIBS})

# Flag compile type feedback
if(DEFINED CMAKE_BUILD_TYPE)
    if(CMAKE_BUILD_TYPE STREQUAL "")
        message(WARNING "\nCMAKE_BUILD_TPYE not specified. Run cmake with -DCMAKE_BUILD_TYPE=Release or -DCMAKE_BUILD_TYPE=Debug")
    else(CMAKE_BUILD_TYPE STREQUAL "")
        string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_upper)
        if(CMAKE_BUILD_TYPE_upper STREQUAL "RELEASE")
            add_definitions(-DRELEASE)
        endif()
    endif()
endif()
    
# Flag to decide convection step size in Convection_2D.cpp
if(${FAST_CONVECTION})
    add_definitions(-DFAST_CONVECTION)
endif()
unset(FAST_CONVECTION CACHE)
