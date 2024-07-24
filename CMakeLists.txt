# CMakeLists.txt
# Buildsystem for VERB4D - See README.md for build instructions
#
# Author: Anthony Miyaguchi (acmiyaguchi@gmail.com)

cmake_minimum_required(VERSION 3.9)
project(VERB4D C CXX)
set(CMAKE_CXX_STANDARD 17)

# Enable/disable data assimilation
set(DATA_ASSIMILATION ON CACHE BOOL "Enable data assimilation")

# Math Constants are not defined in Standard C++. To use them, we must first define _USE_MATH_DEFINES
add_compile_definitions(_USE_MATH_DEFINES)

# Recursively search for header files to include for project files.  This isn't
# necessary to compile, but its nice to have in Visual Studio solutions.
file(GLOB_RECURSE VERB4D_HEADERS *.h)

# Add new source files here and re-run cmake
set(VERB4D_SOURCES
    src/VERB4D_Solver.cpp
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
    
    if(DATA_ASSIMILATION_DEBUG)
        message("Saving DA debug output!")
        add_definitions(-DDATA_ASSIMILATION_DEBUG)
    endif(DATA_ASSIMILATION_DEBUG)

endif(DATA_ASSIMILATION)

# Build and link 
add_executable(VERB4D_Solver ${VERB4D_SOURCES} ${VERB4D_HEADERS})
target_compile_options(VERB4D_Solver PRIVATE -Wall -Wextra -Wpedantic -Werror)

# Add directories for include files
target_include_directories(VERB4D_Solver PUBLIC ${PROJECT_SOURCE_DIR}/src)


# Check if interprocedural optimization (link time optimization) is supported by the compiler
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)

if(NOT DEFINED BLAS_TYPE)
    message("No BLAS_TYPE specified, choosing LAPACK")
    set(BLAS_TYPE LAPACK)
endif()
string(TOUPPER ${BLAS_TYPE} BLAS_TYPE)
if(${BLAS_TYPE} STREQUAL "LAPACK")
    if(WIN32)
        # The libraries for lapack and BLAS for windows are in the lib folder
        set(LAPACK_LIBRARIES ${CMAKE_SOURCE_DIR}/lib/clapack.lib)
        set(BLAS_LIBRARIES   ${CMAKE_SOURCE_DIR}/lib/blas.lib)
        set(F2C_LIBRARIES    ${CMAKE_SOURCE_DIR}/lib/libf2c.lib)
        list(APPEND DEPLIBS ${F2C_LIBRARIES})
    else(WIN32)
        # Installed by distribution's package manager
        find_package(LAPACK REQUIRED)
    endif(WIN32)
    list(APPEND DEPLIBS ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
elseif(${BLAS_TYPE} STREQUAL "OPENBLAS")
    add_library(openblas  SHARED IMPORTED)
    if(WIN32)
        if(NOT DEFINED OPENBLAS_LIBPATH)
            set(OPENBLAS_LIBPATH $ENV{HOMEPATH}/AppData/Local/OpenBLAS/openblas.lib)
            message("No OPENBLAS_LIBPATH specified, choosing ${OPENBLAS_LIBPATH}")
        endif()
        set_target_properties(openblas PROPERTIES IMPORTED_IMPLIB ${OPENBLAS_LIBPATH})
    else()
        if(NOT DEFINED OPENBLAS_LIBPATH)
            set(OPENBLAS_LIBPATH $ENV{HOME}/.local/lib/libopenblas.so)
            message("No OPENBLAS_LIBPATH specified, choosing ${OPENBLAS_LIBPATH}")
        endif()
        set_target_properties(openblas PROPERTIES IMPORTED_LOCATION ${OPENBLAS_LIBPATH})
    endif()
    list(APPEND DEPLIBS openblas)
elseif(${BLAS_TYPE} STREQUAL "MKL")
    set(MKL_THREADING sequential)
    set(MKL_INTERFACE lp64)
    find_package(MKL CONFIG REQUIRED)
    list(APPEND DEPLIBS ${MKL_IMPORTED_TARGETS})
    target_include_directories(VERB4D_Solver PUBLIC ${MKL_INCLUDE})
    add_definitions(-DMKL_FOUND)
else()
    message( FATAL_ERROR "Unknown BLAS_TYPE. Choose either LAPACK, OPENBLAS or MKL" )
endif()

# Find OpenMP
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

# Find Matlab
if(DATA_ASSIMILATION)
    find_package(Matlab REQUIRED COMPONENTS ENG_LIBRARY MAT_LIBRARY)
else(DATA_ASSIMILATION)
    find_package(Matlab COMPONENTS ENG_LIBRARY MAT_LIBRARY)
endif(DATA_ASSIMILATION)

if(Matlab_FOUND)
    include_directories(${Matlab_INCLUDE_DIRS})
    list(APPEND DEPLIBS
        ${Matlab_MAT_LIBRARY}
        ${Matlab_MX_LIBRARY}
        ${Matlab_ENG_LIBRARY})
    # Set the definition for MATLAB_CAPABLE
    add_definitions(-DMATLAB_CAPABLE=true)
endif()

# use interprocidural optimization if supported
if( ipo_supported )
    message(STATUS "IPO enabled")
    set_property(TARGET VERB4D_Solver PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO not supported: <${ipo_error}>")
endif()

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

# Python bindings
if (${PYTHON_BINDINGS})
    include_directories(${PROJECT_SOURCE_DIR}/extern/pybind11/include ${PROJECT_SOURCE_DIR}/src/python_bindings ${PROJECT_SOURCE_DIR}/src)
    add_subdirectory(extern/pybind11)
    pybind11_add_module(verb4d_solver src/python_bindings/Define_python_modules.cpp ${VERB4D_SOURCES})
    target_link_libraries(verb4d_solver PUBLIC ${DEPLIBS})
endif()
unset(PYTHON_BINDINGS CACHE)
