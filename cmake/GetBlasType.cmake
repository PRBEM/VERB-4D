if(NOT DEFINED BLAS_TYPE)
    message("No BLAS_TYPE specified, choosing LAPACK")
    set(BLAS_TYPE LAPACK)
endif()

string(TOUPPER ${BLAS_TYPE} BLAS_TYPE)

# Initialize a local list for dependencies
set(BLAS_LIBRARIES_TO_LINK "")

if(${BLAS_TYPE} STREQUAL "LAPACK")

    if(WIN32)
        # The libraries for lapack and BLAS for windows are in the lib folder
        set(LAPACK_LIBRARIES ${CMAKE_SOURCE_DIR}/lib/clapack.lib)
        set(BLAS_LIBRARIES   ${CMAKE_SOURCE_DIR}/lib/blas.lib)
        set(F2C_LIBRARIES    ${CMAKE_SOURCE_DIR}/lib/libf2c.lib)
        list(APPEND BLAS_LIBRARIES_TO_LINK ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${F2C_LIBRARIES})
    else()
        # Installed by distribution's package manager
        find_package(LAPACK REQUIRED QUIET)
        list(APPEND BLAS_LIBRARIES_TO_LINK ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    endif()

elseif(${BLAS_TYPE} STREQUAL "OPENBLAS")
    if(NOT DEFINED OPENBLAS_LIBPATH)
        if(WIN32)
            set(OPENBLAS_LIBPATH "$ENV{HOMEPATH}/AppData/Local/OpenBLAS/openblas.lib")
        else()
            set(OPENBLAS_LIBPATH "$ENV{HOME}/.local/lib/libopenblas.so")
        endif()
        message(STATUS "No OPENBLAS_LIBPATH specified, using: ${OPENBLAS_LIBPATH}")
    endif()

    if(NOT TARGET OpenBLAS::OpenBLAS)
        add_library(OpenBLAS::OpenBLAS SHARED IMPORTED)
        if(WIN32)
            set_target_properties(OpenBLAS::OpenBLAS PROPERTIES IMPORTED_IMPLIB "${OPENBLAS_LIBPATH}")
        else()
            set_target_properties(OpenBLAS::OpenBLAS PROPERTIES IMPORTED_LOCATION "${OPENBLAS_LIBPATH}")
        endif()
    endif()
    list(APPEND BLAS_LIBRARIES_TO_LINK OpenBLAS::OpenBLAS)

elseif(${BLAS_TYPE} STREQUAL "MKL")
    set(MKL_THREADING sequential CACHE STRING "MKL threading mode")
    set(MKL_INTERFACE lp64 CACHE STRING "MKL interface mode")
    
    find_package(MKL CONFIG REQUIRED)
    list(APPEND BLAS_LIBRARIES_TO_LINK ${MKL_IMPORTED_TARGETS})
    
    # Store MKL include dir in a variable to be used after target creation
    set(BLAS_INCLUDE_DIRS ${MKL_INCLUDE})

else()
    message(FATAL_ERROR "Unknown BLAS_TYPE: ${BLAS_TYPE}. Choose LAPACK, OPENBLAS, or MKL")
endif()

# Export the list of libraries to the main scope
set(BLAS_LIBRARIES_TO_LINK ${BLAS_LIBRARIES_TO_LINK})