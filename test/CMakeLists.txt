# This file is intended for building unit tests

set(SRC_FOLDER "../src/")

# Tests for data assimilation
if(DATA_ASSIMILATION)
    add_executable(test_data_assimilation 
                    ./TestDataAssimilation.cpp
                    ./MatrixOperationsTests.cpp
                    ./DataAssimilationInternalTests.cpp
                    ./PMFInternalTests.cpp
                    ${SRC_FOLDER}/MatrixOperations.cpp
                    ${SRC_FOLDER}/Matrix.cpp
                    ${SRC_FOLDER}/DataAssimilationHelper.cpp
                    ${SRC_FOLDER}/DataAssimilation.cpp 
                    ${SRC_FOLDER}/MatrixOperations.cpp 
                    ${SRC_FOLDER}/CustomDate.cpp 
                    ${SRC_FOLDER}/Parameters.cpp 
                    ${SRC_FOLDER}/PMF.cpp 
                    ${SRC_FOLDER}/Logger.cpp 
                    ${SRC_FOLDER}/Matrix.cpp)

    target_include_directories(test_data_assimilation PRIVATE ./ ${SRC_FOLDER})
    target_compile_options(test_data_assimilation PRIVATE -std=c++11 -g)
    target_link_libraries(test_data_assimilation ${DEPLIBS})
    add_test(NAME TestDataAssimilation COMMAND test_data_assimilation)
endif()
