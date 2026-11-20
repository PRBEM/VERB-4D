include(FetchContent)

find_package(xtensor QUIET)

if(NOT xtensor_FOUND)
    # xtl/xtensor's pinned tags predate CMake 4's removal of pre-3.5 policy
    # support. Relax the minimum only for these fetched subprojects.
    set(CMAKE_POLICY_VERSION_MINIMUM 3.5)

    FetchContent_Declare(xtl GIT_REPOSITORY https://github.com/xtensor-stack/xtl.git GIT_TAG 0.7.6)
    FetchContent_MakeAvailable(xtl)

    FetchContent_Declare(
        xtensor
        GIT_REPOSITORY https://github.com/xtensor-stack/xtensor.git
        GIT_TAG 0.24.7
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
        )
    FetchContent_MakeAvailable(xtensor)

    unset(CMAKE_POLICY_VERSION_MINIMUM)
endif()