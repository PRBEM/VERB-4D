# Regenerates GitVersion.h from GitVersion.h.in with the current git commit hash.
# Run at build time (not just configure time) as a custom target, so the hash stays
# fresh as commits are made without needing to rerun cmake. Expects SRC, DST, SRC_DIR
# to be passed in via -D.

if(DEFINED OVERRIDE AND NOT OVERRIDE STREQUAL "")
    set(GIT_COMMIT_HASH "${OVERRIDE}")
else()
    find_package(Git QUIET)

    set(GIT_COMMIT_HASH "unknown")

    if(GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
            WORKING_DIRECTORY ${SRC_DIR}
            OUTPUT_VARIABLE GIT_COMMIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
            RESULT_VARIABLE GIT_RESULT
        )
        if(NOT GIT_RESULT EQUAL 0 OR GIT_COMMIT_HASH STREQUAL "")
            set(GIT_COMMIT_HASH "unknown")
        else()
            execute_process(
                COMMAND ${GIT_EXECUTABLE} diff --quiet HEAD --
                WORKING_DIRECTORY ${SRC_DIR}
                RESULT_VARIABLE GIT_DIRTY
            )
            if(NOT GIT_DIRTY EQUAL 0)
                set(GIT_COMMIT_HASH "${GIT_COMMIT_HASH}-dirty")
            endif()
        endif()
    endif()
endif()

configure_file(${SRC} ${DST} @ONLY)
