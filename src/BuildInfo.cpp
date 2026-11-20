/*
 * SPDX-FileCopyrightText: 2026 GFZ Helmholtz Centre Potsdam
 * 
 * SPDX-FileContributor: Sahil Jhawar
 *
 * SPDX-License-Identifier: BSD-3
 */

/**
 * \file BuildInfo.cpp
 * \brief Implementation of the --version/--help/--check command line diagnostics.
 */

#include "BuildInfo.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "GitVersion.h"
#include "Logger.h"
#include "MatrixSolver.h"

#ifndef BLAS_TYPE_STRING
#define BLAS_TYPE_STRING "unknown"
#endif

void printBuildInfo(std::ostream &out)
{
    std::string compile_type = "Debug";
#ifdef RELEASE
    compile_type = "Release";
#endif
    out << "VERB4D_Solver (" << compile_type << ") compiled on: " << __DATE__ << ", " << __TIME__ << std::endl;
    std::string git_commit = GIT_COMMIT_HASH;
    if (git_commit.empty() || git_commit == "unknown") {
        out << "Git commit:    unknown (maybe a container build)" << std::endl;
    } else {
        out << "Git commit:    " << git_commit << std::endl;
    }
    out << "BLAS backend:  " << BLAS_TYPE_STRING << std::endl;
#ifdef _OPENMP
    out << "OpenMP:        enabled (max threads: " << omp_get_max_threads() << ")" << std::endl;
#else
    out << "OpenMP:        disabled (serial build)" << std::endl;
#endif
    std::vector<std::string> features;
#ifdef DATA_ASSIMILATION
    features.push_back("DATA_ASSIMILATION");
#endif
#ifdef DATA_ASSIMILATION_DEBUG
    features.push_back("DATA_ASSIMILATION_DEBUG");
#endif
#ifdef USE_PPFV
    features.push_back("USE_PPFV");
#endif
#ifdef FAST_CONVECTION
    features.push_back("FAST_CONVECTION");
#endif
#ifdef SAVE_PSD_LOST_CONV
    features.push_back("SAVE_PSD_LOST_CONV");
#endif
#ifdef LU_CACHING
    features.push_back("LU_CACHING");
#endif
    out << "Build options: ";
    if (features.empty())
    {
        out << "(none)";
    }
    else
    {
        for (size_t i = 0; i < features.size(); i++)
        {
            out << features[i] << (i + 1 < features.size() ? ", " : "");
        }
    }
    out << std::endl;
}

void printUsage(const char *argv0)
{
    std::cout << "Usage: " << argv0 << " [parameter=value ...]" << std::endl
              << std::endl;
    std::cout << "  -h, --help     Show this help and exit." << std::endl;
    std::cout << "  --version      Show build/version information and exit." << std::endl;
    std::cout << "  --check        Run an installation self-check (BLAS/LAPACK, OpenMP) and exit." << std::endl;
    std::cout << "                 Does not require parameters.ini or any input files." << std::endl
              << std::endl;
    std::cout << "Normal operation requires a parameters.ini file" << std::endl;
}

bool runSelfCheck()
{
    std::cout << "VERB4D_Solver installation self-check" << std::endl;
    std::cout << LOGGER_SEPARATOR << std::endl;
    printBuildInfo(std::cout);
    std::cout << LOGGER_SEPARATOR << std::endl;

    bool ok = true;

#ifdef _OPENMP
    std::cout << "[ OK ] OpenMP available, max threads = " << omp_get_max_threads() << std::endl;
#else
    std::cout << "[WARN] Compiled without OpenMP; solver will run single-threaded." << std::endl;
#endif

    // Solve the tridiagonal system A*x = b, A = [[2,1,0],[1,3,1],[0,1,2]], x = [1,1,1],
    // using the same LAPACK banded routines (dgbtrf_/dgbtrs_) MatrixSolver.cpp uses for
    // local diffusion.
    {
        constexpr long n = 3, kl = 1, ku = 1, ldab = 2 * kl + ku + 1;
        double ab[ldab * n] = {0};
        long ipiv[n] = {0};

        auto set = [&](long row, long col, double value)
        {
            long idx = col - row;
            ab[col * ldab + (kl + ku - idx)] = value;
        };
        set(0, 0, 2);
        set(0, 1, 1);
        set(1, 0, 1);
        set(1, 1, 3);
        set(1, 2, 1);
        set(2, 1, 1);
        set(2, 2, 2);

        double b[n] = {3, 5, 3}; // = A * [1,1,1]
        long m = n, kl_ = kl, ku_ = ku, ldab_ = ldab, nrhs = 1, ldb = n, info = 0;
        const char trans = 'N';

        dgbtrf_(&m, &m, &kl_, &ku_, ab, &ldab_, ipiv, &info);
        if (info != 0)
        {
            std::cout << "[FAIL] LAPACK dgbtrf_ (LU factorization) returned info=" << info << std::endl;
            ok = false;
        }
        else
        {
            dgbtrs_(&trans, &m, &kl_, &ku_, &nrhs, ab, &ldab_, ipiv, b, &ldb, &info);
            if (info != 0)
            {
                std::cout << "[FAIL] LAPACK dgbtrs_ (solve) returned info=" << info << std::endl;
                ok = false;
            }
            else
            {
                double max_error = 0;
                for (long i = 0; i < n; i++)
                {
                    max_error = std::max(max_error, std::abs(b[i] - 1.0));
                }
                if (max_error < 1e-9)
                {
                    std::cout << "[ OK ] BLAS/LAPACK (" << BLAS_TYPE_STRING << ") solved a known tridiagonal system correctly." << std::endl;
                }
                else
                {
                    std::cout << "[FAIL] BLAS/LAPACK (" << BLAS_TYPE_STRING << ") returned a wrong result (max error "
                              << max_error << "). The linked library may be incompatible "
                                              "(e.g. mismatched integer width)."
                              << std::endl;
                    ok = false;
                }
            }
        }
    }

    std::cout << LOGGER_SEPARATOR << std::endl;
    std::cout << (ok ? "Self-check PASSED." : "Self-check FAILED.") << std::endl;
    return ok;
}
