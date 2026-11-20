/*
 * SPDX-FileCopyrightText: 2026 GFZ Helmholtz Centre Potsdam
 *
 * SPDX-FileContributor: Sahil Jhawar
 *
 * SPDX-License-Identifier: BSD-3
 */

/**
 * \file BuildInfo.h
 * \brief Version/usage banner and installation self-check, exposed as --version,
 * --help and --check on the VERB4D_Solver command line.
 */

#pragma once
#include <ostream>

void printBuildInfo(std::ostream &out);
void printUsage(const char *argv0);
bool runSelfCheck();
