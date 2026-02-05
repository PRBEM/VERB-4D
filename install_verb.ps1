# SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
#
# SPDX-License-Identifier: BSD-3-Clause

# create a build folder and compile the VERB-4D solver
if (!(Test-Path build)) {
    New-Item build -ItemType Directory
}
Set-Location build
cmake .. -DBLAS_TYPE=OpenBLAS
cmake --build . --target VERB4D_Solver -j --config Release
Set-Location ..