# SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
#
# SPDX-License-Identifier: BSD-3-Clause

# create ./build folder if it does not exist
if [ ! -d ./build ]
then
    mkdir ./build
fi

# run cmake and make to compile VERB4D_Solver
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBLAS_TYPE=OpenBLAS
make -j
cd ..
