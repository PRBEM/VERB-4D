#!/bin/bash

# SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
#
# SPDX-License-Identifier: BSD-3-Clause

#openblas_branch="release-0.3.0"

# create ~/.local/lib if it does not exist
local_lib="$HOME/.local/lib"
if [ ! -d $local_lib ]
then
	echo "creating folder $local_lib"
	mkdir -p $local_lib
fi

# create a build folder inside OpenBLAS if it does not exist
cd ./OpenBLAS
#git checkout $openblas_branch
if [ ! -d ./build ]
then
    mkdir ./build
fi

cd ./build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF -DUSE_THREAD=0 -DUSE_LOCKING=1 -DUSE_OPENMP=1
make -j
cp ./lib/libopenblas.so* $local_lib
cd ..
