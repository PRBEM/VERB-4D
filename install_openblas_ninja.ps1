# SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
#
# SPDX-License-Identifier: BSD-3-Clause

# check if Ninja is installed
try {
    ninja --version
}
catch {
    throw "Ninja build system not found on PATH. Download Ninja binary from https://github.com/ninja-build/ninja/releases and add the location of ninja.exe to your PATH environmental variable."
}
# create a build folder and compile OpenBLAS
Set-Location OpenBLAS
if (!(Test-Path build)) {
    New-Item build -ItemType Directory
}
Set-Location build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_STATIC_LIBS=FALSE -DBUILD_SHARED_LIBS=TRUE -DUSE_THREAD=0 -DUSE_LOCKING=1 -DUSE_OPENMP=1 -GNinja
cmake --build . --target openblas_shared -j --config Release

# copy the OpenBLAS library files to the user's local AppData folder
# and append that location to the user's PATH environment variable
$blas_dir = Join-Path -Path $env:LOCALAPPDATA -ChildPath "OpenBLAS"
if (!(Test-Path $blas_dir)) {
    New-Item $blas_dir -ItemType Directory
}

Copy-Item lib\libopenblas.dll.a $blas_dir\openblas.lib
Copy-Item lib\libopenblas.dll $blas_dir\openblas.dll
Write-Output "copied OpenBLAS library files to $blas_dir"
Set-Location ..\..

$oldPath = (Get-ItemProperty -Path 'Registry::HKEY_CURRENT_USER\Environment' -Name PATH).path
[System.Collections.ArrayList]$path_list = $oldPath.split(';')
if (!$path_list.contains($blas_dir)){
    $Env:Path += ';' + $blas_dir
    [void]$path_list.Add($blas_dir)
    $newPath = $path_list -join ';'
    Set-ItemProperty -Path 'Registry::HKEY_CURRENT_USER\Environment' -Name PATH -Value $newPath
}
