# create a build folder and compile OpenBLAS
$openblas_branch = "release-0.3.0"

Set-Location OpenBLAS
git checkout $openblas_branch
if (!(Test-Path build)) {
    New-Item build -ItemType Directory
}
Set-Location build
cmake .. -DBUILD_STATIC_LIBS=FALSE -DBUILD_SHARED_LIBS=TRUE -DUSE_THREAD=0 -DUSE_LOCKING=1 -DUSE_OPENMP=1 -DC_LAPACK=ON
cmake --build . --target openblas_shared -j --config Release

# copy the OpenBLAS library files to the user's local AppData folder
# and append that location to the user's PATH environment variable
$blas_dir = Join-Path -Path $env:LOCALAPPDATA -ChildPath "OpenBLAS"
if (!(Test-Path $blas_dir)) {
    New-Item $blas_dir -ItemType Directory
}

Copy-Item lib\RELEASE\openblas.lib $blas_dir
Copy-Item lib\RELEASE\openblas.dll $blas_dir
Write-Output "copied OpenBLAS library files to $blas_dir"
Set-Location ..\..

$oldPath = (Get-ItemProperty -Path 'Registry::HKEY_CURRENT_USER\Environment' -Name PATH).path
if (!$oldPath.split(';').contains($blas_dir)){
    $newPath = ($oldPath + $blas_dir) -join ';'
    Set-ItemProperty -Path 'Registry::HKEY_CURRENT_USER\Environment' -Name PATH -Value $newPath
}