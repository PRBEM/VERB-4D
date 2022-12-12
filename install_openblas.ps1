# create a build folder and compile OpenBLAS
Set-Location OpenBLAS
if (!(Test-Path build)) {
    New-Item build -ItemType Directory
}
Set-Location build
cmake .. -DBUILD_STATIC_LIBS=FALSE -DBUILD_SHARED_LIBS=TRUE -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=TRUE -DUSE_THREAD=0 -DUSE_LOCKING=1 -DUSE_OPENMP=1
cmake --build . --target openblas_shared -j --config Release

# copy the OpenBLAS library files to the user's AppData folder and append
# that location to the user's PATH environment variable
$blas_dir = $HOME + "\AppData\Local\OpenBLAS"
if (!(Test-Path $blas_dir)) {
    New-Item $blas_dir -ItemType Directory
}

Copy-Item lib\RELEASE\openblas.lib $blas_dir
Copy-Item lib\RELEASE\openblas.dll $blas_dir
Set-Location ..\..

$oldPath = (Get-ItemProperty -Path 'Registry::HKEY_CURRENT_USER\Environment' -Name PATH).path
if (!$oldPath.split(';').contains($blas_dir)){
    $newPath = ($oldPath + $blas_dir) -join ';'
    Set-ItemProperty -Path 'Registry::HKEY_CURRENT_USER\Environment' -Name PATH -Value $newPath
}