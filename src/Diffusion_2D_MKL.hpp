#include "Matrix.h"
#include "BoundaryConditionType.hpp"
#include <vector>
typedef Matrix2D<double> mat2d;
void initialize_sparse_values(
    const mat2d& v_grid, const mat2d& k_grid,
    BoundaryConditionType v_lower, BoundaryConditionType v_upper,
    BoundaryConditionType k_lower, BoundaryConditionType k_upper,
    const mat2d& Dvv, const mat2d& Dvk, const mat2d& Dkv, const mat2d& Dkk, 
    const mat2d& jacobian, const mat2d& loss, double dt,
    std::vector<double>& values
);

void initialize_sparse_indices(int v_size, int k_size, std::vector<int>& column_indices, std::vector<int>& rows_csr);

void Diffusion_2D_MKL(
    mat2d& psd, const mat2d& v_grid, const mat2d& k_grid,
    BoundaryConditionType v_lower, BoundaryConditionType v_upper,
    BoundaryConditionType k_lower, BoundaryConditionType k_upper,
    const mat2d& Dvv, const mat2d& Dvk, const mat2d& Dkv, const mat2d& Dkk, 
    const mat2d& jacobian, const mat2d& loss, double dt,
    std::vector<int>& column_indices, std::vector<int>& rows_csr
);