#include "Matrix.h"
#include "BoundaryConditionType.hpp"
#include <vector>
#include "mkl_sparse_qr.h"

typedef Matrix2D<double> mat2d;
typedef Matrix1D<double> mat1d;

void initialize_sparse_values(
    const mat2d& v_grid, const mat2d& k_grid,
    BoundaryConditionType v_lower_type, BoundaryConditionType v_upper_type,
    BoundaryConditionType k_lower_type, BoundaryConditionType k_upper_type,
    const mat2d& Dvv, const mat2d& Dvk, const mat2d& Dkv, const mat2d& Dkk, 
    const mat2d& jacobian, const mat2d& loss, double dt,
    std::vector<double>& values
);

void initialize_sparse_indices(int v_size, int k_size, std::vector<int>& column_indices, std::vector<int>& rows_csr);
void initialize_rhs(std::vector<double>& rhs, const Matrix2D<double>& psd, const mat1d& v_lower, const mat1d& v_upper, const mat1d& k_lower, const mat1d& k_upper, const mat2d& source, double dt);
void Diffusion_2D_MKL(
    mat2d& psd, const mat2d& v_grid, const mat2d& k_grid,
    BoundaryConditionType v_lower_type, BoundaryConditionType v_upper_type,
    BoundaryConditionType k_lower_type, BoundaryConditionType k_upper_type,
    const mat1d& v_lower, const mat1d& v_upper, const mat1d& k_lower, const mat1d& k_upper,
    const mat2d& Dvv, const mat2d& Dvk, const mat2d& Dkv, const mat2d& Dkk, 
    const mat2d& jacobian, const mat2d& loss, const mat2d& source, double dt,
    sparse_matrix_t* csrA
);
