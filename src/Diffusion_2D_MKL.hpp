#pragma once
#ifdef MKL_FOUND
#include "Matrix.h"
#include "BoundaryConditionType.hpp"
#include <vector>
#include "mkl_sparse_qr.h"

void initialize_sparse_values(
    const Matrix2D<double>& v_grid, const Matrix2D<double>& k_grid,
    BoundaryConditionType v_lower_type, BoundaryConditionType v_upper_type,
    BoundaryConditionType k_lower_type, BoundaryConditionType k_upper_type,
    const Matrix2D<double>& Dvv, const Matrix2D<double>& Dvk, const Matrix2D<double>& Dkv, const Matrix2D<double>& Dkk, 
    const Matrix2D<double>& jacobian, const Matrix2D<double>& loss, double dt,
    std::vector<double>& values
);

void initialize_sparse_indices(int v_size, int k_size, std::vector<int>& column_indices, std::vector<int>& rows_csr);
void initialize_rhs(
    std::vector<double>& rhs, const Matrix2D<double>& psd, const Matrix1D<double>& v_lower,
    const Matrix1D<double>& v_upper, const Matrix1D<double>& k_lower, const Matrix1D<double>& k_upper,
    const Matrix2D<double>& source, double dt
);
void Diffusion_2D_MKL(
    Matrix2D<double>& psd, const Matrix2D<double>& v_grid, const Matrix2D<double>& k_grid,
    BoundaryConditionType v_lower_type, BoundaryConditionType v_upper_type,
    BoundaryConditionType k_lower_type, BoundaryConditionType k_upper_type,
    const Matrix1D<double>& v_lower, const Matrix1D<double>& v_upper, const Matrix1D<double>& k_lower, const Matrix1D<double>& k_upper,
    const Matrix2D<double>& Dvv, const Matrix2D<double>& Dvk, const Matrix2D<double>& Dkv, const Matrix2D<double>& Dkk, 
    const Matrix2D<double>& jacobian, const Matrix2D<double>& loss, const Matrix2D<double>& source, double dt,
    sparse_matrix_t* csrA
);
#endif