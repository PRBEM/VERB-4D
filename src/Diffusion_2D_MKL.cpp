// SPDX-FileCopyrightText: 2015 UCLA
// SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef MKL_FOUND
#include "Matrix.h"
#include "MatrixSolver.h"
#include "BoundaryConditionType.hpp"
#include "Logger.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include "mkl_sparse_qr.h"
#include "Diffusion_2D_MKL.hpp"

void initialize_sparse_values(
        const Matrix2D<double>& v_grid, const Matrix2D<double>& k_grid,
        BoundaryConditionType v_lower_type, BoundaryConditionType v_upper_type,
        BoundaryConditionType k_lower_type, BoundaryConditionType k_upper_type,
        const Matrix2D<double>& Dvv, const Matrix2D<double>& Dvk, const Matrix2D<double>& Dkv, const Matrix2D<double>& Dkk,
        const Matrix2D<double>& jacobian, const Matrix2D<double>& loss, double dt,
        std::vector<double>& values)
{
    int v_size = v_grid.size_q1;
    int k_size = k_grid.size_q2;
    values.clear();
    values.reserve(
        2 * (2 * (v_size - 2)) // two blocks in the beginning and two in the end, coresponding to lower and upper K boundary
        + (k_size - 2) * (v_size + 2 * (v_size - 1)) // inner rows, main diag block with one main diagonal and two off-diagonals
        + 2 * (k_size - 2) * (3 * (v_size - 2))   // off-diag block with 3 diags of length v_size - 2
    );
    double v_offdiag_lower, v_offdiag_upper;
    double k_offdiag_lower, k_offdiag_upper;
    switch(v_lower_type)
    {
        case BoundaryConditionType::ConstantDerivative:
            v_offdiag_lower = 1;
            break;
        case BoundaryConditionType::ConstantValue:
            v_offdiag_lower = 0;
            break;
        default:
            Logger::error << "Lower diffusion boundary condition in V must be either constant value or constant derivative\n";
            exit(EXIT_FAILURE);
    }
    switch(v_upper_type)
    {
        case BoundaryConditionType::ConstantDerivative:
            v_offdiag_upper = 1;
            break;
        case BoundaryConditionType::ConstantValue:
            v_offdiag_upper = 0;
            break;
        default:
            Logger::error << "Upper diffusion boundary condition in V must be either constant value or constant derivative\n";
            exit(EXIT_FAILURE);
    }
    switch(k_lower_type)
    {
        case BoundaryConditionType::ConstantDerivative:
            k_offdiag_lower = 1;
            break;
        case BoundaryConditionType::ConstantValue:
            k_offdiag_lower = 0;
            break;
        default:
            Logger::error << "Lower diffusion boundary condition in K must be either constant value or constant derivative\n";
            exit(EXIT_FAILURE);
    }
    switch(k_upper_type)
    {
        case BoundaryConditionType::ConstantDerivative:
            k_offdiag_upper = 1;
            break;
        case BoundaryConditionType::ConstantValue:
            k_offdiag_upper = 0;
            break;
        default:
            Logger::error << "Upper diffusion boundary condition in K must be either constant value or constant derivative\n";
            exit(EXIT_FAILURE);
    }

    // lower K boundary, lower V boundary
    // in this edge case, the V boundary condition its prioritized
    values.push_back(1);
    values.push_back(-v_offdiag_lower);

    // lower K boundary, inner V
    for(int i = 1; i < v_size - 1; i++)
    {
        values.push_back(1);
        values.push_back(-k_offdiag_lower);
    }

    // lower K boundary, upper V boundary
    // again, V boundary condition is prioritized
    values.push_back(-v_offdiag_upper);
    values.push_back(1);


    // inner K
    for(int j = 1; j < k_size - 1; j++)
    {
        // lower V boundary
        values.push_back(1);
        values.push_back(-v_offdiag_lower);

        // inner V
        for(int i = 1; i < v_size - 1; i++)
        {
            double dv_left  = v_grid[i][j] - v_grid[i-1][j];
            double dv_right = v_grid[i+1][j] - v_grid[i][j];
            double dk_left  = k_grid[i][j] - k_grid[i][j-1];
            double dk_right = k_grid[i][j+1] - k_grid[i][j];
            // f(i-1,j-1) coeff
            double coeff = 0.25 * (jacobian[i-1][j] * Dvk[i-1][j] + jacobian[i][j-1] * Dkv[i][j-1]) / (jacobian[i][j] * dv_left * dk_left);
            values.push_back(-coeff);

            // f(i,j-1) coeff
            double factor = 0.25 / (jacobian[i][j] * dk_left);
            coeff = factor * (
                2 * (jacobian[i][j]   * Dkk[i][j]   / dk_right + jacobian[i][j-1] * Dkk[i][j-1] / dk_left)
                  + (jacobian[i][j]   * Dvk[i][j]   / dv_right - jacobian[i][j]   * Dvk[i][j]   / dv_left)
                  + (jacobian[i][j-1] * Dkv[i][j-1] / dv_right - jacobian[i][j-1] * Dkv[i][j-1] / dv_left)
            );
            values.push_back(-coeff);

            // f(i+1,j-1) coeff
            coeff = -0.25 * (jacobian[i+1][j] * Dvk[i+1][j] + jacobian[i][j-1] * Dkv[i][j-1]) / (jacobian[i][j] * dv_right * dk_left);
            values.push_back(-coeff);

            

            // f(i-1,j) coeff
            factor = 0.25 / (jacobian[i][j] * dv_left);
            coeff = factor * (
                2 * (jacobian[i][j]   * Dvv[i][j]   / dv_right + jacobian[i-1][j] * Dvv[i-1][j] / dv_left)
                  + (jacobian[i-1][j] * Dvk[i-1][j] / dk_right - jacobian[i-1][j] * Dvk[i-1][j] / dk_left)
                  + (jacobian[i][j]   * Dkv[i][j]   / dk_right - jacobian[i][j]   * Dkv[i][j]   / dk_left)
            );
            values.push_back(-coeff);

            // f(i,j) coeff, 1/dt and loss
            factor = -0.5 / jacobian[i][j];
            coeff = factor * (
                jacobian[i+1][j] * Dvv[i+1][j] / (dv_right * dv_right) + 2 * jacobian[i][j] * Dvv[i][j] / (dv_left * dv_right) + jacobian[i-1][j] * Dvv[i-1][j] / (dv_left * dv_left)
              + jacobian[i][j+1] * Dkk[i][j+1] / (dk_right * dk_right) + 2 * jacobian[i][j] * Dkk[i][j] / (dk_left * dk_right) + jacobian[i][j-1] * Dkk[i][j-1] / (dk_left * dk_left)
              - (1 / (dv_right * dk_right) - 1 / (dv_right * dk_left) - 1 / (dv_left * dk_right) + 1 / (dv_left * dk_left)) * (Dvk[i][j] + Dkv[i][j]) * jacobian[i][j] / 2
            );
            values.push_back(1/dt - loss[i][j] - coeff);
                        
            // f(i+1,j) coeff
            factor = 0.25 / (jacobian[i][j] * dv_right);
            coeff = factor * (
                2 * (jacobian[i+1][j] * Dvv[i+1][j] / dv_right + jacobian[i][j]   * Dvv[i][j]   / dv_left)
                  + (jacobian[i+1][j] * Dvk[i+1][j] / dk_left  - jacobian[i+1][j] * Dvk[i+1][j] / dk_right)
                  + (jacobian[i][j]   * Dkv[i][j]   / dk_left  - jacobian[i][j]   * Dkv[i][j]   / dk_right)
            );
            values.push_back(-coeff);

            // f(i-1,j+1) coeff
            coeff = -0.25 * (jacobian[i-1][j] * Dvk[i-1][j] + jacobian[i][j+1] * Dkv[i][j+1]) / (jacobian[i][j] * dv_left * dk_right);
            values.push_back(-coeff);
            
            // f(i,j+1) coeff
            factor = 0.25 / (jacobian[i][j] * dk_right);
            coeff = factor * (
                2 * (jacobian[i][j+1] * Dkk[i][j+1] / dk_right + jacobian[i][j]   * Dkk[i][j]   / dk_left)
                  +  jacobian[i][j+1] * Dvk[i][j+1] / dv_left  - jacobian[i][j]   * Dvk[i][j]   / dv_right
                  +  jacobian[i][j]   * Dkv[i][j]   / dv_left  - jacobian[i][j+1] * Dkv[i][j+1] / dv_right
            );
            values.push_back(-coeff);

            // f(i+1,j+1) coeff
            coeff = 0.25 * (jacobian[i+1][j] * Dvk[i+1][j] + jacobian[i][j+1] * Dkv[i][j+1]) / (jacobian[i][j] * dv_right * dk_right);
            values.push_back(-coeff);
        }

        // upper V boundary
        values.push_back(-v_offdiag_upper);
        values.push_back(1);
    }
    // upper K boundary, lower V boundary
    values.push_back(1);
    values.push_back(-v_offdiag_lower);

    // upper K boundary, inner V
    for(int i = 1; i < v_size - 1; i++)
    {
        values.push_back(-k_offdiag_upper);
        values.push_back(1);
    }

    // upper K boundary, upper V boundary
    values.push_back(-v_offdiag_upper);
    values.push_back(1);
}

void initialize_sparse_indices(int v_size, int k_size, std::vector<int>& column_indices, std::vector<int>& rows_csr)
{
    column_indices.clear();
    column_indices.reserve(
        2 * (2 * (v_size - 2)) // two blocks in the beginning and two in the end, coresponding to lower and upper K boundary
        + (k_size - 2) * (v_size + 2 * (v_size - 1)) // inner rows, main diag block with one main diagonal and two off-diagonals
        + 2 * (k_size - 2) * (3 * (v_size - 2))
    );
    rows_csr.clear();
    rows_csr.reserve(v_size * k_size + 1);

    // lower K boundary, lower V boundary
    rows_csr.push_back(0);
    column_indices.push_back(0);
    column_indices.push_back(1);
    rows_csr.push_back(rows_csr.back() + 2);

    // lower K boundary, inner V
    for(int i = 1; i < v_size - 1; i++)
    {
        column_indices.push_back(i);
        column_indices.push_back(v_size + i);
        rows_csr.push_back(rows_csr.back() + 2);
    }

    // lower K boundary, upper V boundary
    column_indices.push_back(v_size - 2);
    column_indices.push_back(v_size - 1);
    rows_csr.push_back(rows_csr.back() + 2);

    // inner K
    for(int j = 1; j < k_size - 1; j++)
    {
        // lower V boundary
        column_indices.push_back(j * v_size);
        column_indices.push_back(j * v_size + 1);
        rows_csr.push_back(rows_csr.back() + 2);

        // inner V
        for(int i = 1; i < v_size - 1; i++)
        {
            column_indices.push_back((j-1) * v_size + i-1);
            column_indices.push_back((j-1) * v_size + i);
            column_indices.push_back((j-1) * v_size + i+1);

            column_indices.push_back(j * v_size + i-1);
            column_indices.push_back(j * v_size + i);
            column_indices.push_back(j * v_size + i+1);

            column_indices.push_back((j+1) * v_size + i-1);
            column_indices.push_back((j+1) * v_size + i);
            column_indices.push_back((j+1) * v_size + i+1);

            rows_csr.push_back(rows_csr.back() + 9);
        }
        // upper V boundary
        column_indices.push_back((j + 1) * v_size - 2);
        column_indices.push_back((j + 1) * v_size - 1);
        rows_csr.push_back(rows_csr.back() + 2);
    }
    // upper K boundary, lower V boundary
    column_indices.push_back((k_size - 1) * v_size);
    column_indices.push_back((k_size - 1) * v_size + 1);
    rows_csr.push_back(rows_csr.back() + 2);

    // upper K boundary, inner V
    for(int i = 1; i < v_size - 1; i++)
    {
        column_indices.push_back((k_size - 2) * v_size + i);
        column_indices.push_back((k_size - 1) * v_size + i);
        rows_csr.push_back(rows_csr.back() + 2);
    }

    // upper K boundary, upper V boundary
    column_indices.push_back(k_size * v_size - 2);
    column_indices.push_back(k_size * v_size - 1);
    rows_csr.push_back(rows_csr.back() + 2);
}

void initialize_rhs(std::vector<double>& rhs, const Matrix2D<double>& psd, const Matrix1D<double>& v_lower, const Matrix1D<double>& v_upper, const Matrix1D<double>& k_lower, const Matrix1D<double>& k_upper, const Matrix2D<double>& source, double dt)
{
    rhs.clear();
    rhs.reserve(psd.size_q1 * psd.size_q2);
    for(size_t i = 0; i < psd.size_q1; i++)
    {
        rhs.push_back(k_lower[i]);
    }
    for(size_t j = 1; j < psd.size_q2 - 1; j++)
    {
        rhs.push_back(v_lower[j]);
        for(size_t i = 1; i < psd.size_q1 - 1; i++)
        {
            rhs.push_back(psd[i][j] / dt + source[i][j]);
        }
        rhs.push_back(v_upper[j]);
    }
    for(size_t i = 0; i < psd.size_q1; i++)
    {
        rhs.push_back(k_upper[i]);
    }
}

void Diffusion_2D_MKL(Matrix2D<double>& psd, const Matrix2D<double>& v_grid, const Matrix2D<double>& k_grid,
        BoundaryConditionType v_lower_type, BoundaryConditionType v_upper_type,
        BoundaryConditionType k_lower_type, BoundaryConditionType k_upper_type,
        const Matrix1D<double>& v_lower, const Matrix1D<double>& v_upper, const Matrix1D<double>& k_lower, const Matrix1D<double>& k_upper,
        const Matrix2D<double>& Dvv, const Matrix2D<double>& Dvk, const Matrix2D<double>& Dkv, const Matrix2D<double>& Dkk,
        const Matrix2D<double>& jacobian, const Matrix2D<double>& loss, const Matrix2D<double>& source, double dt,
        sparse_matrix_t* csrA)
{
    int mat_size = psd.size_q1 * psd.size_q2;
    std::vector<double> sparse_values;
    initialize_sparse_values(
        v_grid, k_grid, v_lower_type, v_upper_type, k_lower_type, k_upper_type,
        Dvv, Dvk, Dkv, Dkk, jacobian, loss, dt, sparse_values
    );

    std::vector<double> rhs;
    initialize_rhs(rhs, psd, v_lower, v_upper, k_lower, k_upper, source, dt);

    sparse_status_t status = mkl_sparse_d_qr_factorize(*csrA, sparse_values.data());
    if(status != SPARSE_STATUS_SUCCESS)
	{
		std::cout << "MKL factorize error " << status << '\n';
		exit(EXIT_FAILURE);
	}
    
    Matrix1D<double> dummy(mat_size);
    status = mkl_sparse_d_qr_solve(
        SPARSE_OPERATION_NON_TRANSPOSE, *csrA, sparse_values.data(),
        SPARSE_LAYOUT_COLUMN_MAJOR, 1, &dummy[0], mat_size, rhs.data(), 1
    );
    if(status != SPARSE_STATUS_SUCCESS)
	{
		std::cout << "MKL solve error " << status << '\n';
		exit(EXIT_FAILURE);
	}
    for(size_t i = 0; i < psd.size_q1; i++)
    {
        for(size_t j = 0; j < psd.size_q2; j++)
        {
            psd[i][j] = dummy[j * psd.size_q1 + i];
        }
    }
}
#endif