// SPDX-FileCopyrightText: 2026 GFZ Helmholtz Centre for Geosciences
//
// SPDX-License-Identifier: BSD-3-Clause

#include "Diffusion_2D_PPFV.h"

#include <iostream>

Eigen::VectorXd convertToEigen(const Matrix1D<double> &mat)
{
    Eigen::VectorXd eigen_vec(mat.size_q1);
    for (size_t i = 0; i < mat.size_q1; ++i)
    {
        eigen_vec(i) = mat.plane_array[i];
    }
    return eigen_vec;
}

xt::xtensor<double, 2> convertToXt(const Matrix2D<double> &mat)
{
    xt::xtensor<double, 2> xt_mat = xt::zeros<double>({mat.size_q1, mat.size_q2});
    for (size_t i = 0; i < mat.size_q1; ++i)
    {
        for (size_t j = 0; j < mat.size_q2; ++j)
        {
            xt_mat(i, j) = mat.plane_array[i * mat.size_q2 + j];
        }
    }
    return xt_mat;
}

Matrix2D<double> convertToMatrix(const xt::xtensor<double, 2> &xt_mat)
{
    size_t size_q1 = xt_mat.shape(0);
    size_t size_q2 = xt_mat.shape(1);

    Matrix2D<double> mat(size_q1, size_q2);

    for (size_t i = 0; i < size_q1; ++i)
    {
        for (size_t j = 0; j < size_q2; ++j)
        {
            mat.plane_array[i * size_q2 + j] = xt_mat(i, j);
        }
    }

    return mat;
}

bool Diffusion_2D_PPFV(
    Matrix2D<double> &psd,
    const Matrix1D<double> &x, const Matrix1D<double> &y,
    int &x_size, int &y_size,
    const Matrix1D<double> &x_LBC, const Matrix1D<double> &x_UBC,
    const Matrix1D<double> &y_LBC, const Matrix1D<double> &y_UBC,
    BoundaryConditionType x_LBC_type, BoundaryConditionType x_UBC_type,
    BoundaryConditionType y_LBC_type, BoundaryConditionType y_UBC_type,
    const Matrix2D<double> &Dxx, const Matrix2D<double> &Dyy, const Matrix2D<double> &Dxy,
    const Matrix2D<double> &G, const Matrix2D<double> &Sources, const Matrix2D<double> &Losses, double total_time, double dt)
{

    Eigen::VectorXd x_fvm = convertToEigen(x);
    Eigen::VectorXd y_fvm = convertToEigen(y);
    Eigen::VectorXd xlbc_fvm = convertToEigen(x_LBC);
    Eigen::VectorXd xubc_fvm = convertToEigen(x_UBC);
    Eigen::VectorXd ylbc_fvm = convertToEigen(y_LBC);
    Eigen::VectorXd yubc_fvm = convertToEigen(y_UBC);

    xt::xtensor<double, 2> PSD0_fvm = convertToXt(psd);
    xt::xtensor<double, 2> Dxx_fvm = convertToXt(Dxx);
    xt::xtensor<double, 2> Dyy_fvm = convertToXt(Dyy);
    xt::xtensor<double, 2> Dxy_fvm = convertToXt(Dxy);
    xt::xtensor<double, 2> G_fvm = convertToXt(G);
    xt::xtensor<double, 2> Sources_fvm = convertToXt(Sources);
    xt::xtensor<double, 2> Losses_fvm = convertToXt(Losses);


    const int num_substeps = std::max(1, int(total_time / dt));
    dt = total_time / num_substeps;
    bool success;

    Mesh mesh(x_fvm, y_fvm, dt);
    VERB_ppfv2d eq(mesh, PSD0_fvm, Dxx_fvm, Dxy_fvm, Dyy_fvm, xlbc_fvm, xubc_fvm, ylbc_fvm, yubc_fvm, G_fvm, x_LBC_type, x_UBC_type, y_LBC_type, y_UBC_type, Losses_fvm);
    Solver solver(mesh, &eq);
    
    for (int i = 0; i < num_substeps; ++i)
    {
        solver.update();    
    }
    PSD0_fvm = solver.PSD();
    psd = convertToMatrix(PSD0_fvm);
    return success;
}