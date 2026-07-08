/*
 * SPDX-FileCopyrightText: 2025 Xin Tao
 * SPDX-FileCopyrightText: 2026 GFZ Helmholtz Centre Potsdam
 * 
 * SPDX-FileContributor: Bernhard Haas
 * SPDX-FileContributor: Peng Peng
 * SPDX-FileContributor: Xin Tao
 *
 * SPDX-License-Identifier: BSD-3
 * SPDX-License-Identifier: MIT
 * 
 * This code is adapted from the original version published at https://github.com/xtaohub/Sayram-2D/tree/v1.0.0.
 */

#ifndef SOLVER_H
#define SOLVER_H

#include "common.h"
#include "Mesh.h"
#include "Equation.h"

class Solver
{
public:
  Solver(const Mesh &m_in, Equation *eqp);

  void update();
  const Xtensor2d &PSD() const { return PSD_; }
  const Xtensor2d &PSD_vertices() const { return vertex_f_; }
  double PSD(const Ind &ind) const { return PSD_(ind.i, ind.j); }
  double PSD(int i, int j) const { return PSD_(i, j); }
  void print_matrix(const Xtensor2d &mat, const std::string &name, const int iter) const;
  void print_vertex_matrix(const Xtensor2d &mat, const std::string &name) const;

private:
  const Mesh &m;

  Equation &eq;
  bool already_factorized_ = false;

  // M f = R
  Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> solver;
  SpMat M_;
  std::vector<T> M_coeffs_;

  Xtensor2d PSD_;
  Eigen::VectorXd R_;
  Eigen::VectorXd ftmp_; // used to store results from the Eigen Solver.

  xt::xtensor<Eigen::Matrix2d, 2> Lambda_; // Lambda_(i,j) would be the Lambda Matrix at cell (i,j)
  const Eigen::Matrix2d &Lambda(int i, int j) const { return Lambda_(i, j); }
  const Eigen::Matrix2d &Lambda(const Ind &ind) const { return Lambda_(ind.i, ind.j); }
  void update_Lambda();

  //
  // f at vertices to build a lookup table
  //
  Xtensor2d vertex_f_;

  double vertex_f(int i, int j) const { return vertex_f_(i, j); }
  double vertex_f(const Ind &ind) const { return vertex_f_(ind.i, ind.j); }
  void update_vertex_f();

  void assemble();

  // to calculate coefficients alpha_sigma_i and a_sigma_i
  void a_sigma_func(const Ind &ind, int inbr, Vector2 *a_sigma_i_p, double *a_sigmap);

  // update coefficients in M for cell Ind, and its neighbor.
  // note that coefficients for both cell Ind and its neighbor are updated.
  void update_coeff_inner_pair(const Ind &ind, int inbr);

  // Here: the inbr neighbor is a Dirichlet boundary.
  void update_coeff_dirbc(const Ind &ind, int inbr);

  void calculate_mu(double a_sigma_K, double a_sigma_L, double *mu_Kp, double *mu_Lp)
  {
    double denom = std::abs(a_sigma_K) + std::abs(a_sigma_L) + 2 * gEPS;
    (*mu_Kp) = (std::abs(a_sigma_L) + gEPS) / denom;
    (*mu_Lp) = 1 - (*mu_Kp);
  }

  void init();
};

#endif /* SOLVER_H */
