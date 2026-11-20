// SPDX-FileCopyrightText: 2025 Xin Tao
// SPDX-FileCopyrightText: 2026 GFZ Helmholtz Centre Potsdam
// 
// SPDX-FileContributor: Bernhard Haas
// SPDX-FileContributor: Peng Peng
// SPDX-FileContributor: Xin Tao
//
// SPDX-License-Identifier: BSD-3
// SPDX-License-Identifier: MIT
//
// This code is adapted from the original version published at https://github.com/xtaohub/Sayram-2D/tree/v1.0.0.

#include "Solver.h"
#include <iostream>
#include <iomanip> // for setw
#include <stdexcept>
#include "../BoundaryConditionType.hpp"

Solver::Solver(const Mesh &m_in, Equation *eqp)
    : m(m_in), eq(*eqp)
{

  size_t nx = m.nx();
  size_t ny = m.ny();

  M_.resize(nx * ny, nx * ny);

  PSD_.resize({nx, ny});
  PSD_old_.resize({nx, ny});

  R_.resize(nx * ny);
  ftmp_.resize(nx * ny);

  Lambda_.resize({nx, ny});
  vertex_f_.resize({nx + 1, ny + 1});

  init();
}

void Solver::init()
{

  // iterSolver.setTolerance(1e-8);

  for (size_t i = 0; i < m.nx(); i++)
  {
    for (size_t j = 0; j < m.ny(); j++)
    {
      PSD_(i, j) = eq.init_PSD({i, j});
    }
  }

  // apply boundary conditions

  update_Lambda();
  update_vertex_f();

}

void Solver::print_matrix(const Xtensor2d &mat, const string &name, const int iter) const
{
  if (iter > 0)
  {
    cout << name << iter << endl;
  }
  else
  {
    cout << name << endl;
  }

  const int width = 8;

  for (size_t i = 0; i < m.nx(); i++)
  {
    for (size_t j = 0; j < m.ny(); j++)
    {
      cout << setw(width) << std::fixed << setprecision(5) << mat(i, j) << " ";
    }
    cout << endl;
  }
}

void Solver::print_vertex_matrix(const Xtensor2d &mat, const string &name) const
{
  cout << name << endl;
  const int width = 6;

  for (size_t i = 1; i < m.nx(); i++)
  {
    for (size_t j = 1; j < m.ny(); j++)
    {
      cout << setw(width) << mat(i, j) << " ";
    }
    cout << endl;
  }
}

void Solver::update_Lambda()
{

  for (size_t i = 0; i < m.nx(); i++)
  {
    for (size_t j = 0; j < m.ny(); j++)
    {
      Lambda_(i, j) << eq.Dxx({i, j}) * eq.G({i, j}), eq.Dxy({i, j}) * eq.G({i, j}),
          eq.Dxy({i, j}) * eq.G({i, j}), eq.Dyy({i, j}) * eq.G({i, j});
    }
  }
}

void Solver::a_sigma_func(const Ind &ind, int inbr, Vector2 *a_sigma_i_p, double *a_sigmap)
{

  Edge edge;
  m.get_nbr_edge(ind, inbr, &edge);

  Point K = {m.x(ind.i), m.y(ind.j)};

  Point A = edge.v(0);
  Point B = edge.v(1);

  Ind indA, indB;

  Vector2 vkb = B - K; // as Vector KB (from K point to B)
  Vector2 vka = A - K;

  Vector2 rvkb = {vkb(1), -vkb(0)}; // (x,y) rotated by 90 clockwise, becoming (y, -x)
  Vector2 rvka = {vka(1), -vka(0)};

  double a_sigma_a = edge.length() * (edge.n().transpose() * Lambda(ind) * rvkb)(0) / (vka.transpose() * rvkb)(0);
  double a_sigma_b = edge.length() * (edge.n().transpose() * Lambda(ind) * rvka)(0) / (vkb.transpose() * rvka)(0);

  if (std::isnan(a_sigma_a) || std::isnan(a_sigma_b)) {
    exit(-1);
  }

  (*a_sigma_i_p)(0) = a_sigma_a;
  (*a_sigma_i_p)(1) = a_sigma_b;

  m.indO(A, &indA);
  m.indO(B, &indB);

  *a_sigmap = a_sigma_a * vertex_f(indA) + a_sigma_b * vertex_f(indB);
}

void Solver::update_coeff_inner_pair(const Ind &ind, int inbr)
{ // add coefficient from a inner (no-a-boundary) face

  Ind nbr_ind;

  int rinb;

  double a_sigma_K, a_sigma_L;
  double mu_K, mu_L;

  double B_sigma, B_sigma_abs, B_sigma_plus, B_sigma_minus;
  double A_K, A_L;

  Vector2 a_sigma_i_K;
  Vector2 a_sigma_i_L;

  // calculate A_K
  a_sigma_func(ind, inbr, &a_sigma_i_K, &a_sigma_K);

  // calculate A_L
  m.get_nbr_ind(ind, inbr, &nbr_ind);
  rinb = m.rinbr(inbr);
  a_sigma_func(nbr_ind, rinb, &a_sigma_i_L, &a_sigma_L);

  calculate_mu(a_sigma_K, a_sigma_L, &mu_K, &mu_L);

  B_sigma = mu_L * a_sigma_L - mu_K * a_sigma_K;

  B_sigma_abs = abs(B_sigma);
  B_sigma_plus = (B_sigma_abs + B_sigma) / 2.0;
  B_sigma_minus = (B_sigma_abs - B_sigma) / 2.0;

  A_K = mu_K * (a_sigma_i_K[0] + a_sigma_i_K[1]) + B_sigma_plus / (PSD(ind) + f_eps_);
  A_L = mu_L * (a_sigma_i_L[0] + a_sigma_i_L[1]) + B_sigma_minus / (PSD(nbr_ind) + f_eps_);

  if (std::isnan(A_K) || std::isnan(A_L)) {
    exit(-1);
  }

  size_t K = m.flatten_index(ind), L = m.flatten_index(nbr_ind);

  M_coeffs_.push_back(T(K, K, A_K));
  M_coeffs_.push_back(T(K, L, -A_L));
  M_coeffs_.push_back(T(L, L, A_L));
  M_coeffs_.push_back(T(L, K, -A_K));
}

void Solver::update_coeff_dirbc(const Ind &ind, int inbr)
{ // add coefficient from a Dirichlet face

  double a_sigma_K;

  double B_sigma, B_sigma_abs, B_sigma_plus, B_sigma_minus;
  double A_K;
  Vector2 a_sigma_i_K;

  // calculate A_K
  a_sigma_func(ind, inbr, &a_sigma_i_K, &a_sigma_K);

  B_sigma = -a_sigma_K;
  B_sigma_abs = abs(B_sigma);
  B_sigma_plus = (B_sigma_abs + B_sigma) / 2.0;
  B_sigma_minus = (B_sigma_abs - B_sigma) / 2.0;

  A_K = a_sigma_i_K[0] + a_sigma_i_K[1] + B_sigma_plus / (PSD(ind) + f_eps_);

  size_t K = m.flatten_index(ind);
  M_coeffs_.push_back(T(K, K, A_K));
  R_(K) += B_sigma_minus;
}

void Solver::assemble()
{ // obtain M and R

  size_t i, j;

  // i direction
  // i==0: dirichlet boundary condition
  if (eq.get_x_LBC_type() == BoundaryConditionType::ConstantValue)
  {
    i = 0;
    for (j = 0; j < m.ny(); ++j)
      update_coeff_dirbc({i, j}, m.inbr_im());
  }

  for (i = 1; i < m.nx(); ++i)
    for (j = 0; j < m.ny(); ++j)
      update_coeff_inner_pair({i, j}, m.inbr_im());

  if (eq.get_x_UBC_type() == BoundaryConditionType::ConstantValue)
  {
    i = m.nx() - 1;
    for (j = 0; j < m.ny(); ++j)
      update_coeff_dirbc({i, j}, m.inbr_ip());
  }

  if (eq.get_y_LBC_type() == BoundaryConditionType::ConstantValue)
  {
    j = 0;
    for (i = 0; i < m.nx(); ++i)
      update_coeff_dirbc({i, j}, m.inbr_jm());
  }

  for (i = 0; i < m.nx(); ++i)
    for (j = 1; j < m.ny(); ++j)
      update_coeff_inner_pair({i, j}, m.inbr_jm());

  if (eq.get_y_UBC_type() == BoundaryConditionType::ConstantValue)
  {
    j = m.ny() - 1;
    for (i = 0; i < m.nx(); ++i)
      update_coeff_dirbc({i, j}, m.inbr_jp());
  }

  // no updates for BoundaryConditionType::ConstantDerivative condition

  
  long K;
  double UKK;

  for (size_t i = 0; i < m.nx(); ++i)
  {
    for (size_t j = 0; j < m.ny(); ++j)
    {
      K = m.flatten_index({i, j});
      UKK = eq.G({i, j}) * m.cell_area_dt({i,j});
      M_coeffs_.push_back(T(K, K, UKK));

      M_coeffs_.push_back(T(K, K, -m.cell_area(i,j) * eq.G({i, j}) * eq.Losses({i, j})));

      // M_(K, K) += UKK;
      UKK = eq.G({i, j}) * m.cell_area_dt(i,j);
      R_(K) += UKK * PSD_old_(i, j);
    }
  }

  M_.setFromTriplets(M_coeffs_.begin(), M_coeffs_.end());
}

double Solver::picard_iteration()
{
  // keeps B / (f + f_eps_) finite for f = 0, relative to the PSD so that the scheme does not depend on its units
  double psd_max = 0;
  for (size_t i = 0; i < m.nx(); ++i)
    for (size_t j = 0; j < m.ny(); ++j)
      psd_max = std::max(psd_max, std::abs(PSD_(i, j)));
  f_eps_ = gEPS * std::max(psd_max, std::numeric_limits<double>::min());

  R_.setZero();
  M_coeffs_.clear();

  assemble();

  // the coefficients of M depend on the iterate, so M is factorized again for every solve
  if (!pattern_analyzed_) {
    solver.analyzePattern(M_);
    pattern_analyzed_ = true;
  }
  solver.factorize(M_);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("PPFV: factorization of the diffusion matrix failed");
  }

  ftmp_ = solver.solve(R_);

  double f_max = 0;
  for (Eigen::Index k = 0; k < ftmp_.size(); ++k) f_max = std::max(f_max, std::abs(ftmp_(k)));
  const double floor = gPICARD_FLOOR_REL * f_max + std::numeric_limits<double>::min();

  double change = 0;
  for (size_t i = 0; i < m.nx(); ++i)
  {
    for (size_t j = 0; j < m.ny(); ++j)
    {
      const double f_new = ftmp_(m.flatten_index({i, j}));
      change = std::max(change, std::abs(f_new - PSD_(i, j)) / (std::abs(f_new) + floor));
      PSD_(i, j) = f_new;
    }
  }

  update_vertex_f();
  return change;
}

void Solver::update()
{
  // Backward Euler step with the nonlinear (Picard) iteration of the coefficients. The coefficients of M
  // depend on f through the vertex values and through B / f; built from the old f alone, a cell that is
  // empty at the start of the step (e.g. inside the loss cone) gets B / (f + eps) -> large and stays empty,
  // since the new f is then proportional to the old one. Iterating lets it fill within the step.
  PSD_old_ = PSD_;

  for (int iter = 0; iter < gPICARD_MAX_ITER; ++iter)
  {
    if (picard_iteration() < gPICARD_RTOL) break;
  }
}

void Solver::update_vertex_f()
{

  for (size_t i = 1; i < m.nx(); ++i)
  {
    for (size_t j = 1; j < m.ny(); ++j)
    {
      vertex_f_(i, j) = (PSD_(i - 1, j - 1) + PSD_(i - 1, j) + PSD_(i, j - 1) + PSD_(i, j)) / 4.0;
    }
  }

  eq.apply_bcs(&vertex_f_, PSD_);
}
