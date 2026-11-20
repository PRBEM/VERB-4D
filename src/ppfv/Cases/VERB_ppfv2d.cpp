// SPDX-FileCopyrightText: 2025 Xin Tao
// SPDX-FileCopyrightText: 2026 GFZ Helmholtz Centre Potsdam
// SPDX-FileContributor: Bernhard Haas
// SPDX-FileContributor: Peng Peng
// SPDX-FileContributor: Xin Tao
//
// SPDX-License-Identifier: BSD-3
// SPDX-License-Identifier: MIT
//
// This code is adapted from the original version published at https://github.com/xtaohub/Sayram-2D/tree/v1.0.0.

#include "VERB_ppfv2d.h"

VERB_ppfv2d::VERB_ppfv2d(const Mesh &m_in,
                       const Xtensor2d &PSD,
                       const Xtensor2d &Dxx,
                       const Xtensor2d &Dxy,
                       const Xtensor2d &Dyy,
                       const Eigen::VectorXd &x_LBC,
                       const Eigen::VectorXd &x_UBC,
                       const Eigen::VectorXd &y_LBC,
                       const Eigen::VectorXd &y_UBC,
                       const Xtensor2d &G,
                       const BoundaryConditionType &x_LBC_type,
                       const BoundaryConditionType &x_UBC_type,
                       const BoundaryConditionType &y_LBC_type,
                       const BoundaryConditionType &y_UBC_type,
                       const Xtensor2d &Losses)
    : Equation(m_in), m(m_in),
      PSD(PSD),
      x_LBC_(x_LBC), x_UBC_(x_UBC),
      y_LBC_(y_LBC), y_UBC_(y_UBC),
      x_LBC_type(x_LBC_type),
      x_UBC_type(x_UBC_type),
      y_LBC_type(y_LBC_type),
      y_UBC_type(y_UBC_type)
{

  Equation::G_ = G;
  Equation::Dxx_ = Dxx;
  Equation::Dyy_ = Dyy;
  Equation::Dxy_ = Dxy;
  Equation::Losses_ = Losses;
}

double VERB_ppfv2d::init_PSD(const Ind &ind) const
{
  return calculate_init_PSD(ind);
}

BoundaryConditionType VERB_ppfv2d::get_x_LBC_type() const
{
  return x_LBC_type;
}
BoundaryConditionType VERB_ppfv2d::get_x_UBC_type() const
{
  return x_UBC_type;
}
BoundaryConditionType VERB_ppfv2d::get_y_LBC_type() const
{
  return y_LBC_type;
}
BoundaryConditionType VERB_ppfv2d::get_y_UBC_type() const
{
  return y_UBC_type;
}

void VERB_ppfv2d::apply_bcs(Xtensor2d *vertex_fp, const Xtensor2d &cell_f) const
{
  Xtensor2d &vertex_f = *vertex_fp;
  const std::size_t nx = m.nx(), ny = m.ny();

  // A boundary vertex takes the mean over the (one or two) boundary cells it touches: of the current PSD
  // for a zero normal derivative, of the boundary values (one per cell) for a constant value. The vertex
  // index j (i) runs to ny (nx), the cell index only to ny - 1 (nx - 1).
  auto along_x_edge = [&](std::size_t i_cell, std::size_t j) {
    const std::size_t j0 = j == 0 ? 0 : j - 1, j1 = j == ny ? ny - 1 : j;
    return (cell_f(i_cell, j0) + cell_f(i_cell, j1)) / 2.0;
  };
  auto along_y_edge = [&](std::size_t i, std::size_t j_cell) {
    const std::size_t i0 = i == 0 ? 0 : i - 1, i1 = i == nx ? nx - 1 : i;
    return (cell_f(i0, j_cell) + cell_f(i1, j_cell)) / 2.0;
  };
  auto boundary_mean = [](const Eigen::VectorXd &values, std::size_t vertex, std::size_t n_cells) {
    const std::size_t c0 = vertex == 0 ? 0 : vertex - 1, c1 = vertex == n_cells ? n_cells - 1 : vertex;
    return (values(c0) + values(c1)) / 2.0;
  };

  // i == 0 and m.nx() boundary condition case
  for (std::size_t j = 0; j < ny + 1; ++j)
  {
    if (x_LBC_type == BoundaryConditionType::ConstantValue) vertex_f(0, j) = boundary_mean(x_LBC_, j, ny);
    else if (x_LBC_type == BoundaryConditionType::ConstantDerivative) vertex_f(0, j) = along_x_edge(0, j);

    if (x_UBC_type == BoundaryConditionType::ConstantValue) vertex_f(nx, j) = boundary_mean(x_UBC_, j, ny);
    else if (x_UBC_type == BoundaryConditionType::ConstantDerivative) vertex_f(nx, j) = along_x_edge(nx - 1, j);
  }

  // j == 0 and j == m.ny() boundary
  for (std::size_t i = 0; i < nx + 1; ++i)
  {
    if (y_LBC_type == BoundaryConditionType::ConstantValue) vertex_f(i, 0) = boundary_mean(y_LBC_, i, nx);
    else if (y_LBC_type == BoundaryConditionType::ConstantDerivative) vertex_f(i, 0) = along_y_edge(i, 0);

    if (y_UBC_type == BoundaryConditionType::ConstantValue) vertex_f(i, ny) = boundary_mean(y_UBC_, i, nx);
    else if (y_UBC_type == BoundaryConditionType::ConstantDerivative) vertex_f(i, ny) = along_y_edge(i, ny - 1);
  }
}

