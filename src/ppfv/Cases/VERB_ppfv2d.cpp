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

void VERB_ppfv2d::apply_bcs(Xtensor2d *vertex_fp) const
{
  Xtensor2d &vertex_f = *vertex_fp;

  // i == 0 and m.nx() boundary condition case
  for (std::size_t j = 0; j < m.ny() + 1; ++j)
  { 
    if (x_LBC_type == BoundaryConditionType::ConstantValue) vertex_f(0, j) = x_LBC_(j);
    else if (x_LBC_type == BoundaryConditionType::ConstantDerivative) vertex_f(0, j) = PSD(0, j);

    if (x_UBC_type == BoundaryConditionType::ConstantValue) vertex_f(m.nx(), j) = x_UBC_(j);
    else if (x_UBC_type == BoundaryConditionType::ConstantDerivative)  vertex_f(m.nx(), j) = PSD(m.nx()-1, j);

  }
    // j == 0 and j == m.ny() boundary
    for (std::size_t i = 0; i < m.nx() + 1; ++i)
    {
      if (y_LBC_type == BoundaryConditionType::ConstantValue) vertex_f(i, 0) = y_LBC_(i);
      else if (y_LBC_type == BoundaryConditionType::ConstantDerivative) vertex_f(i, 0) = PSD(i, 0);
      if (y_UBC_type == BoundaryConditionType::ConstantValue) vertex_f(i, m.ny()) = y_UBC_(i);
      else if (y_UBC_type == BoundaryConditionType::ConstantDerivative)  vertex_f(i, m.ny()) = PSD(i, m.ny()-1);
  }

}

