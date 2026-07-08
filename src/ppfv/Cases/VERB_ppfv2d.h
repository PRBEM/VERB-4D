/*
 * SPDX-FileCopyrightText: 2025 Xin Tao
 * SPDX-FileCopyrightText: 2026 GFZ Helmholtz Centre Potsdam
 * SPDX-FileContributor: Bernhard Haas
 * SPDX-FileContributor: Peng Peng
 * SPDX-FileContributor: Xin Tao
 *
 * SPDX-License-Identifier: BSD-3
 * SPDX-License-Identifier: MIT
 * 
 * This code is adapted from the original version published at https://github.com/xtaohub/Sayram-2D/tree/v1.0.0.
 */

#pragma once

#include "../Equation.h"
#include "../../BoundaryConditionType.hpp"

class VERB_ppfv2d : public Equation
{
public:
  VERB_ppfv2d(const Mesh &m_in,
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
             const Xtensor2d &Losses);

  double init_PSD(const Ind &ind) const;
  void apply_bcs(Xtensor2d *vertex_fp) const;
  BoundaryConditionType x_LBC_type;
  BoundaryConditionType x_UBC_type;
  BoundaryConditionType y_LBC_type;
  BoundaryConditionType y_UBC_type;
  BoundaryConditionType get_x_LBC_type() const;
  BoundaryConditionType get_x_UBC_type() const;
  BoundaryConditionType get_y_LBC_type() const;
  BoundaryConditionType get_y_UBC_type() const;

  private:
  const Mesh &m; 
  const Xtensor2d &PSD;
  const Eigen::VectorXd &x_LBC_, &x_UBC_, &y_LBC_, &y_UBC_;
  double calculate_init_PSD(const Ind &ind) const
  {
    // convert from cell indices (ind) to vertex indices

    int idx_left = ind.i;
    int idx_right = ind.i + 1;

    int idx_bottom = ind.j;
    int idx_top = ind.j+1;
  
    return PSD(ind.i, ind.j);
  }

  void init();
};
