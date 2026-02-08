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

#ifndef EQUATION_H_
#define EQUATION_H_

#include "common.h"
#include "Mesh.h"
#include "../BoundaryConditionType.hpp"

//
// Definitions of the Equation are given here
// including D and boundary conditions.
// To handle different cases, we Equation as an abstract class.
//

class Equation{
  public:
    Equation(const Mesh& m) {
      std::size_t nx = m.nx();
      std::size_t ny = m.ny();

      G_.resize({nx,ny});
      Losses_.resize({nx,ny});

      Dxx_.resize({nx,ny});
      Dyy_.resize({nx,ny});
      Dxy_.resize({nx,ny});
    }

    double Losses(const Ind& ind) const {return Losses_(ind.i, ind.j); }

    double G(const Ind& ind) const { return G_(ind.i, ind.j); } 

    double Dxx(const Ind& ind) const { return Dxx_(ind.i, ind.j); }
    double Dyy(const Ind& ind) const { return Dyy_(ind.i, ind.j); } 
    double Dxy(const Ind& ind) const { return Dxy_(ind.i, ind.j); }


    // pure virtual functions to be defined by the CASE of interest. 
    virtual double init_PSD(const Ind& ind) const = 0;
    virtual void apply_bcs(Xtensor2d* vertex_fp) const = 0;
    virtual BoundaryConditionType get_x_LBC_type() const = 0;
    virtual BoundaryConditionType get_x_UBC_type() const = 0;
    virtual BoundaryConditionType get_y_LBC_type() const = 0;
    virtual BoundaryConditionType get_y_UBC_type() const = 0;

  protected:
    Xtensor2d G_; 
    Xtensor2d Losses_;

    Xtensor2d Dxx_; 
    Xtensor2d Dyy_; 
    Xtensor2d Dxy_;
};


#endif /* EQUATION_H */

