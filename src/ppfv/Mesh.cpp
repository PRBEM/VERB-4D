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

#include "Mesh.h"

Mesh::Mesh(const Eigen::VectorXd& x, const Eigen::VectorXd& y,
  double dt)
: x_centers_(x), y_centers_(y), dt_(dt) {

  nx_ = x_centers_.size();
  ny_ = y_centers_.size();

  dt_ = dt; 
  
  dx_.resize(nx_);
  dy_.resize(ny_);

  calculate_mesh_spacing();

  xO_ = x_centers_(0) - dx(0)/2.0;
  yO_ = y_centers_(0) - dy(0)/2.0;

  x_vertices_ = Eigen::VectorXd::Zero(nx()+1);
  y_vertices_ = Eigen::VectorXd::Zero(ny()+1);

  x_vertices_(0) = xO(); 
  y_vertices_(0) = yO(); 

  for (std::size_t i=1; i<nx(); ++i) x_vertices_(i) = (x_centers_(i-1) + x_centers_(i)) / 2.0; 
  for (std::size_t j=1; j<ny(); ++j) y_vertices_(j) = (y_centers_(j-1) + y_centers_(j)) / 2.0; 

  x_vertices_(nx()) = x_vertices_(nx()-1) + dx(nx()-1);
  y_vertices_(ny()) = y_vertices_(ny()-1) + dy(ny()-1);

  nbr_inds.resize({nx(), ny(), nnbrs()});
  edges.resize({nx(), ny(), nnbrs()});

  build_connectivity(); 

  // The reverse inbr number.
  // For example, if the current cell is K, its 0th neighbor is L.
  // Then, for cell L, K is its 2th neighbor.
  rinbr_(0) = 2; 
  rinbr_(1) = 3; 
  rinbr_(2) = 0; 
  rinbr_(3) = 1; 
}

void Mesh::calculate_mesh_spacing()
{
  // Calculate individual dx for each cell

  // Calculate individual dy for each cell
  for (std::size_t i = 0; i < nx_; ++i) {
    if (i == 0) {
      dx_(i) = x_centers_[1] - x_centers_[0];
    } else if (i == nx_ - 1) {
      dx_(i) = x_centers_[i] - x_centers_[i-1];
    } else {
      dx_(i) = (x_centers_[i + 1] - x_centers_[i]);
    }
  }

  // Calculate individual dy for each cell
  for (std::size_t j = 0; j < ny_; ++j) {
    if (j == 0) {
      dy_(j) = y_centers_[1] - y_centers_[0];
    } else if (j == ny_ - 1) {
      dy_(j) = y_centers_[j] - y_centers_[j-1];
    } else {
      dy_(j) = (y_centers_[j + 1] - y_centers_[j]);
    }
  }
}

void Mesh::build_connectivity() {

  int inbr; 

  Point A, B; 
  Vector2 dr; 

  bool is_flipped_x = x_vertices_(0) > x_vertices_(1);
  bool is_flipped_y = y_vertices_(0) > y_vertices_(1);

  for (std::size_t i=0; i<nx(); ++i) {
    for (std::size_t j=0; j<ny(); ++j) {
      // nbr 0
      inbr = 0; 
      nbr_inds(i,j,inbr).i = i-1;
      nbr_inds(i,j,inbr).j = j; 

      A = {x_vertices_(i), y_vertices_(j+1)}; 
      B = {x_vertices_(i), y_vertices_(j)}; 

      if (is_flipped_x) {
        edges(i,j,inbr).set_vs_dir({B,A}, XPOS); 
      } else {
        edges(i,j,inbr).set_vs_dir({A,B}, XNEG); 
      }

      // nbr 1
      inbr = 1;
      nbr_inds(i,j,inbr).i = i;
      nbr_inds(i,j,inbr).j = j+1;

      B = A;
      A = {x_vertices_(i+1), y_vertices_(j+1)};

      if (is_flipped_y) {
        edges(i,j,inbr).set_vs_dir({B,A}, YNEG); 
      } else {
        edges(i,j,inbr).set_vs_dir({A,B}, YPOS); 
      }

      // nbr 2
      inbr = 2; 
      nbr_inds(i,j,inbr).i = i+1;
      nbr_inds(i,j,inbr).j = j;

      B = A;
      A = {x_vertices_(i+1), y_vertices_(j)}; 

      if (is_flipped_x) {
        edges(i,j,inbr).set_vs_dir({B,A}, XNEG); 
      } else {
        edges(i,j,inbr).set_vs_dir({A,B}, XPOS);
      }

      // nbr 3
      inbr = 3;
      nbr_inds(i,j,inbr).i = i;
      nbr_inds(i,j,inbr).j = j-1;

      B = A;
      A = {x_vertices_(i), y_vertices_(j)}; 

      if (is_flipped_y) {
        edges(i,j,inbr).set_vs_dir({B,A}, YPOS); 
      } else {
        edges(i,j,inbr).set_vs_dir({A,B}, YNEG);
      }
    }
  }
}