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

  // Vertices at the mid-points between the centers; the outer cells reach as far outward as inward.
  x_vertices_ = Eigen::VectorXd::Zero(nx()+1);
  y_vertices_ = Eigen::VectorXd::Zero(ny()+1);

  for (std::size_t i=1; i<nx(); ++i) x_vertices_(i) = (x_centers_(i-1) + x_centers_(i)) / 2.0;
  for (std::size_t j=1; j<ny(); ++j) y_vertices_(j) = (y_centers_(j-1) + y_centers_(j)) / 2.0;

  x_vertices_(0) = x_centers_(0) - (x_centers_(1) - x_centers_(0)) / 2.0;
  y_vertices_(0) = y_centers_(0) - (y_centers_(1) - y_centers_(0)) / 2.0;
  x_vertices_(nx()) = x_centers_(nx()-1) + (x_centers_(nx()-1) - x_centers_(nx()-2)) / 2.0;
  y_vertices_(ny()) = y_centers_(ny()-1) + (y_centers_(ny()-1) - y_centers_(ny()-2)) / 2.0;

  xO_ = x_vertices_(0);
  yO_ = y_vertices_(0);

  calculate_mesh_spacing();

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
  // Cell widths between the cell's own vertices, so that the cell volume matches its faces
  // (a forward difference of the centers does not on non-uniform grids). Negative for descending grids.
  for (std::size_t i = 0; i < nx_; ++i) dx_(i) = x_vertices_(i + 1) - x_vertices_(i);
  for (std::size_t j = 0; j < ny_; ++j) dy_(j) = y_vertices_(j + 1) - y_vertices_(j);
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