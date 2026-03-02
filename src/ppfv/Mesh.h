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
 * 
 */

#ifndef MESH_H_
#define MESH_H_

#include <vector>
#include "common.h"
// #include "Parameters.h" 
#include "Edge.h"

struct Ind{
  std::size_t i;
  std::size_t j;
}; 

class Mesh {
  public:
  Mesh(const Eigen::VectorXd& x, const Eigen::VectorXd& y,
    double dt);

    const Eigen::VectorXd& x() const { return x_centers_; }
    const Eigen::VectorXd& y() const { return y_centers_; }

    double x(int i) const { return x_centers_(i); }
    double y(int j) const { return y_centers_(j); }

    double xO() const { return xO_; }
    double yO() const { return yO_; }

    std::size_t nx() const { return nx_; }
    std::size_t ny() const { return ny_; }

    double dx(int i) const { return dx_(i); }
    double dy(int j) const { return dy_(j); }
    double dt() const { return dt_; }

    double cell_area(int i, int j) const { return std::abs(dx(i) * dy(j)); }
    double cell_area_dt(int i, int j) const { return cell_area(i, j) / dt();}
    double cell_area_dt(Ind ind) const { return cell_area(ind.i, ind.j) / dt();}

    void indO(const Point& A, Ind* indp) const { 
      // calculate the i,j coordinate relative to the Origin
      // Note: not the cell index.
      // This function is useful to calculate fA and fB from interpolation

      int i = 0;
      while (A(0) != x_vertices_(i)) {
        i++;
      }

      int j = 0;
      while (A(1) != y_vertices_(j)) {
        j++;
      }

      indp->i = i;
      indp->j = j; 
      
    }

    std::size_t flatten_index(const Ind& ind) const { // map 2d indices to 1
      return ind.j*nx() + ind.i;
    }

    std::size_t nnbrs() const { return 4; } // each cell has 4 nbrs

    // define the neighbor # of 4 adjacent cells
    // im -- (i-1, j, k); jp -- (i, j+1, k)
    // ip -- (i+1, j, k); jm -- (i, j-1, k)
    int inbr_im() const { return 0; }
    int inbr_jp() const { return 1; }
    int inbr_ip() const { return 2; }
    int inbr_jm() const { return 3; }

    int rinbr(int inbr) const { return rinbr_(inbr); }                                    

    void get_nbr_ind(const Ind& ind, int inbr, Ind* nbr_indp) const {
      *nbr_indp = nbr_inds(ind.i,ind.j,inbr); 
    }

    void get_nbr_edge(const Ind& ind, int inbr, Edge* edgep) const {
      *edgep = edges(ind.i,ind.j,inbr); 
    }

  private:
    std::size_t nx_; 
    std::size_t ny_; 
    Eigen::VectorXd dx_;  // Individual dx for each x-cell
    Eigen::VectorXd dy_;  // Individual dy for each y-cell
    double dt_; 

    // coordinate origin: corresponds to i-0.5, j-0.5
    double xO_; 
    double yO_;

    Eigen::VectorXd x_centers_; 
    Eigen::VectorXd y_centers_;

    Eigen::VectorXd x_vertices_; 
    Eigen::VectorXd y_vertices_;

    Eigen::Matrix<int, 4, 1> rinbr_; 

    xt::xtensor<Ind,3> nbr_inds;
    xt::xtensor<Edge,3> edges;

    void build_connectivity();
    void calculate_mesh_spacing();
};

#endif /* MESH_H */

