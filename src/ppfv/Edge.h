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

#ifndef EDGE_H_
#define EDGE_H_

#include "common.h"

enum Direction {XPOS, XNEG, YPOS, YNEG}; 

typedef std::array<Point, 2> Edge_Vertices;

class Edge{
  public:
    static const int NT = 2; // number of points for each edge. 
    const Point& v(int i) const { return vs_[i]; } // vertices
    int v_size() const { return vs_.size(); }                                                   
    double length() const { return length_; }
    const Vector2& n() const { return n_; } // normal vector, depending on the edge direction

    void set_vs_dir(const Edge_Vertices& vs, Direction dir) {
      vs_ = vs;
      length_ = (vs[1] - vs[0]).norm(); 

      switch (dir) {
        case XPOS: // Edge perpendicular to positive x-axis
          n_ = {1, 0};
          break;

        case XNEG:
          n_ = {-1, 0};
          break; 

        case YPOS: // Edge perpendicular to y-axis
          n_ = {0, 1};
          break;

        case YNEG:
          n_ = {0, -1};
          break;

        default: // incorrect dir
          n_ = {0, 0}; 
      }

    }

  private:

    Edge_Vertices vs_; 
    double length_; 
    Vector2 n_; // the unit normal vector
};

#endif /* EDGE_H */

