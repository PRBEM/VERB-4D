/**
 * \file Interpolation.cpp
 * \brief Wrapper around interpolation library, interpolates values into new grid using MonotCubicInterpolator.cpp
 *
 *  Created on: Sep 14, 2013
 *      Author: dimath
 */

#include "Interpolation.h"
#include "MonotCubicInterpolator.hpp"

using namespace std;

/** Wrapper around interpolation library,
 *  Interpolates values into new grid,
 *  Grid size should be the same
 *
 * Using interpolator from Opm::MonotCubicInterpolator
 */
//template <typename T>
Matrix1D<double> Cubic1D(Matrix1D<double> &old_grid, Matrix1D<double> &values, Matrix1D<double> &new_grid) {

	// grids must be same size
	if (new_grid.size_q1 != old_grid.size_q1) {
	    throw("Only interpolation to the same grid size is supported.");
	}

	vector<double> x(old_grid.size_q1);
	vector<double> f(old_grid.size_q1);

	// copying old grid int x and values into f
	for (int i = 0; i < old_grid.size_q1; i++) {
		x.push_back(old_grid[i]);
		f.push_back(values[i]);
	}

	// create an interpolator from x and f
	Opm::MonotCubicInterpolator interpolator(x, f);

	// put the values generated from the interpolator into the newly created grid
	for (int i = 0; i < new_grid.size_q1; i++) {
		values[i] = interpolator.evaluate(new_grid[i]);
	}

	return Matrix1D<double>();
}
