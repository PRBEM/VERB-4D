/*
 * Interpolation.cpp
 *
 *  Created on: Sep 14, 2013
 *      Author: dimath
 */

#include "Interpolation.h"
#include "MonotCubicInterpolator.hpp"

using namespace std;

/** Wrapper around interpolation library
 *  Interpolates values into new grid
 *  Grid size should be the same
 *
 */
//template <typename T>
Matrix1D<double> Cubic1D(Matrix1D<double> &old_grid, Matrix1D<double> &values, Matrix1D<double> &new_grid) {

	if (new_grid.size_q1 != old_grid.size_q1) {
	    throw("Only interpolation to the same grid size is supported.");
	}

	vector<double> x(old_grid.size_q1);
	vector<double> f(old_grid.size_q1);

	for (int i = 0; i < old_grid.size_q1; i++) {
		x.push_back(old_grid[i]);
		f.push_back(values[i]);
	}

	Opm::MonotCubicInterpolator interpolator(x, f);

	for (int i = 0; i < new_grid.size_q1; i++) {
		values[i] = interpolator.evaluate(new_grid[i]);
	}

	return Matrix1D<double>();
}
