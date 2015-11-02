/**
 * \file Interpolation.h
 * \brief Wrapper around interpolation library, interpolates values into new grid using MonotCubicInterpolator.cpp
 *
 *  Created on: Sep 14, 2013
 *      Author: dimath
 */
#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include "Matrix.h"

using namespace std;


/** One-dimensional interpolation
*/
//template <typename T>
Matrix1D<double> Cubic1D(Matrix1D<double> &old_grid, Matrix1D<double> &f, Matrix1D<double> &new_grid);


#endif /* INTERPOLATION_H_ */
