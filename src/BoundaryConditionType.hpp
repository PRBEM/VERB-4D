// SPDX-FileCopyrightText: 2015 UCLA
// SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * \file BoundaryConditionType.hpp
 *
 * \brief Boundary condition types for numerical solvers
 *
 * This file defines the boundary condition types used throughout the VERB4D solver
 * for handling domain boundaries in diffusion and convection equations. The boundary
 * conditions are essential for well-posed mathematical problems and physical realism.
 *
 * Created on: May 28, 2011
 * Author: dimath
 */

#pragma once

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

/**
 * @brief Enumeration of boundary condition types for numerical solvers
 * 
 * @details This enum defines the types of boundary conditions that can be applied
 * at domain boundaries in the VERB4D solver. Each type corresponds to a different
 * mathematical constraint on the solution at the boundary.
 * 
 * **Mathematical Formulations:**
 * 
 * For a function f(x) at boundary point x_b:
 * 
 * - **Periodic:** \f$ f(x_{left}) = f(x_{right}) \f$
 * - **ConstantValue (Dirichlet):** \f$ f(x_b) = f_{BC} \f$
 * - **ConstantDerivative (Neumann):** \f$ \frac{\partial f}{\partial x}\big|_{x_b} = g_{BC} \f$
 * 
 * **Physical Interpretations:**
 * 
 * - **Periodic:** Used for cyclic domains
 * - **ConstantValue:** Fixed values at boundaries 
 * - **ConstantDerivative:** Fixed gradient at boundaries
 * 
 * **Usage in Solvers:**
 * 
 * These boundary conditions are implemented in the matrix assembly routines
 * of diffusion and convection solvers, where they modify the coefficient
 * matrices to enforce the appropriate mathematical constraints.
 */
enum class BoundaryConditionType { 
    /**
     * @brief Periodic boundary conditions
     * 
     * @details The solution values at opposite boundaries are equal, creating
     * a cyclic domain. Mathematically: f(x_left) = f(x_right).
     * 
     * **Implementation:** Special handling in matrix assembly to connect
     * boundary points across the domain.
     */
    Periodic, 
    
    /**
     * @brief Constant value (Dirichlet) boundary conditions
     * 
     * @details The solution value is fixed at the boundary point.
     * Mathematically: f(x_boundary) = f_BC where `f_BC` is the prescribed value.
     * 
     * **Implementation:** Matrix row is replaced with identity equation
     * f[boundary] = boundary_value.
     */
    ConstantValue, 
    
    /**
     * @brief Constant derivative (Neumann) boundary conditions
     * 
     * @details The normal derivative is fixed at the boundary point.
     * Mathematically: \f$ \frac{\partial f}{\partial n}|_{boundary} = g_{BC} \f$ where `g_BC` is the prescribed gradient.
     * 
     * **Implementation:** Finite difference approximation of derivative
     * using neighboring points.
     * 
     * @warning Current implementation only works correctly for zero derivative (g_BC = 0)
     */
    ConstantDerivative 
};
/** @brief String representation for periodic boundary condition type */
constexpr char PERIODIC_STR[]         = "BCT_PERIODIC";

/** @brief String representation for constant value boundary condition type */
constexpr char CONST_VALUE_STR[]      = "BCT_CONSTANT_VALUE";

/** @brief String representation for constant derivative boundary condition type */
constexpr char CONST_DERIVATIVE_STR[] = "BCT_CONSTANT_DERIVATIVE";

/**
 * @brief Output stream operator for BoundaryConditionType
 * 
 * @details Converts a BoundaryConditionType enum value to its string representation
 * for output to streams (e.g., std::cout, file output, logging).
 * 
 * @param os [in,out] Output stream to write to
 * @param type [in] BoundaryConditionType enum value to convert
 * 
 * @return Reference to the output stream for chaining
 * 
 * @throws std::invalid_argument if an invalid BoundaryConditionType is encountered
 * 
 * @code
 * BoundaryConditionType bc = BoundaryConditionType::ConstantValue;
 * std::cout << "Boundary condition: " << bc << std::endl;
 * // Output: "Boundary condition: BCT_CONSTANT_VALUE"
 * @endcode
 */
inline std::ostream& operator<<(std::ostream& os, const BoundaryConditionType& type)
{
    switch(type) {
        case BoundaryConditionType::Periodic:
            return os << PERIODIC_STR;
        case BoundaryConditionType::ConstantValue:
            return os << CONST_VALUE_STR;
        case BoundaryConditionType::ConstantDerivative:
            return os << CONST_DERIVATIVE_STR;
        default:
            throw(std::invalid_argument("Encountered invalid BoundaryConditionType!"));
    }
}

/**
 * @brief Input stream operator for BoundaryConditionType
 * 
 * @details Parses a string representation of a boundary condition type from an input
 * stream and converts it to the corresponding BoundaryConditionType enum value.
 * This is useful for reading boundary condition types from configuration files
 * or user input.
 * 
 * @param ss [in,out] Input stream to read from
 * @param type [out] BoundaryConditionType enum value to store the result
 * 
 * @return Reference to the input stream for chaining
 * 
 * @throws std::invalid_argument if an invalid boundary condition string is encountered
 * 
 * **Recognized strings:**
 * - "BCT_PERIODIC" → BoundaryConditionType::Periodic
 * - "BCT_CONSTANT_VALUE" → BoundaryConditionType::ConstantValue  
 * - "BCT_CONSTANT_DERIVATIVE" → BoundaryConditionType::ConstantDerivative
 * 
 * @code
 * std::stringstream ss("BCT_CONSTANT_VALUE");
 * BoundaryConditionType bc;
 * ss >> bc;
 * // bc now contains BoundaryConditionType::ConstantValue
 * @endcode
 */
inline std::basic_istream<char>& operator>>(std::basic_istream<char>& ss, BoundaryConditionType& type)
{
    std::string s;
    ss >> s;
    if (s == CONST_VALUE_STR) {
		type = BoundaryConditionType::ConstantValue;
	}
	else if (s == PERIODIC_STR) {
		type = BoundaryConditionType::Periodic;
	}
	else if (s == CONST_DERIVATIVE_STR) {
		type = BoundaryConditionType::ConstantDerivative;
	}
    else throw(std::invalid_argument("Encountered invalid BoundaryConditionType!"));
    return ss;
}
