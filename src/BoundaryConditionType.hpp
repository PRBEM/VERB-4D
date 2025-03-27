// SPDX-FileCopyrightText: 2015 UCLA
// SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

enum class BoundaryConditionType { Periodic, ConstantValue, ConstantDerivative };
constexpr char PERIODIC_STR[]         = "BCT_PERIODIC";
constexpr char CONST_VALUE_STR[]      = "BCT_CONSTANT_VALUE";
constexpr char CONST_DERIVATIVE_STR[] = "BCT_CONSTANT_DERIVATIVE";

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
