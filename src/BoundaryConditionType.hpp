#ifndef BOUNDARY_CONDITION_TYPE_H_
#define BOUNDARY_CONDITION_TYPE_H_

#include <string>
#include <stdexcept>

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

inline std::stringstream& operator>>(std::stringstream& ss, BoundaryConditionType& type)
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
#endif
