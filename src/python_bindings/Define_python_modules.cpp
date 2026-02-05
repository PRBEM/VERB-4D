#include "BoundaryConditionType.hpp"
#include "Convection_2D.h"
#include "Diffusion_2D.h"
#include "Matrix_conversion.h"
#include "const_reference_wrapper.h"

namespace py = pybind11;

PYBIND11_MODULE(verb4d_solver, m) {
    py::enum_<BoundaryConditionType>(m, "BoundaryConditionType", py::arithmetic(), "Boundary condition type for diffusion and convection.")
        .value("Periodic", BoundaryConditionType::Periodic, "Periodic boundary condition")
        .value("ConstantValue", BoundaryConditionType::ConstantValue, "Constant value boundary condition")
        .value("ConstantDerivative", BoundaryConditionType::ConstantDerivative, "Constant derivative boundary condition")
        .export_values();

    m.doc() = "VERB4D Solver python bindings";

    m.def("Convection_1D", const_reference_wrapper(&Convection_1D_ULTIMATE_QUICKEST6));
    m.def("Convection_2D", const_reference_wrapper(&Convection_2D));
    m.def("Diffusion_2D", const_reference_wrapper(&Diffusion_2D));
}