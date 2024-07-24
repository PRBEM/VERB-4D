#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Matrix.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {
template <typename T>
struct type_caster<Matrix1D<T>> {
   public:
    PYBIND11_TYPE_CASTER(Matrix1D<T>, _("Matrix1D<T>"));

    // Conversion part 1 (Python -> C++)
    bool load(py::handle src, bool convert) {
        if (!convert and !py::array_t<T>::check_(src))
            return false;

        auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
        if (!buf)
            return false;

        auto dims = buf.ndim();
        if (dims != 1)
            return false;

        value = Matrix1D<T>(buf.shape()[0], buf.data());

        return true;
    }

    // Conversion part 2 (C++ -> Python)
    static py::handle cast(const Matrix1D<T>& src, py::return_value_policy policy, py::handle parent) {
        std::vector<size_t> shape(1);
        std::vector<size_t> strides(1);

        shape[0] = src.size_q1;
        strides[0] = sizeof(T);

        py::array a(std::move(shape), std::move(strides), src.matrix_array);

        return a.release();
    }
};

template <typename T>
struct type_caster<Matrix2D<T>> {
   public:
    PYBIND11_TYPE_CASTER(Matrix2D<T>, _("Matrix2D<T>"));

    // Conversion part 1 (Python -> C++)
    bool load(py::handle src, bool convert) {
        if (!convert and !py::array_t<T>::check_(src))
            return false;

        auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
        if (!buf)
            return false;

        auto dims = buf.ndim();
        if (dims != 2)
            return false;

        value = Matrix2D<T>(buf.shape(), buf.data());

        return true;
    }

    // Conversion part 2 (C++ -> Python)
    static py::handle cast(const Matrix2D<T>& src, py::return_value_policy policy, py::handle parent) {
        std::vector<size_t> shape(2);
        std::vector<size_t> strides(2);

        shape[0] = src.size_q1;
        shape[1] = src.size_q2;
        strides[0] = src.size_q2 * sizeof(T);
        strides[1] = sizeof(T);

        py::array a(std::move(shape), std::move(strides), src.plane_array);

        return a.release();
    }
};

}  // namespace detail
}  // namespace pybind11