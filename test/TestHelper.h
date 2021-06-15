#pragma once

#include "../src/Matrix.h"
#include <initializer_list>
#include <iostream>
#include <limits>
#include <sstream>
#include "test_runner.h"

#define PRINT(matrix) \
{\
    ostringstream os; \
    os << #matrix; \
    print((matrix), os.str()); \
}\

template <typename T>
void print(const Matrix2D<T>& m, const string& name) {
    std::cout << name << ':' << std::endl;
    std::cout << '{' << endl;
    for (auto i = 0; i < m.size_q1; ++i) {
        for (auto j = 0; j < m.size_q2; ++j) {
            if (j != 0) {
                std::cout << ' ';
            }
            std::cout << m[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << '}' << std::endl;
}

template <typename T>
Matrix1D<T> initialize(std::initializer_list<T> lst) {
    Matrix1D<T> result(lst.size());
    for (auto i = 0; i < lst.size(); ++i) {
        result[i] = *(lst.begin() + i);
    }
    return result;
}

template <typename T>
Matrix2D<T> initialize(
    std::initializer_list<std::initializer_list<T>> lst) {

    int size_q1 = lst.size();
    if (size_q1 == 0) {
        return {};
    }
    int size_q2 = lst.begin()->size();
    for (auto it = std::next(lst.begin()); it != lst.end(); ++it) {
        if (size_q2 != it->size()) {
            std::cout << 
                "Error! Cannot initialize Matrix2D. Wrong initializer list size." 
                << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    Matrix2D<T> result(size_q1, size_q2);
    for (auto i = 0; i < size_q1; ++i) {
        auto lst_1d = *(lst.begin() + i);
        for (auto j = 0; j < size_q2; ++j) { 
            result[i][j] = *(lst_1d.begin() + j);
        }
    }
    return result;
}

template <typename T>
Matrix3D<T> initialize(
    std::initializer_list<
        std::initializer_list<
            std::initializer_list<T>>> lst) {

    int size_q1 = lst.size();
    if (size_q1 == 0) {
        return {};
    }
    int size_q2 = lst.begin()->size();
    if (size_q2 == 0) {
        return {};
    }
    int size_q3 = lst.begin()->begin()->size();

    for (auto it = std::next(lst.begin()); it != lst.end(); ++it) {
        if (size_q2 != it->size()) {
            std::cout << 
                "Error! Cannot initialize Matrix3D. Wrong initializer list size." 
                << std::endl;
            exit(EXIT_FAILURE);
        }
        for (auto jt = std::next(it->begin()); jt != it->end(); ++jt) {
            if (size_q3 != jt->size()) {
                std::cout << 
                    "Error! Cannot initialize Matrix3D. Wrong initializer list size." 
                    << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    Matrix3D<T> result(size_q1, size_q2, size_q3);
    for (auto i = 0; i < size_q1; ++i) {
        auto lst_2d = *(lst.begin() + i);
        for (auto j = 0; j < size_q2; ++j) { 
            auto lst_1d = *(lst_2d.begin() + j);
            for (auto k = 0; k < size_q3; ++k) { 
                result[i][j][k] = *(lst_1d.begin() + k);
            }
        }
    }
    return result;
}

template <typename T>
void TestEqual(const Matrix1D<T>& a, const Matrix1D<T>& b) {
    ASSERT_EQUAL(a.size_q1, b.size_q1);
    for (auto i = 0; i < a.size_q1; ++i) {
        ASSERT_EQUAL(a[i], b[i]);
    }
}

template <typename T>
void TestEqual(const Matrix2D<T>& a, const Matrix2D<T>& b) {
    ASSERT_EQUAL(a.size_q1, b.size_q1);
    ASSERT_EQUAL(a.size_q2, b.size_q2);
    for (auto i = 0; i < a.size_q1; ++i) {
        for (auto j = 0; j < a.size_q2; ++j) {
            ASSERT_EQUAL(a[i][j], b[i][j]);
        }
    }
}

template <typename T>
void TestEqual(const Matrix3D<T>& a, const Matrix3D<T>& b) {
    ASSERT_EQUAL(a.size_q1, b.size_q1);
    ASSERT_EQUAL(a.size_q2, b.size_q2);
    ASSERT_EQUAL(a.size_q3, b.size_q3);
    for (auto i = 0; i < a.size_q1; ++i) {
        for (auto j = 0; j < a.size_q2; ++j) {
            for (auto k = 0; k < a.size_q3; ++k) {
                ASSERT_EQUAL(a[i][j][k], b[i][j][k]);
            }
        }
    }
}

template <typename T>
void TestEqualTol(T a, T b, T tol = std::numeric_limits<T>::epsilon() * 3
) {
    ASSERT(std::abs(a - b) < tol);
}

template <typename T>
void TestEqualTol(
    const Matrix2D<T>& a, const Matrix2D<T>& b, 
    T tol = std::numeric_limits<T>::epsilon() * 3
) {
    ASSERT_EQUAL(a.size_q1, b.size_q1);
    ASSERT_EQUAL(a.size_q2, b.size_q2);
    for (auto i = 0; i < a.size_q1; ++i) {
        for (auto j = 0; j < a.size_q2; ++j) {
            //ASSERT(std::abs(a[i][j] - b[i][j]) < tol);
            TestEqualTol(a[i][j], b[i][j], tol);
        }
    }
}

template <typename T>
void TestEqualTol(const Matrix1D<T>& a, const Matrix1D<T>& b,
    T tol = std::numeric_limits<T>::epsilon() * 3
) {
    ASSERT_EQUAL(a.size_q1, b.size_q1);
    for (auto i = 0; i < a.size_q1; ++i) {
        //ASSERT(std::abs(a[i] - b[i]) < tol);
        TestEqualTol(a[i], b[i], tol);
    }
}

template <typename T>
void TestEqualTol(const Matrix3D<T>& a, const Matrix3D<T>& b,
    T tol = std::numeric_limits<T>::epsilon() * 3
) {
    ASSERT_EQUAL(a.size_q1, b.size_q1);
    ASSERT_EQUAL(a.size_q2, b.size_q2);
    ASSERT_EQUAL(a.size_q3, b.size_q3);
    for (auto i = 0; i < a.size_q1; ++i) {
        for (auto j = 0; j < a.size_q2; ++j) {
            for (auto k = 0; k < a.size_q3; ++k) {
                //ASSERT(std::abs(a[i][j][k] - b[i][j][k]) < tol);
                TestEqualTol(a[i][j][k], b[i][j][k]);
            }
        }
    }
}

