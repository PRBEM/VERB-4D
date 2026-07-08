/*
 * SPDX-FileCopyrightText: 2025 GFZ Helmholtz Centre for Geosciences
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

#include <functional>

template <typename T, typename R, typename... Args>
auto const_reference_wrapper(R (*fun_to_wrap)(T, Args...)) {
    return [fun_to_wrap](T x, Args... args) {
        fun_to_wrap(x, args...);
        return x;
    };
}
