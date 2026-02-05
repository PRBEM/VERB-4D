#include <functional>

template <typename T, typename R, typename... Args>
auto const_reference_wrapper(R (*fun_to_wrap)(T, Args...)) {
    return [fun_to_wrap](T x, Args... args) {
        fun_to_wrap(x, args...);
        return x;
    };
}
