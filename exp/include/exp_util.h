#ifndef _EXP_UTIL_H
#define _EXP_UTIL_H

#include <numeric>

template<typename T>
bool almost_equal(T a, T b, T tol = std::numeric_limits<T>::epsilon()) {
    return abs(a - b) <= tol;
}

#endif // _EXP_UTIL_H
