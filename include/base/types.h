#ifndef _TYPES_H
#define _TYPES_H

#include <complex>

using Complex = std::complex<double>;

constexpr Complex _I{0.0, 1.0};
constexpr Complex _MI{0.0, -1.0};

#ifdef USE_64_INT
    using sz_t = unsigned long long;
    static_assert(sizeof(sz_t) == 8);
#else // !USE_64_INT
    using sz_t = unsigned int;
    static_assert(sizeof(sz_t) == 4);
#endif // USE_64_INT

using ssz_t = typename std::make_signed<sz_t>::type;

#endif // _CMPLX_H