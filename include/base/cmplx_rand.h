#ifndef _CMPLX_RAND_H
#define _CMPLX_RAND_H

#include <random>
#include <cmath>

#include "base/types.h"

class CNormalRand {
public:
    using seed_t = typename std::mt19937::result_type;
    explicit CNormalRand(seed_t seed) : eng(seed), dis(0.0, M_SQRT1_2) {}
    Complex operator()() { return Complex{dis(eng), dis(eng)}; }
private:
    std::mt19937 eng;
    std::normal_distribution<double> dis;
};

#endif // _CMPLX_RAND_H