#ifndef _EXP_PULSE_H
#define _EXP_PULSE_H

#include "base/types.h"

template<typename T> T pow2(T x) { return x * x; }

template<sz_t N>
class FCos {
public:
    FCos(double T, double* params = nullptr) : T(T) {
        for (sz_t i = 1; i <= N; ++i) {
            omegas[i - 1] = 2 * M_PI * (double)i / T;
        }
        if (params) reset_params(params);
    }

    constexpr sz_t n_params() const {
        return N;
    }

    void reset_params(double* params) {
        for (sz_t i = 0; i < N; ++i) {
            coeff[i] = params[i] / 2.;
        }
    }

    double operator()(double t) const {
        double sum = 0;
        for (sz_t i = 0; i < N; ++i) {
            sum += coeff[i] * (1 - cos(omegas[i] * t));
        }
        return sum;
    }

    double coeff[N];
    double omegas[N];
    const double T;
};

class GaussianZero { // Not the same with Python, A in Python is double that in C++
public:
    GaussianZero(double T, double sigma, double A) 
        : T(T), sigma(sigma), A(A / 2),
         sigma2(sigma * sigma),
         alpha(exp(-T * T / 8 / sigma2)),
         beta(sqrt(2 * M_PI * sigma2) * erf(T / sqrt(8) / sigma)) {}

    constexpr sz_t n_params() const {
        return 1;
    }

    void reset_params(double* params) {
        A = params[0] / 2;
    }

    double operator()(double t) const {
        double g = exp(-pow2(t - T / 2) / 2 / sigma2);
        return A * (g - alpha) / (beta - T * alpha);
    }

    const double T;
    const double sigma;
    double A;

    const double sigma2;
    const double alpha;
    const double beta;
};

#endif // _EXP_PULSE_H
