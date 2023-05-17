#ifndef _ZVODE_H
#define _ZVODE_H

#include <complex>

#ifdef USE_64_INT
    #include <stdint.h>
    using FINT = int64_t;
#else // !USE_64_INT
    using FINT = int;
#endif // USE_64_INT

using zvode_f = void (*) (FINT* neq, double* t, std::complex<double>* y, std::complex<double>* ydot, void* rpar, void* ipar);

using zvode_jac = void (*) (FINT* neq, double* t, std::complex<double>* y, FINT* ml, FINT* mu, std::complex<double>* pd, FINT* nrowpd, void* rpar, void* ipar);

extern "C" void zvode_(zvode_f f, FINT* neq, std::complex<double>* y, double* t, double* tout,
                       FINT* itol, double* rtol, double* atol,
                       FINT* itask, FINT* istate, FINT* iopt,
                       std::complex<double>* zwork, FINT* lzw, double* rwork, FINT* lrw, FINT* iwork, FINT* liw,
                       zvode_jac jac, FINT* mf, void* rpar, void* ipar);

#endif // _ZVODE_H