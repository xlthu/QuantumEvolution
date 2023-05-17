#ifndef _BLAS_H
#define _BLAS_H

#include "base/types.h"
#include <mkl.h>
static_assert(sizeof(MKL_INT) == sizeof(sz_t));

inline void blas_init(int num_threads = 1) {
    mkl_set_num_threads(num_threads);
}

#endif // _BLAS_H
