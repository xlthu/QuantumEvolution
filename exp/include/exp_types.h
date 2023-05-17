#ifndef _EXP_TYPES_H
#define _EXP_TYPES_H

#include <blaze/Blaze.h>

#include "base/types.h"

using ket_t = blaze::DynamicVector<Complex,blaze::columnVector>;
using bra_t = blaze::DynamicVector<Complex,blaze::rowVector>;
using dm_t = blaze::DynamicMatrix<Complex>;

#endif // _EXP_TYPES_H
