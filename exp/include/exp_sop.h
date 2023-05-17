#ifndef _EXP_SOP_H
#define _EXP_SOP_H

#include <vector>
#include "op/sop.h"

SOp embed(sz_t N, const std::vector<SOp>& ops, const std::vector<sz_t>& targets);

SOp sumLdagL(const std::vector<SOp>& L, const std::vector<SOp>& Ldag);

SOp rz(double theta);

SOp rx(double theta);

SOp ry(double theta);

#endif // _EXP_SOP_H
