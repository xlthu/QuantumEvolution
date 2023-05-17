#ifndef _EXP_LINDBLAD_H
#define _EXP_LINDBLAD_H

#include <vector>
#include "base/types.h"

class Op;

struct Lindblad {
    Op* L;
    Op* Ldag;
};

Lindblad amplitude_damping(sz_t target, double T1);

Lindblad amplitude_damping_s(double T1);

Lindblad phase_damping(sz_t target, double T2);

Lindblad phase_damping_s(double T2);

std::vector<Lindblad> relaxation(sz_t target, double T1, double T2);

std::vector<Lindblad> relaxation_s(double T1, double T2);

std::vector<Lindblad> depolaring(sz_t target, double gamma);

std::vector<Lindblad> depolaring_s(double gamma);

#endif // _EXP_LINDBLAD_H
