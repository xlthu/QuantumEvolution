#ifndef _EXP_STATE_H
#define _EXP_STATE_H

#include "exp_types.h"
#include "state/state.h"

inline dm_t ket2dm(const ket_t& ket) {
    bra_t bra = ctrans(ket);
    return ket * bra;
}

inline double state_fidelity(const State& a, const State& b) {
    return abs(a.inner(b));
}

double trace_distance(const dm_t& A, const dm_t& B);

double dm_fidelity(const dm_t& A, const dm_t& B);

template<typename Container>
inline sz_t bit_rep(sz_t n_qubits, const Container& qubits) {
    sz_t b = 0;
    for (auto q : qubits)
        b |= (1 << (n_qubits - 1 - i));
    return b;
}

#endif // _EXP_STATE_H
