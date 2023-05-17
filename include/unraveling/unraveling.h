#ifndef _UNRAVELING_H
#define _UNRAVELING_H

#include "state/state.h"
#include "ode/ode.h"

class Unraveling : public ODE {
public:
    virtual ~Unraveling() = default;
    // Evolve psi1 from t1 to t2
    virtual void solve(ODESolver* solver, State& psi1, double t1, double t2) = 0;
    
    virtual void new_trajectory() = 0;
};

#endif // _UNRAVELING_H