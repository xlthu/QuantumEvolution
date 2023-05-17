#ifndef _ODE_SOLVER_H
#define _ODE_SOLVER_H

#include "state/state.h"
#include "state/state_pool.h"

class ODE {
public:
    virtual ~ODE() = default;
    virtual void derivative(State& dy, const State& y, double t) = 0;
};

class ODESolver {
public:
    virtual ~ODESolver() = default;
    // Evolve psi1 according ode from t1 to t2
    virtual void solve(ODE* ode, State& psi1, double t1, double t2) = 0;

    virtual void init_one_step(ODE* ode, const State& psi1, double t1, double t2) = 0;
    virtual double solve_one_step(ODE* ode, State& psi1, double t1, double t2) = 0;
};

#endif // _ODE_SOLVER_H