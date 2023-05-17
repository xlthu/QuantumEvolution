#ifndef _RK45_H
#define _RK45_H

#include <limits>
#include "ode/ode.h"

class RK45Solver : public ODESolver {
public:
    explicit RK45Solver(StatePool& pool) : pool(pool) {}

    void solve(ODE* ode, State& psi1, double t1, double t2) override;

    void init_one_step(ODE* ode, const State& psi1, double t1, double t2) override;
    double solve_one_step(ODE* ode, State& psi1, double t1, double t2) override;

    RK45Solver& set_atol(double atol) { this->atol = atol; return *this; }
    RK45Solver& set_suggested_first_step_size(double h_suggested) { this->h_suggested = h_suggested; return *this; }
    RK45Solver& set_min_step_size(double h_min) { this->h_min = h_min; return *this; }
    RK45Solver& set_max_step_size(double h_max) { this->h_max = h_max; return *this; }
    RK45Solver& set_max_nsteps(std::size_t max_nsteps) { this->max_nsteps = max_nsteps; return *this; }

private:
    double atol = 1e-8;
    double h_suggested = 0.01;
    double h_min = 0.0;
    double h_max = std::numeric_limits<double>::infinity();
    std::size_t max_nsteps = 10000;

    StatePool& pool;

    void rkqs(State& yout, const State& y, double t, double t_crtl, double& h, double& hnext, ODE* ode);
    void rkck(State& yout, State& yerr, const State& y, double t, double h, ODE* ode);

    // For one step
    double h_next;
};

#endif // _RK45_H
