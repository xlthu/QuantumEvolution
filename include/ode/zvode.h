#ifndef _ZVODE_SOLVER_H
#define _ZVODE_SOLVER_H

#include "ode/ode.h"
#include "zvode/zvode.h"

class ZVODESolver : public ODESolver {
public:
    explicit ZVODESolver(StatePool& pool) : pool(pool) {}

    void solve(ODE* ode, State& psi1, double t1, double t2) override;
    void init_one_step(ODE* ode, const State& psi1, double t1, double t2) override;
    double solve_one_step(ODE* ode, State& psi1, double t1, double t2) override;

    enum class Method {
        Adams = 1,
        BDF = 2,
    };

    enum class Jacobian {
        No = 0,
        GeneratedFull = 2,
        GeneratedDiagonal = 3,
    };

    ZVODESolver& set_atol(double atol) { this->atol = atol; return *this; }
    ZVODESolver& set_rtol(double rtol) { this->rtol = rtol; return *this; }
    ZVODESolver& set_suggested_first_step_size(double h_suggested) { this->h_suggested = h_suggested; return *this; }
    ZVODESolver& set_min_step_size(double h_min) { this->h_min = h_min; return *this; }
    ZVODESolver& set_max_step_size(double h_max) { this->h_max = h_max; return *this; }
    ZVODESolver& set_max_nsteps(FINT max_nsteps) { this->max_nsteps = max_nsteps; return *this; }

    ZVODESolver& set_method(Method method) { this->method = method; return *this; }
    ZVODESolver& set_max_order(FINT max_order) {
        Assert_msg(
            ((method == Method::Adams) && max_order <= 12) || max_order <= 5,
            "Max order <= 12 for Adams, <=5 for BDF"
        );
        this->max_order = max_order;
        return *this;
    }
    ZVODESolver& set_jacobian(Jacobian jacobian) { this->jacobian = jacobian; return *this; }

private:
    double atol = 1e-8;
    double rtol = 1e-6;
    double h_suggested = 0.0;
    double h_min = 0.0;
    double h_max = 0.0;
    FINT max_nsteps = 10000;

    Method method = Method::Adams;
    FINT max_order = 12;
    Jacobian jacobian = Jacobian::No;

    StatePool& pool;

    // Workspace and settings
    FINT calc_mf();

    std::vector<std::complex<double>> zwork;
    FINT calc_lzw(FINT neq);

    std::vector<double> rwork;
    FINT calc_lrw(FINT neq);

    std::vector<FINT> iwork;
    FINT calc_liw(FINT neq);

    // Error
    void on_error(FINT istate);

    // For one step
    FINT neq;
    FINT itol;
    FINT itask;
    FINT istate;
    FINT iopt;
    FINT lzw;
    FINT lrw;
    FINT liw;
    FINT mf;
};

#endif // _ZVODE_SOLVER_H
