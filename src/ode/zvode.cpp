#include "ode/zvode.h"
#include "unraveling/unraveling.h"
#include <cstring>

static void _zvode_f(FINT* neq, double* t, std::complex<double>* y, std::complex<double>* ydot, void* rpar, void* ipar) {
    ODE* ode = static_cast<ODE*>(rpar);
    State* for_shadow = static_cast<State*>(ipar);

    State sy = State::shadow(y, *for_shadow);
    State sydot = State::shadow(ydot, *for_shadow);

    ode->derivative(sydot, sy, *t);
}

void ZVODESolver::solve(ODE* ode, State& psi1, double t1, double t2) {
    init_one_step(ode, psi1, t1, t2);
    itask = 1; // 1: calc up to t2
    solve_one_step(ode, psi1, t1, t2);
}

void ZVODESolver::init_one_step(ODE* ode, const State& psi1, double t1, double t2) {
    neq = psi1.total_dims();
    itol = 1; // Scalar atol and scalar rtol

    itask = 5; // Task type:
    // 1: calc up to t2
    // 2: calc one step
    // 5: calc one step without passing TCRIT

    istate = 1; // State:
    // 1: to init
    // 2: continue

    iopt = 1; // Optional input is used

    // Workspace and settings
    lzw = calc_lzw(neq);
    zwork.resize(lzw);

    lrw = calc_lrw(neq);
    rwork.resize(lrw);
    memset(rwork.data(), 0, sizeof(double) * 20); // std::fill(rwork.begin(), rwork.begin() + 20, 0.);
    rwork[0] = t2; // TCRIT, the time that the sovler will not pass
    rwork[4] = h_suggested;
    rwork[5] = h_max;
    rwork[6] = h_min;

    liw = calc_liw(neq);
    iwork.resize(liw);
    memset(iwork.data(), 0, sizeof(FINT) * 30); // std::fill(iwork.begin(), iwork.begin() + 30, 0);
    iwork[4] = max_order;
    iwork[5] = max_nsteps;
    // iwork[6] = 0; // Maximum number of messages printed -- Keep default

    mf = calc_mf();
}

double ZVODESolver::solve_one_step(ODE* ode, State& psi1, double t1, double t2) {
    zvode_(_zvode_f, &neq, psi1.data(), &t1, &t2,
           &itol, &rtol, &atol,
           &itask, &istate, &iopt,
           zwork.data(), &lzw, rwork.data(), &lrw, iwork.data(), &liw,
           nullptr, &mf, ode, &psi1
          );

    if (istate < 0) on_error(istate);

    return rwork[10];
}

FINT ZVODESolver::calc_mf() {
    return 10 * static_cast<FINT>(method) + static_cast<FINT>(jacobian);
}

FINT ZVODESolver::calc_lzw(FINT neq) {
    FINT lwm = 0;
    switch (jacobian) {
        case Jacobian::GeneratedFull: lwm = 2 * neq * neq; break;
        case Jacobian::GeneratedDiagonal: lwm = neq; break;
        default: break;
    }
    return neq * (max_order + 3) + lwm;
}

FINT ZVODESolver::calc_lrw(FINT neq) {
    return 20 + neq;
}

FINT ZVODESolver::calc_liw(FINT neq) {
    return (jacobian == Jacobian::No || jacobian == Jacobian::GeneratedDiagonal) ? 30 : 30 + neq;
}

void ZVODESolver::on_error(FINT istate) {
    switch (istate) {
        case -1:
            Error("[ZVODE] Max nsteps exceeded.");
            break;
        case -2:
            Error("[ZVODE] Too small atol/rtol.");
            break;
        case -3:
            Error("[ZVODE] Illegal input.");
            break;
        case -4:
            Error("[ZVODE] Repeated error test failures. (The problem may have a singularity, or the input may be inappropriate.)");
            break;
        case -5:
            Error("[ZVODE] Repeated convergence test failures. (Perhaps bad Jacobian supplied or wrong choice of MF or tolerances.)");
            break;
        case -6:
            Error("[ZVODE] Error weight became zero during the integration. (Solution component i vanished that ATOL or ATOL[i] = 0.)");
            break;
        default:
            Error("[ZVODE] Unknown state: " << istate);
            break;
    }
}
