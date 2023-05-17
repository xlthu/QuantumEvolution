#include "unraveling/qsd.h"

#include <chrono>

void QSD::derivative(State& dpsi, const State& psi, double t) {
    H->apply(dpsi, psi, t); // dpsi = -i H(t) |psi>
    dpsi *= _MI;

    if (n_lindblads()) {
        auto tmp1_g = pool.allocate_similar(psi);
        State& tmp1 = tmp1_g.state;

        auto tmp2_g = pool.allocate_similar(psi);
        State& tmp2 = tmp2_g.state;

        double sum = 0;
        for (std::size_t i = 0; i < n_lindblads(); ++i) {
            L[i]->apply(tmp1, psi, t); // tmp1 = L |psi>
            Complex e = psi.inner(tmp1); // e = <L> = <psi| L |psi>

            dpsi.axpy(conj(e), tmp1); // dpsi += <L>^ L |psi>

            Ldag[i]->apply(tmp2, tmp1, t); // tmp2 = L^ L |psi>
            dpsi.axpy(-0.5, tmp2); // dpsi += -1/2 L^ L |psi>

            sum += norm(e); // sum += <L>^ <L>
        }

        dpsi.axpy(-0.5 * sum, psi); // dpsi += -1/2 sum |psi>
    }
}

bool QSD::apply_stochastic(State& y_t_h, const State& y, double t, double h) {
    if (!n_lindblads()) return false;

    auto tmp1_g = pool.allocate_similar(y);
    State& tmp1 = tmp1_g.state;

    double sqrt_h = std::sqrt(h);
    for (std::size_t i = 0; i < n_lindblads(); ++i) {
        L[i]->apply(tmp1, y, t); // tmp1 = L |psi>
        Complex e = y.inner(tmp1); // e = <L> = <psi| L |psi>

        Complex dxi = sqrt_h * rnd();

        y_t_h.axpy(dxi, tmp1); // |dpsi> += L |psi> dxi
        y_t_h.axpy(-dxi * e, y); // |dpsi> += -<L> |psi> dxi
    }

    y_t_h.normalize();

    return true;
}

void QSD::solve(ODESolver* solver, State& psi1, double t1, double t2) {
    if (!n_lindblads()) {
        solver->solve(this, psi1, t1, t2);
        return;
    }

    auto psi_last_g = pool.allocate_similar(psi1);
    State& psi_last = psi_last_g.state;

    if (h_stoch > t2 - t1) Error("The stochastic step size is too large.");
    if (h_stoch == 0) Error("Evolve with linblads, but the stochastic step size is not set!");

    while (t1 < t2) {
        // One step
        psi_last = psi1;
        solver->solve(this, psi1, t1, t1 + h_stoch);
        // stochastic
        apply_stochastic(psi1, psi_last, t1, h_stoch);
        t1 += h_stoch;
    }
}

// QSDOpt

void QSDOpt::derivative(State& dpsi, const State& psi, double t) {
    H->apply(dpsi, psi, t); // dpsi = -i H(t) |psi>
    dpsi *= _MI;

    if (n_lindblads()) {
        auto tmp1_g = pool.allocate_similar(psi);
        State& tmp1 = tmp1_g.state;

        double sum = 0;
        for (std::size_t i = 0; i < n_lindblads(); ++i) {
            L[i]->apply(tmp1, psi, t); // tmp1 = L |psi>
            Complex e = psi.inner(tmp1); // e = <L> = <psi| L |psi>

            dpsi.axpy(conj(e), tmp1); // dpsi += <L>^ L |psi>

            sum += norm(e); // sum += <L>^ <L>
        }

        sum_LdagL.axpy_apply(dpsi, -0.5, psi, t); // dpsi += -1/2 sum L^L |psi>
        dpsi.axpy(-0.5 * sum, psi); // dpsi += -1/2 sum |psi>
    }
}