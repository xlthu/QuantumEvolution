#ifndef _QSD_H
#define _QSD_H

#include "base/cmplx_rand.h"
#include "base/assertion.h"
#include "state/state_pool.h"
#include "op/op.h"
#include "op/sop.h"
#include "unraveling/unraveling.h"

class QSD : public Unraveling {
public:
    QSD(Op* H, std::vector<Op*> L, std::vector<Op*> Ldag, StatePool& pool, CNormalRand& rnd)
        : pool(pool), rnd(rnd) {
        set_H(H);
        set_lindblads(std::move(L), std::move(Ldag));
    }

    std::size_t n_lindblads() const { return L.size(); }

    void solve(ODESolver* solver, State& psi1, double t1, double t2) override;
    void derivative(State& dy, const State& y, double t) override;

    void set_stochastic_step_size(double h_stoch) { this->h_stoch = h_stoch; }

    void set_H(Op* H) { this->H = H; }
    void set_lindblads(std::vector<Op*> L, std::vector<Op*> Ldag) {
        Assert(L.size() == Ldag.size());
        this->L = std::move(L); this->Ldag = std::move(Ldag);
    }

    void new_trajectory() override {}

protected:
    double h_stoch = 0.01;

    Op* H;
    std::vector<Op*> L;
    std::vector<Op*> Ldag;

    StatePool& pool;

    CNormalRand& rnd;

    bool apply_stochastic(State& y_t_h, const State& y, double t, double h);
};

class QSDOpt : public QSD {
public:
    QSDOpt(Op* H, std::vector<Op*> L, std::vector<Op*> Ldag, SOp sum_LdagL, StatePool& pool, CNormalRand& rnd) : QSD(H, std::move(L), std::move(Ldag), pool, rnd), sum_LdagL(std::move(sum_LdagL)) {}

    void derivative(State& dy, const State& y, double t) override;

    void set_lindblads(std::vector<Op*> L, std::vector<Op*> Ldag, SOp sum_LdagL) {
        set_lindblads(std::move(L), std::move(Ldag));
        this->sum_LdagL = std::move(sum_LdagL);
    }

protected:
    SOp sum_LdagL;

    using QSD::set_lindblads;
};

#endif // _QSD_H
