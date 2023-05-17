#ifndef _JUMP_H
#define _JUMP_H

#include <random>

#include "base/assertion.h"
#include "op/op.h"
#include "op/sop.h"
#include "unraveling/unraveling.h"

class JumpInfo {
public:
    double time;
    std::size_t idx;
};

class Jump : public Unraveling {
public:
    Jump(Op* H, std::vector<Op*> L, std::vector<Op*> Ldag, StatePool& pool, std::mt19937& eng)
        : pool(pool), eng(eng), rnd(0.0, 1.0) {
        set_H(H);
        set_lindblads(std::move(L), std::move(Ldag));
        new_trajectory();
    }

    std::size_t n_lindblads() const { return L.size(); }
    const std::vector<JumpInfo>& jumps() const { return jump_info; }
    void clear_jumps() { jump_info.clear(); }

    void solve(ODESolver* solver, State& psi1, double t1, double t2) override;

    void derivative(State& dy, const State& y, double t) override;

    void set_norm_time_atol(double norm_time_atol) { this->norm_time_atol = norm_time_atol; }
    void set_norm_rtol(double norm_rtol) { this->norm_rtol = norm_rtol; }
    void set_max_norm_nsteps(std::size_t max_norm_nsteps) { this->max_norm_nsteps = max_norm_nsteps; }

    void set_H(Op* H) { this->H = H; }
    void set_lindblads(std::vector<Op*> L, std::vector<Op*> Ldag) {
        Assert(L.size() == Ldag.size());
        this->L = std::move(L); this->Ldag = std::move(Ldag);
        cum_probs.resize(this->L.size());
    }

    void new_trajectory() override {
        target_norm2 = rnd(eng);
        clear_jumps();
    }

protected:
    Op* H;
    std::vector<Op*> L;
    std::vector<Op*> Ldag;

    StatePool& pool;

    std::mt19937& eng;
    std::uniform_real_distribution<double> rnd;

    void locate_jump_time(State& psi_prev, double t_prev, double norm2_prev, State& psi, double& t, double norm2_now, double target_norm2, ODESolver* solver);
    double norm_time_atol = 1e-6;
    double norm_rtol = 1e-3;
    std::size_t max_norm_nsteps = 5;
    
    void jump(State& psi, double t);
    double target_norm2;
    std::vector<double> cum_probs;

    std::vector<JumpInfo> jump_info;
};

class JumpOpt : public Jump {
public:
    JumpOpt(Op* H, std::vector<Op*> L, std::vector<Op*> Ldag, SOp sum_LdagL, StatePool& pool, std::mt19937& eng) : Jump(H, std::move(L), std::move(Ldag), pool, eng), sum_LdagL(std::move(sum_LdagL)) {}

    void derivative(State& dy, const State& y, double t) override;

    void set_lindblads(std::vector<Op*> L, std::vector<Op*> Ldag, SOp sum_LdagL) {
        set_lindblads(std::move(L), std::move(Ldag));
        this->sum_LdagL = std::move(sum_LdagL);
    }

protected:
    SOp sum_LdagL;

    using Jump::set_lindblads;
};

#endif // _JUMP_H
