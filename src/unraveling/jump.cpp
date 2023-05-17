#include "unraveling/jump.h"

#include <cmath>

template<typename T> T pow2(T x) { return x * x; }

void Jump::derivative(State& dy, const State& y, double t) {
    H->apply(dy, y, t); // dpsi = -i H(t) |psi>
    dy *= _MI;

    if (n_lindblads()) {
        auto tmp1_g = pool.allocate_similar(y);
        State& tmp1 = tmp1_g.state;

        auto tmp2_g = pool.allocate_similar(y);
        State& tmp2 = tmp2_g.state;

        for (std::size_t i = 0; i < n_lindblads(); ++i) {
            L[i]->apply(tmp1, y, t); // tmp1 = L |psi>
            Ldag[i]->apply(tmp2, tmp1, t); // tmp2 = L^ L |psi>

            dy.axpy(-0.5, tmp2); // dpsi += -1/2 L^ L |psi>
        }
    }
}

void Jump::solve(ODESolver* solver, State& psi, double t1, double t2) {
    if (!n_lindblads()) {
        solver->solve(this, psi, t1, t2);
        return;
    }

    auto psi_prev_g = pool.allocate_similar(psi);
    State& psi_prev = psi_prev_g.state;

    solver->init_one_step(this, psi, t1, t2);
    double target_norm2 = rnd(eng);

    double norm2_prev = pow2(psi.norm());
    while (t1 < t2) {
        psi_prev = psi;

        double h = solver->solve_one_step(this, psi, t1, t2);
        double norm2_now = pow2(psi.norm());
        // std::cout << norm2_now << std::endl;
        t1 += h;
        // psi_prev (at time t1 - h, with norm2_prev) -> psi (at time t1, with norm2_now)

        if (norm2_now <= target_norm2) { // Jump
            locate_jump_time(psi_prev, t1 - h, norm2_prev, psi, t1, norm2_now, target_norm2, solver);
            // Now, we have psi at time t1 with target_norm2 waiting to jump

            // std::cout << "jump" << std::endl;
            jump(psi, t1);
            norm2_now = 1.;
            
            // Re-init
            solver->init_one_step(this, psi, t1, t2);
            target_norm2 = rnd(eng);
        }

        norm2_prev = norm2_now;
    }

    psi.normalize();
}

void Jump::locate_jump_time(State& psi_prev, double t_prev, double norm2_prev, State& psi, double& t, double norm2_now, double target_norm2, ODESolver* solver) {
    double t_final = t;

    for (std::size_t step = 0; step < max_norm_nsteps; ++step) {
        if (t_final - t_prev < norm_time_atol) return;

        double t_guess = t_prev +
                         log(norm2_prev / target_norm2) / log(norm2_prev / norm2_now) * (t_final - t_prev);
        if (t_guess < t_prev + norm_time_atol) t_guess = t_prev + norm_time_atol;

        psi = psi_prev;
        solver->solve(this, psi, t_prev, t_guess);
        t = t_guess;

        double norm2_guess = pow2(psi.norm());

        if (abs(target_norm2 - norm2_guess) < norm_rtol * target_norm2) {
            return;

        } else if (norm2_guess < target_norm2) {
            // t_guess is still > t_jump
            t_final = t_guess;
            norm2_now = norm2_guess;
            
        } else {
            // t_guess < t_jump
            t_prev = t_guess;
            psi_prev = psi;
            norm2_prev = norm2_guess;
        }
    }

    Error("Norm tolerance for jumping is not reached." << Error_endl
          << "Increase accuracy of ODE solver or max. norm nsteps."
         );
}

void Jump::jump(State& psi, double t) {
    auto tmp_g = pool.allocate_similar(psi);
    State& tmp = tmp_g.state;

    // Calc probs
    double sum = 0;
    for (std::size_t i = 0; i < L.size(); ++i) {
        L[i]->apply(tmp, psi, t); // tmp = L |psi>
        sum += pow2(tmp.norm());
        cum_probs[i] = sum;
    }

    // Determine which to jump
    double r = rnd(eng) * sum;
    std::size_t i = 0;
    for (; i < L.size(); ++i) {
        if (cum_probs[i] > r) break;
    }

    // Jump
    tmp = psi;
    L[i]->apply(psi, tmp, t);
    psi.normalize();

    // Record
    jump_info.push_back(JumpInfo{
        .time = t,
        .idx = i,
    });
}


// JumpOpt

void JumpOpt::derivative(State& dy, const State& y, double t) {
    H->apply(dy, y, t); // dpsi = -i H(t) |psi>
    dy *= _MI;

    if (n_lindblads()) {
        sum_LdagL.axpy_apply(dy, -0.5, y, t);
    }
}