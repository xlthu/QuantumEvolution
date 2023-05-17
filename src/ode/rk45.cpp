#include "ode/rk45.h"

constexpr double SAFETY = 0.9;
constexpr double PGROW = -0.2;
constexpr double PSHRNK = -0.25;
constexpr double ERRCON = 1.89e-4;

void RK45Solver::rkck(State& yout, State& yerr, const State& y, double t, double h, ODE* ode) {
    // Given time t, state y(t) and step size h, calc result yout and the error yerr

    constexpr double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875, b21 = 0.2,
                     b31 = 3.0 / 8.0, b32 = 9.0 / 40.0, b41 = 4.0, b42 = -4.0, b43 = 1.2,
                     b51 = -55.0 / 81.0, b52 = -25.0 / 9.0, b53 = -175.0 / 81.0, b54 = 35.0 / 27.0,
                     b61 = -1631.0 / 11264.0, b62 = 35.0 / 256.0, b63 = -115.0 / 7168.0,
                     b64 = 1265.0 / 4096.0, b65 = 253.0 / 4096.0, c1 = 37888.0 / 11417.0,
                     c3 = 5120.0 / 529.0, c4 = 10240.0 / 19481.0, c6 = 512.0 / 1771.0,
                     dc5 = -554.0 / 1771.0;
    constexpr double dc1 = -831.0 / 18944.0, dc3 = 831.0 / 17920.0,
                     dc4 = -831.0 / 5120.0, dc6 = 277.0 / 2048.0;

    auto dydt_g = pool.allocate_similar(y);
    State& dydt = dydt_g.state;

    auto ak2_g = pool.allocate_similar(y);
    State& ak2 = ak2_g.state;

    auto ak3_g = pool.allocate_similar(y);
    State& ak3 = ak3_g.state;

    auto ak4_g = pool.allocate_similar(y);
    State& ak4 = ak4_g.state;

    auto ak5_g = pool.allocate_similar(y);
    State& ak5 = ak5_g.state;

    yout = y;
    ode->derivative(dydt, yout, t);             // First step
    yout.axpy(h * b21, dydt);

    ode->derivative(ak2, yout, t + a2 * h);         // Second step
    yout = y;
    yout.axpy(h * b21 * b31, dydt);
    yout.axpy(h * b32, ak2);

    ode->derivative(ak3, yout, t + a3 * h);         // Third step
    yout = y;
    yout.axpy(h * b21 * b31 * b41, dydt);
    yout.axpy(h * b32 * b42, ak2);
    yout.axpy(h * b43, ak3);

    ode->derivative(ak4, yout, t + a4 * h);         // Fourth step
    yout = y;
    yout.axpy(h * b21 * b31 * b41 * b51, dydt);
    yout.axpy(h * b32 * b42 * b52, ak2);
    yout.axpy(h * b43 * b53, ak3);
    yout.axpy(h * b54, ak4);

    ode->derivative(ak5, yout, t + a5 * h);         // Fifth step
    yout = y;
    yout.axpy(h * b21 * b31 * b41 * b51 * b61, dydt);
    yout.axpy(h * b32 * b42 * b52 * b62, ak2);
    yout.axpy(h * b43 * b53 * b63, ak3);
    yout.axpy(h * b54 * b64, ak4);
    yout.axpy(h * b65, ak5);

    ode->derivative(ak2, yout, t + a6 * h);         // Sixth step

    // Accumulate increments with proper weights
    yout = y;
    yout.axpy(h * b21 * b31 * b41 * b51 * b61 * c1, dydt);
    yout.axpy(h * b43 * b53 * b63 * c3, ak3);
    yout.axpy(h * b54 * b64 * c4, ak4);
    yout.axpy(h * c6, ak2);

    // Estimate error as difference between fourth and fifth order methods.
    yerr.mul(h * b21 * b31 * b41 * b51 * b61 * c1 * dc1, dydt);
    yerr.axpy(h * b43 * b53 * b63 * c3 * dc3, ak3);
    yerr.axpy(h * b54 * b64 * c4 * dc4, ak4);
    yerr.axpy(h * b65 * dc5, ak5);
    yerr.axpy(h * c6 * dc6, ak2);
}

void RK45Solver::rkqs(State& yout, const State& y, double t, double t_crtl, double& h, double& hnext, ODE* ode) {
    // One step
    // t: the time from which one step is taken
    // h: the suggested step size -> the actual step size
    // y: the state at time t
    // yout : the state at time t + h
    // hnext: only used to record the next suggested step size

    auto yerr_g = pool.allocate_similar(y);
    State& yerr = yerr_g.state;

    // Adjust h
    if (h > h_max) h = h_max;
    if (h < h_min) h = h_min;
    auto h_ctrl = t_crtl - t;
    if (h > h_ctrl) h = h_ctrl;

    while (true) {
        rkck(yout, yerr, y, t, h, ode);

        double errmax = yerr.norm() / atol;
        if (errmax > 1.0) {
            // Reject the tried step size h
            // Calculate the next trial step size and update h
            
            h = SAFETY * h * pow(errmax, PSHRNK);

            if (t + h == t) {
                Error("Stepsize underflow in rkqs" << Error_endl
                    << "    errmax = " << errmax << Error_endl
                    << "    h = " << h << Error_endl
                    << "    atol = " << atol
                );
            }
        } else {
            // Accept the step size h
            // Calculate the next suggested step size

            // hnext
            if (errmax > ERRCON) hnext = SAFETY * h * pow(errmax, PGROW);
            else hnext = 5.0 * h;

            break;
        }
    }
}

void RK45Solver::solve(ODE* ode, State& psi, double t1, double t2) {
    auto psi_tmp_g = pool.allocate_similar(psi);
    State& psi_tmp = psi_tmp_g.state;

    State* y = &psi;
    State* yout = &psi_tmp;

    init_one_step(ode, *y, t1, t2);

    for (std::size_t step = 1; step <= max_nsteps; ++step) {
        // Solve for the current step
        double h = h_next;
        rkqs(*yout, *y, t1, t2, h, h_next, ode);

        // Procceed to the next step
        t1 += h;
        std::swap(y, yout);

        if (t1 >= t2) {
            // Finished
            if (y != &psi) psi = *y;
            return;
        }
    }
    
    Error("Max. ode steps exceeded.");
}

void RK45Solver::init_one_step(ODE* ode, const State& psi1, double t1, double t2) {
    Assert(t2 > t1);

    h_next = h_suggested;
}

double RK45Solver::solve_one_step(ODE* ode, State& psi1, double t1, double t2) {
    auto tmp_g = pool.allocate_similar(psi1);
    State& tmp = tmp_g.state;

    tmp = psi1;
    double h = h_next;
    rkqs(psi1, tmp, t1, t2, h, h_next, ode);
    return h;
}
