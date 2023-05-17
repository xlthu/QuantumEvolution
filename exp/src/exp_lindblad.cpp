#include "exp_lindblad.h"

#include "op/sop.h"
#include "op/lindblad.h"

Lindblad amplitude_damping(sz_t target, double T1) {
    double c = 1. / sqrt(T1);

    return {
        .L = new ScaledSigmaPlus<double>(target, c),
        .Ldag = new ScaledSigmaMinus<double>(target, c),
    };
}

Lindblad amplitude_damping_s(double T1) {
    double c = 1. / sqrt(T1);

    return {
        .L = new SOp(c * SSigmaPlus),
        .Ldag = new SOp(c * SSigmaMinus),
    };
}

Lindblad phase_damping(sz_t target, double T2) {
    double c = 1. / sqrt(T2 * 2);

    return {
        .L = new ScaledSigmaZ<double>(target, c),
        .Ldag = new ScaledSigmaZ<double>(target, c),
    };
}

Lindblad phase_damping_s(double T2) {
    double c = 1. / sqrt(T2 * 2);

    return {
        .L = new SOp(c * SSigmaZ),
        .Ldag = new SOp(c * SSigmaZ),
    };
}

std::vector<Lindblad> relaxation(sz_t target, double T1, double T2) {
    if (T1 && T2) {
        Assert(2 * T1 >= T2);
        T2 = 1. / (1. / T2 - 1. / 2. / T1);
        return {amplitude_damping(target, T1), phase_damping(target, T2)};
    } else if (T1) {
        return {amplitude_damping(target, T1)};
    } else if (T2) {
        return {phase_damping(target, T2)};
    }

    return {};
}

std::vector<Lindblad> relaxation_s(double T1, double T2) {
    if (T1 && T2) {
        Assert(2 * T1 >= T2);
        T2 = 1. / (1. / T2 - 1. / 2. / T1);
        return {amplitude_damping_s(T1), phase_damping_s(T2)};
    } else if (T1) {
        return {amplitude_damping_s(T1)};
    } else if (T2) {
        return {phase_damping_s(T2)};
    }

    return {};
}

std::vector<Lindblad> depolaring(sz_t target, double gamma) {
    double c = sqrt(gamma);

    return {
        {
            .L = new ScaledSigmaX<double>(target, c),
            .Ldag = new ScaledSigmaX<double>(target, c),
        },
        {
            .L = new ScaledSigmaY<double>(target, c),
            .Ldag = new ScaledSigmaY<double>(target, c),
        },
        {
            .L = new ScaledSigmaZ<double>(target, c),
            .Ldag = new ScaledSigmaZ<double>(target, c),
        },
    };
}

std::vector<Lindblad> depolaring_s(double gamma) {
    double c = sqrt(gamma);

    return {
        {
            .L = new SOp(c * SSigmaX),
            .Ldag = new SOp(c * SSigmaX),
        },
        {
            .L = new SOp(c * SSigmaY),
            .Ldag = new SOp(c * SSigmaY),
        },
        {
            .L = new SOp(c * SSigmaZ),
            .Ldag = new SOp(c * SSigmaZ),
        },
    };
}
