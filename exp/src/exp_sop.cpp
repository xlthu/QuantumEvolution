#include "exp_sop.h"

SOp embed(sz_t N, const std::vector<SOp>& ops, const std::vector<sz_t>& targets) {
    Assert(ops.size() == targets.size());
    std::vector<const SOp*> ids(N, &SId2);

    for (std::size_t i = 0; i < ops.size(); ++i) {
        ids[targets[i]] = &ops[i];
    }

    return SOp(ids);
}

SOp sumLdagL(const std::vector<SOp>& L, const std::vector<SOp>& Ldag) {
    Assert(L.size() == Ldag.size());

    SOp s = Ldag[0] * L[0];
    for (sz_t i = 1; i < L.size(); ++i) {
        s += Ldag[i] * L[i];
    }
    return s;
}

SOp rz(double theta) {
    theta /= 2.;
    return SOp{
        {
            {Complex{cos(-theta), sin(-theta)}, 0},
            {0, Complex{cos(theta), sin(theta)}},
        },
        SOpProperty{
            .symmetric = Symmetric::Normal,
            .hermitian = false,
            .diagonal = true,
            .triangular = Triangular::Not
        }
    };
}

SOp rx(double theta) {
    theta /= 2.;
    return SOp{
        {
            {Complex{cos(theta), 0}, Complex{0, -sin(theta)}},
            {Complex{0, -sin(theta)}, Complex{cos(theta), 0}},
        },
        SOpProperty{
            .symmetric = Symmetric::Normal,
            .hermitian = true,
            .diagonal = false,
            .triangular = Triangular::Not
        }
    };
}

SOp ry(double theta) {
    theta /= 2.;
    return SOp{
        {
            {Complex{cos(theta), 0}, Complex{-sin(theta), 0}},
            {Complex{sin(theta), 0}, Complex{cos(theta), 0}},
        },
        SOpProperty{
            .symmetric = Symmetric::Not,
            .hermitian = true,
            .diagonal = false,
            .triangular = Triangular::Not
        }
    };
}
