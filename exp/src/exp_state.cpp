#include "exp_state.h"

#include "base/assertion.h"

static dm_t sqrt_dm(const dm_t& dm) {
    Assert(dm.rows() == dm.columns());
    auto n = dm.rows();

    blaze::DynamicVector<Complex> w(n);
    blaze::DynamicMatrix<Complex> V(n, n);

    blaze::eigen(dm, w, V);

    w = blaze::sqrt(w);
    
    blaze::DynamicMatrix<Complex> W(n, n);
    blaze::diagonal(W) = w;

    return blaze::inv(V) * W * V;
}
double trace_distance(const dm_t& A, const dm_t& B) {
    dm_t M = A - B;
    auto MM = ctrans(M) * M;
    auto e = sqrt(abs(eigen(MM)));
    return 0.5 * sum(e);
}

double dm_fidelity(const dm_t& A, const dm_t& B) {
    auto sqrt_A = sqrt_dm(A);
    return sum(sqrt(abs(eigen(sqrt_A * B * sqrt_A))));
}
