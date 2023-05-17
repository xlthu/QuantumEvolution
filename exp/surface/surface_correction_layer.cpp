#include "surface_correction_layer.h"

#include "exp_print.h"

void CorrectionLayer::init(Sys* sys) {
    // Sort
    sort(x.begin(), x.end(), [](sz_t a, sz_t b) { return a > b; });
    sort(id.begin(), id.end(), [](sz_t a, sz_t b) { return a > b; });

    // Sys
    this->sys = sys;
}

void CorrectionLayer::diagnose(std::ostream& out, const std::string& indent) {
    out << indent << "Type: " << (is_Z_only() ? "Z" : "X") << std::endl;
    out << indent << " X: " << x << std::endl;
    out << indent << " id: " << id << std::endl;
    out << indent << " Z: " << z << std::endl;
}

void CorrectionLayer::apply_layer(State& s, ODESolver* solver, bool noise_free) {
    if (is_Z_only()) {
        for (auto q : z) sys->sz.on(q).inplace_apply(s, 0);

    } else {
        if (noise_free) {
            for (auto q : x) sys->sx.on(q).inplace_apply(s, 0);

        } else {
            sys->unr.set_H(this);
            sys->unr.solve(solver, s, 0, sys->T_single);
        }
    }
}

void CorrectionLayer::apply(State& out, const State& s, double t) {
    out = 0;

    auto tmpg = sys->pool.allocate_similar(s);
    State& tmp = tmpg.state;

    // X
    double p_r180 = sys->p_r180(t);
    for (auto q : x) {
        sys->sx.on(q).apply(tmp, s, t);
        out.axpy(p_r180, tmp);
    }

    // ID
    double p_id = sys->p_id(t);
    for (auto q : id) {
        sys->sx.on(q).apply(tmp, s, t);
        out.axpy(p_id, tmp);
    }
}
