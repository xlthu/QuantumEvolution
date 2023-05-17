#include "op/prim_op.h"

#include "base/assertion.h"

void PrimOp::apply(State& out, const State& s, double t) {
    out = s;
    inplace_apply(out, t);
}

void PrimOp::inplace_apply(State& s, double t) {
    Assert(the_freedom < s.n_freedoms());

    sz_t skip = s.skip(the_freedom);
    Freedom v{s.data(), skip};

    sz_t next_skip = s.skip(the_freedom - 1);

    for (sz_t j = 0; j < s.total_dims(); j += next_skip) {
        for (sz_t i = 0; i < skip; ++i) {
            v.data = s.data() + i + j;
            apply_inplace(v, t);
        }
    }
}