#include "op/prim2_op.h"

void Prim2Op::apply(State& out, const State& s, double t) {
    out = s;
    inplace_apply(out, t);
}

void Prim2Op::inplace_apply(State& s, double t) {
    Assert(the_freedom1 < s.n_freedoms() && the_freedom2 < s.n_freedoms());

    Freedom2 v{s.data(), s.skip(the_freedom1), s.skip(the_freedom2)};

    sz_t min_freedom = the_freedom1, max_freedom = the_freedom2;
    if (min_freedom > max_freedom) {
        min_freedom = the_freedom2; max_freedom = the_freedom1;
    }

    sz_t max_skip = s.skip(min_freedom);
    sz_t next_max_skip = s.skip(min_freedom - 1);
    sz_t min_skip = s.skip(max_freedom);
    sz_t next_min_skip = s.skip(max_freedom - 1);

    for (sz_t i = 0; i < s.total_dims(); i += next_max_skip) {
        for (sz_t j = 0; j < max_skip; j += next_min_skip) {
            for (sz_t k = 0; k < min_skip; ++k) {
                v.data = s.data() + i + j + k;
                apply_inplace(v, t);
            }
        }
    }
}