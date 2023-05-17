#ifndef _OP_H
#define _OP_H

#include "state/state.h"

class Op {
public:
    virtual ~Op() = default;

    // out = op(s;t)
    virtual void apply(State& out, const State& s, double t) = 0;

protected:
    void not_implemented(const char* msg);
};

#endif // _OP_H