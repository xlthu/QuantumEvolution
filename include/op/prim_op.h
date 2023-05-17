#ifndef _PRIM_OP_H
#define _PRIM_OP_H

#include "op/op.h"

class Freedom {
public:
    Freedom(Complex* data, sz_t skip) : data(data), skip(skip) {}

    Complex& operator()(sz_t i) {
        return data[i * skip];
    }

    const Complex& operator()(sz_t i) const {
        return data[i * skip];
    }

    Complex* data;
    sz_t skip;
};

class PrimOp : public Op {
public:
    explicit PrimOp(sz_t the_freedom) : the_freedom(the_freedom) {}

    // out = op(s;t)
    void apply(State& out, const State& s, double t) override;

    // s = op(s; t)
    void inplace_apply(State& s, double t);

    sz_t the_freedom;

protected:
    // v = Op(v)
    virtual void apply_inplace(Freedom& v, double t = 0) = 0;
};

class SigmaX : public PrimOp {
public:
    explicit SigmaX(sz_t the_freedom) : PrimOp(the_freedom) {}

    SigmaX& on(sz_t the_freedom) { this->the_freedom = the_freedom; return *this; }

protected:
    void apply_inplace(Freedom& v, double t = 0) override {
        std::swap(v(0), v(1));
    }
};

class SigmaY : public PrimOp {
public:
    explicit SigmaY(sz_t the_freedom) : PrimOp(the_freedom) {}

    SigmaY& on(sz_t the_freedom) { this->the_freedom = the_freedom; return *this; }

protected:
    void apply_inplace(Freedom& v, double t = 0) override {
        Complex tmp = v(0);
        v(0) = _MI * v(1);
        v(1) = _I * tmp;
    }
};

class SigmaZ : public PrimOp {
public:
    explicit SigmaZ(sz_t the_freedom) : PrimOp(the_freedom) {}

    SigmaZ& on(sz_t the_freedom) { this->the_freedom = the_freedom; return *this; }

protected:
    void apply_inplace(Freedom& v, double t = 0) override {
        v(1) *= -1.;
    }
};

class SigmaPlus : public PrimOp {
public:
    explicit SigmaPlus(sz_t the_freedom) : PrimOp(the_freedom) {}

    SigmaPlus& on(sz_t the_freedom) { this->the_freedom = the_freedom; return *this; }

protected:
    void apply_inplace(Freedom& v, double t = 0) override {
        v(0) = v(1);
        v(1) = 0.;
    }
};

class SigmaMinus : public PrimOp {
public:
    explicit SigmaMinus(sz_t the_freedom) : PrimOp(the_freedom) {}

    SigmaMinus& on(sz_t the_freedom) { this->the_freedom = the_freedom; return *this; }

protected:
    void apply_inplace(Freedom& v, double t = 0) override {
        v(1) = v(0);
        v(0) = 0.;
    }
};

#endif // _PRIM_OP_H
