#ifndef _PRIM2_OP_H
#define _PRIM2_OP_H

#include "op/op.h"

class Freedom2 {
public:
    Freedom2(Complex* data, sz_t skip1, sz_t skip2) : data(data), skip1(skip1), skip2(skip2) {}

    Complex& operator()(sz_t i, sz_t j) {
        return data[i * skip1 + j * skip2];
    }

    const Complex& operator()(sz_t i, sz_t j) const {
        return data[i * skip1 + j * skip2];
    }

    Complex* data;
    sz_t skip1;
    sz_t skip2;
};

class Prim2Op : public Op {
public:
    Prim2Op(sz_t the_freedom1, sz_t the_freedom2) : the_freedom1(the_freedom1), the_freedom2(the_freedom2) {}

    // out = op(s;t)
    void apply(State& out, const State& s, double t) override;

    // s = op(s; t)
    void inplace_apply(State& s, double t);

    sz_t the_freedom1;
    sz_t the_freedom2;

protected:
    // v = Op(v)
    virtual void apply_inplace(Freedom2& v, double t = 0) = 0;
};

class SigmaXX : public Prim2Op {
public:
    SigmaXX(sz_t the_freedom1, sz_t the_freedom2) : Prim2Op(the_freedom1, the_freedom2) {}

    SigmaXX& on(sz_t the_freedom1, sz_t the_freedom2) {
        this->the_freedom1 = the_freedom1;
        this->the_freedom2 = the_freedom2;
        return *this;
    }

protected:
    void apply_inplace(Freedom2& v, double t = 0) override {
        Complex tmp = v(0, 0);
        v(0, 0) = v(1, 1);
        v(1, 1) = tmp;

        tmp = v(1, 0);
        v(1, 0) = v(0, 1);
        v(0, 1) = tmp;
    }
};

class SigmaZX : public Prim2Op {
public:
    SigmaZX(sz_t the_freedom1, sz_t the_freedom2) : Prim2Op(the_freedom1, the_freedom2) {}

    SigmaZX& on(sz_t the_freedom1, sz_t the_freedom2) {
        this->the_freedom1 = the_freedom1;
        this->the_freedom2 = the_freedom2;
        return *this;
    }

protected:
    void apply_inplace(Freedom2& v, double t = 0) override {
        Complex tmp = v(0, 0);
        v(0, 0) = v(0, 1);
        v(0, 1) = tmp;

        tmp = v(1, 0);
        v(1, 0) = -1. * v(1, 1);
        v(1, 1) = -1. * tmp;
    }
};

class SigmaZZ : public Prim2Op {
public:
    SigmaZZ(sz_t the_freedom1, sz_t the_freedom2) : Prim2Op(the_freedom1, the_freedom2) {}

    SigmaZZ& on(sz_t the_freedom1, sz_t the_freedom2) {
        this->the_freedom1 = the_freedom1;
        this->the_freedom2 = the_freedom2;
        return *this;
    }

protected:
    void apply_inplace(Freedom2& v, double t = 0) override {
        v(0, 1) = -v(0, 1);
        v(1, 0) = -v(1, 0);
    }
};

#endif // _PRIM2_OP_H
