#ifndef _LINDBLAD_H
#define _LINDBLAD_H

#include "op/prim_op.h"

template<typename SCALAR_TYPE=double>
class ScaledSigmaX : public PrimOp {
public:
    ScaledSigmaX(sz_t the_freedom, SCALAR_TYPE scale) : PrimOp(the_freedom), scale(scale) {}

    ScaledSigmaX& on(sz_t the_freedom) { this->the_freedom = the_freedom; return *this; }

protected:
    SCALAR_TYPE scale;
    void apply_inplace(Freedom& v, double t = 0) override {
        Complex tmp = v(1);
        v(1) = scale * v(0);
        v(0) = scale * tmp;
    }
};

template<typename SCALAR_TYPE=double>
class ScaledSigmaY : public PrimOp {
public:
    ScaledSigmaY(sz_t the_freedom, SCALAR_TYPE scale) : PrimOp(the_freedom), scale(scale) {}

    ScaledSigmaY& on(sz_t the_freedom) { this->the_freedom = the_freedom; return *this; }

protected:
    SCALAR_TYPE scale;
    void apply_inplace(Freedom& v, double t = 0) override {
        Complex tmp = v(0);
        v(0) = _MI * scale * v(1);
        v(1) = _I * scale * tmp;
    }
};

template<typename SCALAR_TYPE=double>
class ScaledSigmaZ : public PrimOp {
public:
    ScaledSigmaZ(sz_t the_freedom, SCALAR_TYPE scale) : PrimOp(the_freedom), scale(scale) {}

    ScaledSigmaZ& on(sz_t the_freedom) { this->the_freedom = the_freedom; return *this; }

protected:
    SCALAR_TYPE scale;
    void apply_inplace(Freedom& v, double t = 0) override {
        v(0) *= scale;
        v(1) *= (-1. * scale);
    }
};

template<typename SCALAR_TYPE=double>
class ScaledSigmaPlus : public PrimOp {
public:
    ScaledSigmaPlus(sz_t the_freedom, SCALAR_TYPE scale) : PrimOp(the_freedom), scale(scale) {}

    ScaledSigmaPlus& on(sz_t the_freedom) { this->the_freedom = the_freedom; return *this; }

protected:
    SCALAR_TYPE scale;
    void apply_inplace(Freedom& v, double t = 0) override {
        v(0) = scale * v(1);
        v(1) = 0.;
    }
};

template<typename SCALAR_TYPE=double>
class ScaledSigmaMinus : public PrimOp {
public:
    ScaledSigmaMinus(sz_t the_freedom, SCALAR_TYPE scale) : PrimOp(the_freedom), scale(scale) {}

    ScaledSigmaMinus& on(sz_t the_freedom) { this->the_freedom = the_freedom; return *this; }
protected:
    SCALAR_TYPE scale;
    void apply_inplace(Freedom& v, double t = 0) override {
        v(1) = scale * v(0);
        v(0) = 0.;
    }
};

#endif // _LINDBLAD_H
