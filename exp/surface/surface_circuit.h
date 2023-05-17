#ifndef _SURFACE_CIRCUIT_H
#define _SURFACE_CIRCUIT_H

#include "op/prim_op.h"
#include "op/prim2_op.h"

#include "nlohmann/json.hpp"

class HGate: public PrimOp {
public:
    explicit HGate(sz_t the_freedom) : PrimOp(the_freedom) {}

    HGate& on(sz_t the_freedom) { this->the_freedom = the_freedom; return *this; }

protected:
    void apply_inplace(Freedom& v, double t = 0) override {
        Complex tmp = M_SQRT1_2 * (v(0) - v(1));
        v(0) = M_SQRT1_2 * (v(0) + v(1));
        v(1) = tmp;
    }
};

class CXGate: public Prim2Op {
public:
    CXGate(sz_t the_freedom1, sz_t the_freedom2) : Prim2Op(the_freedom1, the_freedom2) {}

    CXGate& on(sz_t the_freedom1, sz_t the_freedom2) {
        this->the_freedom1 = the_freedom1;
        this->the_freedom2 = the_freedom2;
        return *this;
    }

protected:
    void apply_inplace(Freedom2& v, double t = 0) override {
        std::swap(v(1, 0), v(1, 1));
    }
};

class Ins{
public:
    enum class Type {
        None = 0,
        H = 1,
        CX = 2,
    };

    Ins(Ins::Type type, sz_t target1, sz_t target2 = 0) : type(type), target1(target1), target2(target2) {}

    Ins::Type type;
    sz_t target1;
    sz_t target2;
};

class Circuit {
public:
    Circuit() = default;
    
    void diagnose(std::ostream& out, const std::string& indent = "");

    void run(State& s);

    static Circuit parse_circuit(nlohmann::json& def);

    // No Copy
    Circuit(const Circuit&) = delete;
    Circuit& operator=(const Circuit&) = delete;

    // Move
    Circuit(Circuit&&) = default;
    Circuit& operator=(Circuit&&) = default;

private:
    HGate h{0};
    CXGate cx{0, 0};

    std::vector<Ins> ins_list;
};

#endif // _SURFACE_CIRCUIT_H
