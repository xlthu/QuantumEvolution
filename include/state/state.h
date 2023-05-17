#ifndef _STATE_H
#define _STATE_H

#include <vector>
#include <cstring>

#include "base/types.h"
#include "base/assertion.h"

#include "blaze/math/DynamicVector.h"

class State {
public:
    // using VecType = std::vector<Complex>;
    using VecType = blaze::DynamicVector<Complex,blaze::columnVector>;

    State() = default; // Undefined behaviour if used.

    // Basis state
    explicit State(std::vector<sz_t> dims, sz_t basis = 0);
    // Product state
    explicit State(const std::vector<const State*>& states);

    // Copy
    State(const State& s) : State(s.dims) { *this = s; }
    State& operator=(const State& s) {
        Assert(dims == s.dims);
        std::memcpy(data(), s.data(), total_dims() * sizeof(Complex));
        return *this;
    }

    // Move
    State(State&& s) = default;
    State& operator=(State&& s) = default;

    // Shadow
    State(Complex* data_ptr, std::vector<sz_t> dims)
        : data_ptr(data_ptr), shadow_mode(true), dims(std::move(dims)) { init_skips(); }
    static State shadow(Complex* data_ptr, const State& s) { return State(data_ptr, s.dims, s.skips); }

    // Similar
    void make_similar(const State& s) {
        dims = s.dims;
        skips = s.skips;
        resize_amps(s.total_dims());
    }

    // Operator
    State& operator+=(const State& s) { return axpy(1.0, s); }
    State operator+(const State& s) const { return State(*this) += s; }
    State& operator-=(const State& s) { return axpy(-1.0, s); }
    State operator-(const State& s) const { return State(*this) -= s; }
    State& operator*=(Complex a);
    State operator*(Complex a) const { return State(*this) *= a; }
    friend State operator*(Complex a, const State& s) { return State(s) * a; }

    State& mul(Complex a, const State& x); // this = a*x;
    State& axpy(Complex a, const State& x); // this += a * x;
    State& axpby(Complex a, const State& x, Complex b); // this = a * x + b * this

    State& operator=(Complex a) { std::fill(begin(), end(), a); return *this; } // fill

    // Operation
    Complex inner(const State& s) const; // <this|s>
    double norm() const;
    void normalize();

    // Measure
    std::vector<sz_t> measure(const std::vector<sz_t>& frees, std::mt19937& eng);
    sz_t measure2(sz_t qubits, std::mt19937& eng); // qubits: [0, ..., n-1]

    // Accessor
    Complex* data() { return data_ptr; }
    const Complex* data() const { return data_ptr; }
    Complex& operator[](sz_t i) { return data()[i]; }
    const Complex& operator[](sz_t i) const { return data()[i]; }
    Complex* begin() { return data(); }
    Complex* end() { return data() + total_dims(); }

    const VecType& vector() const { return amps; }

    // Property
    sz_t total_dims() const { return skips[0]; } // use after init_skips()
    sz_t n_freedoms() const { return dims.size(); }
    sz_t dim(sz_t i) const { return dims[i]; }
    sz_t skip(sz_t i) const { return skips[i + 1]; }

private:
    Complex* data_ptr = nullptr;
    bool shadow_mode = false;

    VecType amps;
    std::vector<sz_t> dims;
    std::vector<sz_t> skips;

    // Shadow
    State(Complex* data_ptr, std::vector<sz_t> dims, std::vector<sz_t> skips)
        : data_ptr(data_ptr), shadow_mode(true), dims(std::move(dims)), skips(std::move(skips)) {}
    
    // Init
    void init_skips();
    void resize_amps(sz_t n);
};

#endif // _STATE_H