#include "state/state.h"
#include "base/assertion.h"
#include "base/blas.h"

State::State(std::vector<sz_t> dims, sz_t basis) : dims(std::move(dims)) {
    init_skips();
    Assert(basis < total_dims());

    resize_amps(total_dims());
    amps[basis] = 1;
}

State::State(const std::vector<const State*>& states) {
    dims.reserve(states.size());
    for (auto& s : states) {
        Assert(s->n_freedoms() == 1 && s->total_dims() >= 2);
        dims.push_back(s->dim(0));
    }
    init_skips();

    resize_amps(total_dims());
    for (sz_t i = 0; i < total_dims(); ++i) {
        sz_t k = i;
        Complex amp = 1.;
        for (sz_t j = 0; j < n_freedoms(); ++j) {
            sz_t l = k / skip(j);
            amp *= (*states[j])[l];
            k %= skip(j);
        }
        amps[i] = amp;
    }
}

State& State::operator*=(Complex a) {
    cblas_zscal(total_dims(), &a, data(), 1);
    return *this;
}

State& State::mul(Complex a, const State& x) {
    Complex zero{0.0, 0.0};
    cblas_zaxpby(total_dims(), &a, x.data(), 1, &zero, data(), 1);
    return *this;
}

State& State::axpy(Complex a, const State& x) {
    cblas_zaxpy(total_dims(), &a, x.data(), 1, data(), 1);
    return *this;
}

State& State::axpby(Complex a, const State& x, Complex b) {
    cblas_zaxpby(total_dims(), &a, x.data(), 1, &b, data(), 1);
    return *this;
}

Complex State::inner(const State& s) const {
    Complex res;
    cblas_zdotc_sub(total_dims(), data(), 1, s.data(), 1, &res);
    return res;
}

double State::norm() const {
    return cblas_dznrm2(total_dims(), data(), 1);
}

void State::normalize() {
    double n = norm();
    if (n > 1e-8) {
        *this *= (1. / n);
    }
}

static std::uniform_real_distribution<double> measure_rnd{0.0, 1.0};
template<typename T> T pow2(T x) { return x * x; }

class BigNum {
public:
    std::vector<sz_t> ns;
    const std::vector<sz_t>& dims;

    BigNum(const std::vector<sz_t>& dims, std::size_t num) : ns(dims.size(), 0), dims(dims) {
        reset(num);
    }

    void reset(std::size_t num) {
        for (std::size_t i = dims.size(); i-- > 0;) {
            ns[i] = num % dims[i];
            num /= dims[i];
        }
    }

    void add_one() {
        for (std::size_t i = dims.size(); i-- > 0;) {
            if (ns[i] < dims[i] - 1) {
                ++ns[i];
                break;
            } else ns[i] = 0;
        }
    }

    sz_t operator[](std::size_t i) const {
        return ns[i];
    }
};

std::vector<sz_t> State::measure(const std::vector<sz_t>& frees, std::mt19937& eng) {
    // Random a prob
    double norm2 = pow2(this->norm());
    double r = measure_rnd(eng) * norm2;

    // Measure the whole state
    sz_t meas = 0;
    for (double sum = 0; meas < total_dims(); ++meas) {
        sum += std::norm(data()[meas]);
        if (sum > r) break;
    }
    BigNum meas_big{dims, meas};

    // Only care about frees
    auto match = [&meas_big, &frees](const BigNum& num) -> bool {
        for (auto idx : frees) {
            if (num[idx] != meas_big[idx]) return false;
        }
        return true;
    };

    // Collapse
    BigNum i_big{dims, 0};
    for (sz_t i = 0; i < total_dims(); ++i, i_big.add_one()) {
        if (!match(i_big)) data()[i] = 0.;
    }
    normalize();
    
    // Return measure result
    std::vector<sz_t> ret;
    for (auto idx : frees) ret.push_back(meas_big[idx]);
    return ret;
};

sz_t State::measure2(sz_t qubits, std::mt19937& eng) {
    // Random a prob
    double norm2 = pow2(this->norm());
    double r = measure_rnd(eng) * norm2;

    // Measure the whole state
    sz_t meas = 0;
    for (double sum = 0; meas < total_dims(); ++meas) {
        sum += std::norm(data()[meas]);
        if (sum > r) break;
    }

    // Only care about frees
    meas &= qubits;
    auto match = [&meas, &qubits](sz_t idx) -> bool {
        return meas == (idx & qubits);
    };

    // Collapse
    for (sz_t i = 0; i < total_dims(); ++i) {
        if (!match(i)) data()[i] = 0.;
    }
    normalize();
    
    // Return measure result
    return meas;
}

void State::init_skips() {
    skips.resize(dims.size() + 1);

    sz_t total_dims = 1;
    for (auto i = dims.size(); i-- > 0;) {
        skips[i + 1] = total_dims;
        total_dims *= dims[i];
    }
    skips[0] = total_dims;
}

void State::resize_amps(sz_t n) {
    Assert(!shadow_mode);

    amps.resize(n);
    data_ptr = amps.data();
}
