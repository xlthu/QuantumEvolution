#include <iostream>
#include <chrono>
#include <bitset>

#include "qe.h"
#include "exp_timer.h"
#include "exp_sop.h"
#include "exp_types.h"
#include "exp_pulse.h"

using namespace std;
using namespace chrono;

void test1() {
    blas_init(1);

    sz_t N = 17;

    State s1{std::vector<sz_t>{2}};
    s1[0] = 2;
    s1[1] = 3;
    s1.normalize();

    State psi1{std::vector<const State*>(N, &s1)};
    State psi_tmp(psi1);
    State psi_p(psi1);
    State psi_p2(psi1);
    State psi_s(psi1);

    std::vector<SigmaX> X;
    for (sz_t i = 0; i < N; ++i) X.emplace_back(i);
    std::vector<SigmaZ> Z;
    for (sz_t i = 0; i < N; ++i) Z.emplace_back(i);

    std::vector<SigmaZX> ZX;
    for (sz_t i = 1; i < N; ++i) ZX.emplace_back(i - 1, i);

    std::vector<SOp> ZXs;
    for (sz_t i = 1; i < N; ++i) ZXs.push_back(embed(N, {SSigmaZ, SSigmaX}, {i-1, i}));

    // Z[0].apply(psi_tmp, psi1, 0);
    // X[1].apply(psi_p, psi_tmp, 0);

    // ZX[0].apply(psi_p2, psi1, 0);

    // ZXs[0].apply(psi_s, psi1, 0);

    // cout << max(abs(psi_p.vector() - psi_s.vector())) << endl;
    // cout << max(abs(psi_s.vector() - psi_p2.vector())) << endl;

    sz_t M = 100;
    {
        Timer t{"prim"};
        
        for (sz_t m = 0; m < M; ++m)
        for (sz_t i = 1; i < N; ++i) {
            Z[i - 1].apply(psi_tmp, psi1, 0);
            X[i].apply(psi_p, psi_tmp, 0);
        }
    }

    {
        Timer t{"prim2"};

        for (sz_t m = 0; m < M; ++m)
        for (auto& zx : ZX) {
            zx.apply(psi_p2, psi1, 0);
        }
    }

    for (auto& zxs : ZXs) zxs.apply(psi_s, psi1, 0);
    {
        Timer t{"sop"};

        for (sz_t m = 0; m < M; ++m)
        for (auto& zxs : ZXs) {
            zxs.apply(psi_s, psi1, 0);
        }
    }

    cout << max(abs(psi_p.vector() - psi_s.vector())) << endl;
    cout << max(abs(psi_s.vector() - psi_p2.vector())) << endl;
}

void test2() {
    const std::vector<std::pair<sz_t, sz_t>> edges{
        {10, 0}, {10, 3},
        {0, 11}, {3, 11}, {3, 13}, {6, 13}, {6, 16},
        {11, 1}, {11, 4}, {13, 4}, {13, 7}, {16, 7},
        {1, 9}, {1, 12}, {4, 12}, {4, 14}, {7, 14},
        {9, 2}, {12, 2}, {12, 5}, {14, 5}, {14, 8},
        {5, 15}, {8, 15}
    };

    sz_t N = 17;

    State s1{std::vector<sz_t>{2}};
    s1[0] = 2;
    s1[1] = 3;
    s1.normalize();

    State psi1{std::vector<const State*>(N, &s1)};
    State tmp{psi1};
    State psi_s{psi1};
    State psi_s2{psi1};
    State psi_v{psi1};

    SOp zz;
    for (size_t i = 0; i < edges.size(); ++i) {
        auto& e = edges[i];
        if (i == 0) zz = embed(N, {SSigmaZ, SSigmaZ}, {e.first, e.second});
        else zz += embed(N, {SSigmaZ, SSigmaZ}, {e.first, e.second});
    }

    cout << zz.property().diagonal << endl;

    zz.apply(tmp, psi1, 0);

    sz_t M = 100;
    tmp = 0;
    {
        Timer t{"sop1"};
        for (sz_t i = 0; i < M; ++i) {
            zz.apply(tmp, psi1, 0);
            psi_s += tmp;
        }
    }

    tmp = 0;
    {
        Timer t{"sop2"};
        for (sz_t i = 0; i < M; ++i) {
            zz.axpy_apply(psi_s2, 1.0, psi1, 0);
        }
    }

    // vmlSetMode(VML_LA | VML_FTZDAZ_OFF | VML_ERRMODE_NOERR);
    // {
    //     Timer t{"v"};
    //     for (sz_t i = 0; i < M; ++i)
    //         vzMul(zz.values.size(), zz.values.data(), psi1.data(), psi_v.data());
    // }

    // {
    //     Timer t{"manual"};
    //     for (sz_t i = 0; i < M; ++i) {
    //         auto size = zz.values.size();
    //         auto a = zz.values.data();
    //         auto b = psi1.data();
    //         auto r = psi_v.data();

    //         for (sz_t j = 0; j < size; ++j) {
    //             r[i] = a[i] * b[i];
    //         }
    //     }
    // }

    cout << max(abs(psi_s.vector() - psi_s2.vector())) << endl;
}

static double numpy_rand(std::mt19937& eng) {
    int a = eng() >> 5;
    int b = eng() >> 6;
    return (a * 67108864.0 + b) / 9007199254740992.0;
}

void test_numpy_rand() {
    std::mt19937 eng{12345};
    for (int i = 0; i < 10; ++i) {
        cout << numpy_rand(eng) << endl;
    }
}

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

void teset_BigNum(const std::vector<sz_t>& dims) {
    sz_t total_dims = 1;
    for (auto d : dims) total_dims *= d;
    
    auto output = [](const BigNum& b) {
        for (auto n : b.ns) cout << n << "  ";
        cout << endl;
    };

    cout << "Create:" << endl;
    for (sz_t i = 0; i < total_dims; ++i) {
        BigNum b{dims, i};
        output(b);
    }

    cout << endl << "Add:" << endl;
    BigNum b{dims, 0};
    for (sz_t i = 0; i < total_dims; ++i, b.add_one()) {
        output(b);
    }
}

template<typename T>
ostream& operator<<(ostream& out, const std::vector<T>& arr) {
    out << "[";

    if (!arr.empty()) {
        out << arr[0];
        for (std::size_t i = 1; i < arr.size(); ++i) out << ", " << arr[i];
    }

    out << "]";

    return out;
}

template<typename Tout, typename Tin, typename Gen>
std::vector<Tout> generate(const std::vector<Tin>& in, Gen gen) {
    std::vector<Tout> out;
    for (auto& i : in) out.push_back(gen(i));
    return out;
}

void test_measure_1() {
    std::mt19937 eng;

    State s1({2}, 0);
    s1[0] = 2;
    s1[1] = 3;
    State s2({2}, 1);

    int N = 3;

    for (int i = 0; i < N; ++i) {
        cout << "[" << i << "]Measure s1 " << s1.vector() << endl;
        State tmp = s1;
        auto meas = tmp.measure({0} ,eng);
        cout << "M: " << meas << endl;
        cout << "V: " << tmp.vector() << endl;
    }

    for (int i = 0; i < N; ++i) {
        cout << "[" << i << "]Measure s2 " << s2.vector() << endl;
        State tmp = s2;
        auto meas = tmp.measure({0} ,eng);
        cout << "M: " << meas << endl;
        cout << "V: " << tmp.vector() << endl;
    }

    State s(std::vector<const State*>{&s1, &s2});

    for (int i = 0; i < N; ++i) {
        cout << "[" << i << "]Measure s1 in " << s.vector() << endl;
        State tmp = s;
        auto meas = tmp.measure({0} ,eng);
        cout << "M: " << meas << endl;
        cout << "V: " << tmp.vector() << endl;
    }

    for (int i = 0; i < N; ++i) {
        cout << "[" << i << "]Measure s2 in " << s.vector() << endl;
        State tmp = s;
        auto meas = tmp.measure({1} ,eng);
        cout << "M: " << meas << endl;
        cout << "V: " << tmp.vector() << endl;
    }

    for (int i = 0; i < N; ++i) {
        cout << "[" << i << "]Measure s1,s2 in " << s.vector() << endl;
        State tmp = s;
        auto meas = tmp.measure({0, 1} ,eng);
        cout << "M: " << meas << endl;
        cout << "V: " << tmp.vector() << endl;
    }

    for (int i = 0; i < N; ++i) {
        cout << "[" << i << "]Measure s2,s1 in " << s.vector() << endl;
        State tmp = s;
        auto meas = tmp.measure({1, 0} ,eng);
        cout << "M: " << meas << endl;
        cout << "V: " << tmp.vector() << endl;
    }
}

void test_measure_2() {
    std::mt19937 eng;

    State s1({2});
    s1[0] = 1;
    s1[1] = sqrt(10);
    State s2({2});
    s2[0] = 2;
    s2[1] = 3;

    int N = 100000;

    std::vector<sz_t> res{0, 0};
    for (int i = 0; i < N; ++i) {
        State tmp = s1;
        auto meas = tmp.measure({0}, eng);
        res[meas[0]] += 1;

        ket_t ref(2);
        ref[meas[0]] = Complex{1,0};
        Assert(ref == tmp.vector());
    }
    cout << "Measure s1 " << s1.vector() << endl;
    cout << "Result: " << res << endl;
    cout << "Prob:" << generate<double, sz_t>(res, [&N](sz_t r) { return (double)r / N; }) << endl;
    cout << endl;

    res = {0, 0};
    for (int i = 0; i < N; ++i) {
        State tmp = s2;
        auto meas = tmp.measure({0}, eng);
        res[meas[0]] += 1;

        ket_t ref(2);
        ref[meas[0]] = Complex{1,0};
        Assert(ref == tmp.vector());
    }
    cout << "Measure s2 " << s2.vector() << endl;
    cout << "Result: " << res << endl;
    cout << "Prob:" << generate<double, sz_t>(res, [&N](sz_t r) { return (double)r / N; }) << endl;
    cout << endl;

    State s(std::vector<const State*>{&s1, &s2});

    res = {0, 0};
    for (int i = 0; i < N; ++i) {
        State tmp = s;
        auto meas = tmp.measure({0}, eng);
        res[meas[0]] += 1;

        State ref({2, 2}, 0);
        ref[0] = 0;
        ref[2 * meas[0]] = s[2 * meas[0]];
        ref[2 * meas[0] + 1] = s[2 * meas[0] + 1];
        ref.normalize();

        Assert(ref.vector() == tmp.vector());
    }
    cout << "Measure s1 in " << s.vector() << endl;
    cout << "Result: " << res << endl;
    cout << "Prob:" << generate<double, sz_t>(res, [&N](sz_t r) { return (double)r / N; }) << endl;
    cout << endl;

    res = {0, 0};
    for (int i = 0; i < N; ++i) {
        State tmp = s;
        auto meas = tmp.measure({1}, eng);
        res[meas[0]] += 1;

        State ref({2, 2}, 0);
        ref[0] = 0;
        ref[meas[0]] = s[meas[0]];
        ref[2 + meas[0]] = s[2 + meas[0]];
        ref.normalize();

        Assert(ref.vector() == tmp.vector());
    }
    cout << "Measure s2 in " << s.vector() << endl;
    cout << "Result: " << res << endl;
    cout << "Prob:" << generate<double, sz_t>(res, [&N](sz_t r) { return (double)r / N; }) << endl;
    cout << endl;

    res = {0, 0, 0, 0};
    for (int i = 0; i < N; ++i) {
        State tmp = s;
        auto meas = tmp.measure({0, 1}, eng);
        res[meas[0] * 2 + meas[1]] += 1;

        State ref({2, 2}, 0);
        ref[0] = 0;
        ref[meas[0] * 2 + meas[1]] = 1;

        Assert(ref.vector() == tmp.vector());
    }
    cout << "Measure s1, s2 in " << s.vector() << endl;
    cout << "Result: " << res << endl;
    cout << "Prob:" << generate<double, sz_t>(res, [&N](sz_t r) { return (double)r / N; }) << endl;
    cout << endl;
}

void test_measure2_1() {
    std::mt19937 eng;

    State s1({2}, 0);
    s1[0] = 2;
    s1[1] = 3;
    State s2({2}, 1);
    s2[0] = 1;
    s2[1] = 1;

    int N = 3;

    for (int i = 0; i < N; ++i) {
        cout << "[" << i << "]Measure s1 " << s1.vector() << endl;
        State tmp = s1;
        auto meas = tmp.measure2(0b1,eng);
        cout << "M: " << std::bitset<1>(meas) << endl;
        cout << "V: " << tmp.vector() << endl;

        ket_t ref(2);
        ref[meas] = Complex{1,0};
        Assert(ref == tmp.vector());
    }

    for (int i = 0; i < N; ++i) {
        cout << "[" << i << "]Measure s2 " << s2.vector() << endl;
        State tmp = s2;
        auto meas = tmp.measure2(0b1,eng);
        cout << "M: " << std::bitset<1>(meas) << endl;
        cout << "V: " << tmp.vector() << endl;

        ket_t ref(2);
        ref[meas] = Complex{1,0};
        Assert(ref == tmp.vector());
    }

    State s(std::vector<const State*>{&s1, &s2});

    for (int i = 0; i < N; ++i) {
        cout << "[" << i << "]Measure s1 in " << s.vector() << endl;
        State tmp = s;
        auto meas = tmp.measure2(0b10,eng);
        cout << "M: " << std::bitset<2>(meas) << endl;
        cout << "V: " << tmp.vector() << endl;

        State ref({2, 2}, 0);
        ref[0] = 0;
        ref[meas] = s[meas];
        ref[meas + 1] = s[meas + 1];
        ref.normalize();

        Assert(ref.vector() == tmp.vector());
    }

    for (int i = 0; i < N; ++i) {
        cout << "[" << i << "]Measure s2 in " << s.vector() << endl;
        State tmp = s;
        auto meas = tmp.measure2(0b01,eng);
        cout << "M: " << std::bitset<2>(meas) << endl;
        cout << "V: " << tmp.vector() << endl;

        State ref({2, 2}, 0);
        ref[0] = 0;
        ref[meas] = s[meas];
        ref[2 + meas] = s[2 + meas];
        ref.normalize();

        Assert(ref.vector() == tmp.vector());
    }

    for (int i = 0; i < N; ++i) {
        cout << "[" << i << "]Measure s1,s2 in " << s.vector() << endl;
        State tmp = s;
        auto meas = tmp.measure2(0b11,eng);
        cout << "M: " << std::bitset<2>(meas) << endl;
        cout << "V: " << tmp.vector() << endl;

        State ref({2, 2}, 0);
        ref[0] = 0;
        ref[meas] = 1.;

        Assert(ref.vector() == tmp.vector());
    }
}

void test_measure2_2() {
    std::mt19937 eng;

    State s1({2});
    s1[0] = 1;
    s1[1] = sqrt(10);
    State s2({2});
    s2[0] = 2;
    s2[1] = 3;

    int N = 100000;

    std::vector<sz_t> res = {0, 0};
    for (int i = 0; i < N; ++i) {
        State tmp = s1;
        auto meas = tmp.measure2(0b1,eng);
        res[meas] += 1;

        ket_t ref(2);
        ref[meas] = Complex{1,0};
        Assert(ref == tmp.vector());
    }
    cout << "Measure s1 " << s1.vector() << endl;
    cout << "Result: " << res << endl;
    cout << "Prob:" << generate<double, sz_t>(res, [&N](sz_t r) { return (double)r / N; }) << endl;
    cout << endl;

    res = {0, 0};
    for (int i = 0; i < N; ++i) {
        State tmp = s2;
        auto meas = tmp.measure2(0b1,eng);
        res[meas] += 1;

        ket_t ref(2);
        ref[meas] = Complex{1,0};
        Assert(ref == tmp.vector());
    }
    cout << "Measure s2 " << s2.vector() << endl;
    cout << "Result: " << res << endl;
    cout << "Prob:" << generate<double, sz_t>(res, [&N](sz_t r) { return (double)r / N; }) << endl;
    cout << endl;

    State s(std::vector<const State*>{&s1, &s2});

    res = {0, 0};
    for (int i = 0; i < N; ++i) {
        State tmp = s;
        auto meas = tmp.measure2(0b10,eng);
        Assert((meas & 0b1) == 0);
        meas >>= 1;
        res[meas] += 1;

        State ref({2, 2}, 0);
        ref[0] = 0;
        ref[2 * meas] = s[2 * meas];
        ref[2 * meas + 1] = s[2 * meas + 1];
        ref.normalize();

        Assert(ref.vector() == tmp.vector());
    }
    cout << "Measure s1 in " << s.vector() << endl;
    cout << "Result: " << res << endl;
    cout << "Prob:" << generate<double, sz_t>(res, [&N](sz_t r) { return (double)r / N; }) << endl;
    cout << endl;

    res = {0, 0};
    for (int i = 0; i < N; ++i) {
        State tmp = s;
        auto meas = tmp.measure2(0b01,eng);
        Assert((meas & 0b10) == 0);
        res[meas] += 1;

        State ref({2, 2}, 0);
        ref[0] = 0;
        ref[meas] = s[meas];
        ref[2 + meas] = s[2 + meas];
        ref.normalize();

        Assert(ref.vector() == tmp.vector());
    }
    cout << "Measure s2 in " << s.vector() << endl;
    cout << "Result: " << res << endl;
    cout << "Prob:" << generate<double, sz_t>(res, [&N](sz_t r) { return (double)r / N; }) << endl;
    cout << endl;

    res = {0, 0, 0, 0};
    for (int i = 0; i < N; ++i) {
        State tmp = s;
        auto meas = tmp.measure2(0b11,eng);
        res[meas] += 1;

        State ref({2, 2}, 0);
        ref[0] = 0;
        ref[meas] = 1.;

        Assert(ref.vector() == tmp.vector());
    }
    cout << "Measure s1, s2 in " << s.vector() << endl;
    cout << "Result: " << res << endl;
    cout << "Prob:" << generate<double, sz_t>(res, [&N](sz_t r) { return (double)r / N; }) << endl;
    cout << endl;
}

void test_gaussianzero() {
    auto g = GaussianZero{1, 1./4, M_PI};

    sz_t N = 50;
    for (sz_t i = 0; i < N; ++i) {
        double x = i * 1. / (N - 1);
        cout << x << " : " << g(x) << endl;
    }
}

int main() {
    test_gaussianzero();
    return 0;
}