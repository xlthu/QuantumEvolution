#ifndef _SURFACE_SYS_H
#define _SURFACE_SYS_H

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <functional>

#include "state/state_pool.h"

#include "op/prim_op.h"
#include "op/prim2_op.h"
#include "op/sop.h"

#include "unraveling/jump.h"

#include "exp_pulse.h"

#include "surface_types.h"

class Sys {
public:
    // Qubits
    sz_t n_qubits;
    std::vector<sz_t> data_qubits;
    std::vector<sz_t> ancilla_qubits;
    std::vector<Coupling> topo;

    // ZZ
    SOp zz;
    bool zz_enabled = false;

    void reset_zz_strength(const std::vector<double>& zz_strength);

    // Lindblad
    std::vector<Op*> L;
    std::vector<Op*> Ldag;
    SOp sum_LdagL;

    void reset_relaxation(const std::vector<double>& T1, const std::vector<double>& T2);

    // Pulse
    virtual double p_r90(double t) = 0;
    virtual double p_r180(double t) = 0;
    virtual double p_id(double t) = 0;
    virtual double apply_zx(State& out, const State& in, sz_t ctl, sz_t trg, double t) = 0;

    static constexpr sz_t T_rz = 0;
    static constexpr sz_t T_single = 1;
    const sz_t T_rzx;

    // Evolution
    JumpOpt unr;

    // Hamiltonian terms
    SigmaX sx{0};
    SigmaY sy{0};
    SigmaZ sz{0};
    SigmaZX szx{0, 0};

    // Over Rotation
    bool is_over_rotation_enabled = false;
    std::normal_distribution<double> over_rotation_dis;

    // Misc
    StatePool& pool;
    std::mt19937& eng;
    static constexpr char* type() { return "Sys"; };

    Sys(sz_t T_rzx, StatePool& pool, std::mt19937& eng);

    virtual ~Sys();

    void diagnose(std::ostream& out, const std::string& indent = "");
};

template<typename... Args>
class _SysFactory {
public:
    using SysCreateFunc = std::function<Sys*(Args...)>;

    template<typename SysType>
    class Register {
    public:
        Register() {
            _SysFactory<Args...>::registry.insert({
                SysType::type(),
                [](Args&&... args) -> Sys* {
                    return new SysType{std::forward<Args>(args)...};
                }
            });
        };
    };

    static Sys* create(const std::string& type, Args... args) {
        auto it = registry.find(type);
        if (it == registry.end()) return nullptr;
        else return it->second(std::forward<Args>(args)...);
    }

    static std::vector<std::string> allowed_types() {
        std::vector<std::string> ret;
        for (auto& r : registry) ret.push_back(r.first);
        return ret;
    }

private:
    static std::unordered_map<std::string, SysCreateFunc> registry;
};

template<typename... Args>
std::unordered_map<std::string, typename _SysFactory<Args...>::SysCreateFunc> _SysFactory<Args...>::registry;

using SysFactory = _SysFactory<sz_t, StatePool&, std::mt19937&>;

class OptSys : public Sys {
public:
    OptSys(sz_t T_rzx, double p1, double p2, double p3, StatePool& pool, std::mt19937& eng);
    
    // Pulse
    FCos<5> _p_r90{T_single};
    double p_r90(double t) override { return _p_r90(t); };

    FCos<5> _p_r180{T_single};
    double p_r180(double t) override { return _p_r180(t); }

    FCos<5> _p_id{T_single};
    double p_id(double t) override { return _p_id(t); };

    FCos<1> p_zx_0;
    FCos<5> p_zx_1;
    FCos<1> p_zx_2;

    double apply_zx(State& out, const State& in, sz_t ctl, sz_t trg, double t) override;
};

class Opt2Sys : public OptSys {
public:
    Opt2Sys(sz_t T_rzx, StatePool& pool, std::mt19937& eng);
    static constexpr char* type() { return "Opt2"; };
};

class Opt3Sys : public OptSys {
public:
    Opt3Sys(sz_t T_rzx, StatePool& pool, std::mt19937& eng);
    static constexpr char* type() { return "Opt3"; };
};

class GauSys : public Sys {
public:
    GauSys(sz_t T_rzx, StatePool& pool, std::mt19937& eng)
        : Sys(T_rzx, pool, eng), _p_zx(T_rzx, T_rzx / 4., M_PI / 2.) {}
    static constexpr char* type() { return "Gau"; };

    // Pulse
    GaussianZero _p_r90{T_single, T_single / 4., M_PI / 2.};
    double p_r90(double t) override { return _p_r90(t); };

    GaussianZero _p_r180{T_single, T_single / 4., M_PI};
    double p_r180(double t) override { return _p_r180(t); }

    double p_id(double t) override { return 0; };

    GaussianZero _p_zx;
    double apply_zx(State& out, const State& in, sz_t ctl, sz_t trg, double t) override;
};

#endif // _SURFACE_SYS_H
