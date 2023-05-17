#include "surface_sys.h"

#include "exp_io.h"
#include "exp_sop.h"
#include "exp_lindblad.h"

Sys::Sys(sz_t T_rzx, StatePool& pool, std::mt19937& eng) : T_rzx(T_rzx), unr(nullptr, L, Ldag, sum_LdagL, pool, eng), pool(pool), eng(eng) {}

Sys::~Sys() {
    for (auto& l : L) delete l;
    for (auto& ldag : Ldag) delete ldag;
}

void Sys::reset_zz_strength(const std::vector<double>& zz_strength) {
    zz_enabled = false;
    for (std::size_t i = 0; i < topo.size(); ++i) {
        auto& e = topo[i];
        auto s = zz_strength[i];
        if (s) {
            auto tmp = s * embed(n_qubits, {SSigmaZ, SSigmaZ}, {e.first, e.second});
            if (!zz_enabled) zz = tmp;
            else zz += tmp;
            zz_enabled = true;
        }
    }
}

void Sys::reset_relaxation(const std::vector<double>& T1, const std::vector<double>& T2) {
    // Lindblad
    for (auto& l : L) delete l;
    for (auto& ldag : Ldag) delete ldag;
    L.clear(); Ldag.clear();
    
    std::vector<SOp> sL;
    std::vector<SOp> sLdag;
    for (sz_t i = 0; i < n_qubits; ++i) {
        std::vector<Lindblad> Lindblads = relaxation(i, T1[i], T2[i]);
        std::vector<Lindblad> sLindblads = relaxation_s(T1[i], T2[i]);

        for (auto& ld : Lindblads) {
            L.push_back(ld.L);
            Ldag.push_back(ld.Ldag);
        }

        for (auto& sld : sLindblads) {
            SOp& l = *dynamic_cast<SOp*>(sld.L);
            SOp& ldag = *dynamic_cast<SOp*>(sld.Ldag);

            sL.push_back(embed(n_qubits, {l}, {i}));
            sLdag.push_back(embed(n_qubits, {ldag}, {i}));

            delete sld.L;
            delete sld.Ldag;
        }
    }

    if (!sL.empty()) sum_LdagL = sumLdagL(sL, sLdag);

    unr.set_lindblads(L, Ldag, sum_LdagL);
}

void Sys::diagnose(std::ostream& out, const std::string& indent) {
    out << indent << "Type: " << type() << std::endl;
    out << indent << "zz enabled: " << zz_enabled << std::endl;
    out << indent << "zz dims: [" << zz.matrix().rows() << ", " << zz.matrix().columns() << "] with nnz " << zz.matrix().nonZeros() << std::endl;
    out << indent << "#L: " << L.size() << std::endl;
    out << indent << "#Ldag: " << Ldag.size() << std::endl;
    out << indent << "sum_LdagL dims: [" << sum_LdagL.matrix().rows() << ", " << sum_LdagL.matrix().columns() << "]" << std::endl;
}

/////////////////////////////////////////////////// OptSys

OptSys::OptSys(sz_t T_rzx, double p1, double p2, double p3, StatePool& pool, std::mt19937& eng)
    : Sys(T_rzx, pool, eng),
      p_zx_0(p1 / (p1 + p2 + p3) * T_rzx),
      p_zx_1(p2 / (p1 + p2 + p3) * T_rzx),
      p_zx_2(p3 / (p1 + p2 + p3) * T_rzx) {
    // Pulse parameter
    constexpr char* r90_params_bin = "./pulse/r90.bin";
    constexpr char* r180_params_bin = "./pulse/r180.bin";
    constexpr char* id_params_bin = "./pulse/I.bin";

    // Pulse
    _p_r90.reset_params(load_vec(r90_params_bin).data());
    _p_r180.reset_params(load_vec(r180_params_bin).data());
    _p_id.reset_params(load_vec(id_params_bin).data());
}

double OptSys::apply_zx(State& out, const State& in, sz_t ctl, sz_t trg, double t) {
    if (t < p_zx_0.T) {
        szx.on(ctl, trg).apply(out, in, t);
        return p_zx_0(t);

    } else if (t < p_zx_0.T + p_zx_1.T) {
        sx.on(trg).apply(out, in, t);
        return p_zx_1(t);

    } else {
        szx.on(ctl, trg).apply(out, in, t);
        return p_zx_2(t);
    }
}

/////////////////////////////////////////////////// Opt2Sys

SysFactory::Register<Opt2Sys> _reg_opt2_sys;

Opt2Sys::Opt2Sys(sz_t T_rzx, StatePool& pool, std::mt19937& eng)
    : OptSys(T_rzx, 1., 3., 1., pool, eng) {

    constexpr char* zx_params_bin = "./pulse/zx2.bin";
    auto zx_params = load_vec(zx_params_bin);
    for (auto& p : zx_params) p *= (2. / T_rzx);
    
    sz_t n_params = 0;
    p_zx_0.reset_params(zx_params.data() + n_params); n_params += p_zx_0.n_params();
    p_zx_1.reset_params(zx_params.data() + n_params); n_params += p_zx_1.n_params();
    p_zx_2.reset_params(zx_params.data() + n_params); n_params += p_zx_2.n_params();
    Assert(n_params == zx_params.size());
}

/////////////////////////////////////////////////// Opt3Sys

SysFactory::Register<Opt3Sys> _reg_opt3_sys;

Opt3Sys::Opt3Sys(sz_t T_rzx, StatePool& pool, std::mt19937& eng)
    : OptSys(T_rzx, 1., 4., 1., pool, eng) {

    constexpr char* zx_params_bin = "./pulse/zx3.bin";
    
    double pi_2 = M_PI_2 * 3. / T_rzx;
    p_zx_0.reset_params(&pi_2);
    
    auto zx_params = load_vec(zx_params_bin);
    for (auto& p : zx_params) p *= (3. / T_rzx);
    p_zx_1.reset_params(zx_params.data());
    
    p_zx_2.reset_params(&pi_2);
}

/////////////////////////////////////////////////// GauSys

SysFactory::Register<GauSys> _reg_gau_sys;

double GauSys::apply_zx(State& out, const State& in, sz_t ctl, sz_t trg, double t) {
    szx.on(ctl, trg).apply(out, in, t);
    return _p_zx(t);
}
