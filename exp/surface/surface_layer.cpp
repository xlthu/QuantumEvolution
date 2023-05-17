#include "surface_layer.h"

#include <string>

#include "exp_util.h"
#include "exp_sop.h"
#include "exp_print.h"

void Layer::init(Sys* sys) {
    if (is_rz) Assert_msg(rz.property().diagonal, "rz must be diagonal.");

    // Sort and unique ID
    sort(rx.begin(), rx.end(), [](sz_t a, sz_t b) { return a < b; });
    
    sort(ry.begin(), ry.end(), [](sz_t a, sz_t b) { return a < b; });

    // auto orig_size = id.size();
    sort(id.begin(), id.end(), [](sz_t a, sz_t b) { return a < b; });
    // id.erase(unique(id.begin(), id.end()), id.end());

    // if (is_2q()) {
        // Assert_msg(orig_size == sys->T_rzx / sys->T_single * id.size(), "Incomplete 2Q layer detected.");
        // Assert_msg(!rzx.empty(), "2Q layer without rzx detected.");
    // }

    sort(rzx.begin(), rzx.end(), [](const Coupling& a, const Coupling& b) { 
        sz_t am = (a.first < a.second) ? a.first : a.second;
        sz_t aM = (a.first < a.second) ? a.second : a.first;
        sz_t bm = (b.first < b.second) ? b.first : b.second;
        sz_t bM = (b.first < b.second) ? b.second : b.first;
        
        return am < bm || (am == bm && aM < bM);
    });

    // Sys
    this->sys = sys;

    // Over rotation
    or_rx.resize(rx.size(), 1.);
    or_ry.resize(ry.size(), 1.);
    or_rzx.resize(rzx.size(), 1.);
    or_id.resize(duration * id.size(), 1.);
}

void Layer::setup_over_roration() {
    if (sys->is_over_rotation_enabled) {
        for (auto& o : or_rx) o = sys->over_rotation_dis(sys->eng);
        for (auto& o : or_ry) o = sys->over_rotation_dis(sys->eng);
        for (auto& o : or_rzx) o = sys->over_rotation_dis(sys->eng);
        for (auto& o : or_id) o = sys->over_rotation_dis(sys->eng);
    }
}

void Layer::diagnose(std::ostream& out, const std::string& indent) {
    out << indent << "Type: " << (is_rz ? "RZ" : "SIM") << std::endl;
    out << indent << " rz dims: [" << rz.matrix().rows() << ", " << rz.matrix().columns() << "]" << std::endl;
    out << indent << " rx: " << rx << std::endl;
    out << indent << " ry: " << ry << std::endl;
    out << indent << " id: " << id << std::endl;
    out << indent << " rzx: " << rzx << std::endl;
}

void Layer::apply_layer(State& s, ODESolver* solver) {
    if (is_rz) {
        auto tmpg = sys->pool.allocate_similar(s);
        State& tmp = tmpg.state;
        rz.apply(tmp, s, 0);
        s = tmp;
    } else {
        setup_over_roration();
        sys->unr.set_H(this);
        sys->unr.solve(solver, s, 0, duration);
    }
}

void Layer::apply(State& out, const State& s, double t) {
    out = 0;

    auto tmpg = sys->pool.allocate_similar(s);
    State& tmp = tmpg.state;

    // 1. Pulse
    if (is_2q()) {
        // Rzx(PI/2)
        for (std::size_t i = 0; i < rzx.size(); ++i) {
            double p_zx = sys->apply_zx(tmp, s, rzx[i].first, rzx[i].second, t);
            out.axpy(p_zx * or_rzx[i], tmp);
        }
    } else {
        double p_r90 = (rx.empty() & ry.empty()) ? 0 : sys->p_r90(t);

        // Rx(-PI/2)
        for (std::size_t i = 0; i < rx.size(); ++i) {
            sys->sx.on(rx[i]).apply(tmp, s, t);
            out.axpy(-p_r90 * or_rx[i], tmp);
        }

        // Ry(PI/2)
        for (std::size_t i = 0; i < ry.size(); ++i) {
            sys->sy.on(ry[i]).apply(tmp, s, t);
            out.axpy(p_r90 * or_ry[i], tmp);
        }
    }

    // Id
    sz_t step = sz_t(t);
    double p_id = sys->p_id(t - step);

    if (step >= duration) step = duration - 1;
    for (std::size_t i = 0; i < id.size(); ++i) {
        sys->sx.on(id[i]).apply(tmp, s, t);
        out.axpy(p_id * or_id[step * id.size() + i], tmp);
    }

    // 2. ZZ
    if (sys->zz_enabled) sys->zz.axpy_apply(out, 1., s, t);
}


////////////////////////////////////////////////// Parser

using nlohmann::json;

Layer Layer::parse_rz(Sys* sys, json& layer) {
    Layer ly;
    ly.is_rz = true;

    std::vector<SOp> rzs;
    std::vector<sz_t> targets;

    for (auto& ins : layer) {
        Assert(ins["name"].get<std::string>() == "rz");
        sz_t q = ins["qubits"][0].get<sz_t>();
        double theta = ins["params"][0].get<double>();
        
        rzs.push_back(::rz(theta));
        targets.push_back(q);
    }

    ly.rz = embed(sys->n_qubits, rzs, targets);

    ly.duration = sys->T_rz;
    ly.init(sys);
    return ly;
}

Layer Layer::parse_sim(Sys* sys, json& layer) {
    Layer ly;
    ly.is_rz = false;

    for (auto& ins : layer) {
        std::string name = ins["name"].get<std::string>();
        sz_t q = ins["qubits"][0].get<sz_t>();

        if (name == "id") {
            ly.id.push_back(q);
        } else {
            double theta = ins["params"][0].get<double>();

            if (name == "rx") {
                Assert(almost_equal(theta, -M_PI_2));
                ly.rx.push_back(q);

            } else if (name == "ry") {
                Assert(almost_equal(theta, M_PI_2));
                ly.ry.push_back(q);

            } else if (name == "rzx") {
                Assert(almost_equal(theta, M_PI_2));
                sz_t q1 = ins["qubits"][1].get<sz_t>();

                ly.rzx.emplace_back(q, q1);
            } else {
                Error("Unknown instruction: " << ins);
            }
        }
    }
    
    ly.duration = (ly.is_2q() ? sys->T_rzx : sys->T_single);
    ly.init(sys);
    return ly;
}

Layer Layer::parse_meas(Sys* sys, nlohmann::json& def) {
    Layer ly;
    ly.is_rz = false;

    for (auto& q_j : def["id"]) {
        ly.id.push_back(q_j.get<sz_t>());
    }
    
    ly.duration = def["T"].get<sz_t>();
    Assert_msg(ly.duration % sys->T_single == 0, "T_meas must be a multiple of T_single(" << sys->T_single << ")");

    ly.init(sys);
    return ly;
}

std::vector<Layer> Layer::parse_layers(Sys* sys, json& def){
    Assert(def.is_array());

    std::vector<Layer> layers;

    // std::cout << "#Layers: " << def.size() << std::endl;
    for (auto& layer : def) {
        Assert(layer.is_array());

        // Type check
        bool is_rz = (layer[0]["name"].get<std::string>() == "rz");
        if (is_rz) for (auto& ins : layer) Assert(ins["name"].get<std::string>() == "rz");

        // std::cout << "[Layer " << layers.size() << "] Type: " << (is_rz ? "RZ" : "SIM") << std::endl;

        if (is_rz) layers.push_back(Layer::parse_rz(sys, layer));
        else layers.push_back(Layer::parse_sim(sys, layer));
    }

    return layers;
}