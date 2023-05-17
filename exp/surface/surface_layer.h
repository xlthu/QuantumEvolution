#ifndef _SURFACE_LAYER_H
#define _SURFACE_LAYER_H

#include <iostream>
#include <vector>
#include <random>

#include "op/sop.h"
#include "surface_sys.h"

#include "nlohmann/json.hpp"

class Layer: public Op {
public:
    void diagnose(std::ostream& out, const std::string& indent = "");

    void apply_layer(State& s, ODESolver* solver);

    void apply(State& out, const State& s, double t) override;

    static std::vector<Layer> parse_layers(Sys* sys, nlohmann::json& def);

    static Layer parse_meas(Sys* sys, nlohmann::json& def);

    // No Copy
    Layer(const Layer&) = delete;
    Layer& operator=(const Layer&) = delete;

    // Move
    Layer(Layer&&) = default;
    Layer& operator=(Layer&&) = default;

    SOp rz;
    std::vector<sz_t> rx; // Rx(-PI/2)
    std::vector<sz_t> ry; // Ry(PI/2)
    std::vector<sz_t> id; // Id
    std::vector<Coupling> rzx; // Rzx(PI/2)

    // Over rotation, rz is perfect
    std::vector<double> or_rx;
    std::vector<double> or_ry;
    std::vector<double> or_id;
    std::vector<double> or_rzx;

    bool is_rz = false;
    bool is_2q() const { return !rzx.empty(); }

    sz_t duration;

    Sys* sys;

private:
    void init(Sys* sys);
    
    void setup_over_roration();

    static Layer parse_rz(Sys* sys, nlohmann::json& layer);

    static Layer parse_sim(Sys* sys, nlohmann::json& layer);

    Layer() = default;
};

#endif // _SURFACE_LAYER_H
