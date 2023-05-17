#ifndef _SURFACE_CORRECTION_LAYER_H
#define _SURFACE_CORRECTION_LAYER_H

#include <iostream>
#include <vector>

#include "op/op.h"
#include "surface_sys.h"

class CorrectionLayer: public Op {
public:
    void diagnose(std::ostream& out, const std::string& indent = "");

    void apply_layer(State& s, ODESolver* solver, bool noise_free = false);

    void apply(State& out, const State& s, double t) override;

    // No Copy
    CorrectionLayer(const CorrectionLayer&) = delete;
    CorrectionLayer& operator=(const CorrectionLayer&) = delete;

    // Move
    CorrectionLayer(CorrectionLayer&&) = default;
    CorrectionLayer& operator=(CorrectionLayer&&) = default;

    std::vector<sz_t> x; // X
    std::vector<sz_t> id; // Id
    std::vector<sz_t> z; // Z

    bool is_Z_only() const { return !z.empty(); }

    Sys* sys;

private:
    void init(Sys* sys);

    friend class Decoder;

    CorrectionLayer() = default;
};

#endif // _SURFACE_CORRECTION_LAYER_H
