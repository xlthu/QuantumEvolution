#ifndef _SURFACE_DECODER_H
#define _SURFACE_DECODER_H

#include <stdint.h>
#include <unordered_map>
#include <array>
#include <bitset>

#include "base/types.h"
#include "base/assertion.h"

#include "surface_correction_layer.h"
#include "nlohmann/json.hpp"

class Syndrome {
public:
    static constexpr sz_t MAX_LEN = 96;
    static sz_t SYNDROME_LEN;

    using SType = std::bitset<MAX_LEN>;
    SType syndrome;

    void shift_in(sz_t m) {
        Assert((m & SYNDROME_MASK) == m);
        syndrome <<= SYNDROME_LEN;
        syndrome |= SType(m);
    }

    bool operator==(const Syndrome& other) const {
        return syndrome == other.syndrome;
    }

    static void set_syndrome_len(sz_t syndrome_len) {
        SYNDROME_LEN = syndrome_len;
        
        SYNDROME_MASK = 0;
        for (sz_t i = 0; i < syndrome_len; ++i) {
            SYNDROME_MASK |= (1 << i);
        }
    }

    static Syndrome from_results(const std::vector<sz_t>& results) {
        Syndrome res;
        for (auto r : results) res.shift_in(r);
        return res;
    }

private:
    static sz_t SYNDROME_MASK;
};

std::ostream& operator<<(std::ostream& out, const Syndrome& syn);

inline void output_meas(sz_t meas) {
    for (sz_t mask = (1 << (Syndrome::SYNDROME_LEN - 1)); mask != 0; mask >>= 1)
        std::cout << ((meas & mask) ? '1' : '0');
}

namespace std {  
struct hash<Syndrome> {
    std::hash<typename Syndrome::SType> hash_fn;
    std::size_t operator()(const Syndrome& s) const { return hash_fn(s.syndrome); }
};
}

class Decoder {
public:
    Decoder() = default;
    ~Decoder();

    // No Copy
    Decoder(const Decoder&) = delete;
    Decoder& operator=(const Decoder&) = delete;

    // Move
    Decoder(Decoder&&) = default;
    Decoder& operator=(Decoder&&) = default;

    static Decoder parse_decoder(Sys* sys, nlohmann::json& def, sz_t n_rounds);

    void diagnose(std::ostream& out, const std::string& indent = "");

    sz_t measX;
    std::unordered_map<Syndrome, CorrectionLayer*> decodeX;

    sz_t measZ;
    std::unordered_map<Syndrome, CorrectionLayer*> decodeZ;

private:
    static void parse(Sys* sys, nlohmann::json& def, sz_t& meas, std::unordered_map<Syndrome, CorrectionLayer*>& decode, bool is_Z_correction, sz_t n_rounds);
};


#endif // _SURFACE_DECODER_H
