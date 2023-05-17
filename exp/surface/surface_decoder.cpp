#include "surface_decoder.h"

sz_t Syndrome::SYNDROME_LEN = 0;
sz_t Syndrome::SYNDROME_MASK = 0;

std::ostream& operator<<(std::ostream& out, const Syndrome& syn) {
    auto& s = syn.syndrome;

    sz_t n = ((Syndrome::MAX_LEN - 1) / Syndrome::SYNDROME_LEN) * Syndrome::SYNDROME_LEN;

    for (sz_t i = Syndrome::MAX_LEN; i-- > 0;) {
        out << s[i];
        if (i == n) {
            out << ' ';
            n -= Syndrome::SYNDROME_LEN;
        }
    }

    return out;
}

Decoder::~Decoder() {
    for (auto item : decodeX) if (item.second) delete item.second;
    for (auto item : decodeZ) if (item.second) delete item.second;
}

void Decoder::parse(Sys* sys, nlohmann::json& def, sz_t& meas, std::unordered_map<Syndrome, CorrectionLayer*>& decode, bool is_Z_correction, sz_t n_rounds) {
    meas = def["measure"].get<sz_t>();

    auto& syndrome = def["syndrome"];
    auto& correction = def["correction"];

    Assert(syndrome.is_array());
    Assert(correction.is_array());
    Assert(syndrome.size() == correction.size());

    auto syn_it = syndrome.begin();
    auto corr_it = correction.begin();

    while (syn_it != syndrome.end()) {
        // syndrome
        auto& syn_j = *syn_it;
        Assert(syn_j.is_array() && syn_j.size() == n_rounds);

        std::vector<sz_t> results;
        for (auto& s_j : syn_j) results.push_back(s_j.get<sz_t>());
        Syndrome syn = Syndrome::from_results(results);

        // correction
        auto& corr_j = *corr_it;
        Assert(corr_j.is_array() && corr_j.size() == 2);
        CorrectionLayer* ly = nullptr;

        if (is_Z_correction) {
            Assert(corr_j[0].is_array());
            Assert(corr_j[1].is_array() && corr_j[1].size() == 0);

            if (!corr_j[0].empty()) {
                ly = new CorrectionLayer();
                for (auto& q : corr_j[0]) ly->z.push_back(q.get<sz_t>());
            }

        } else {
            Assert(corr_j[0].is_array() && corr_j[1].is_array());

            if (!corr_j[0].empty()) {
                ly = new CorrectionLayer();
                for (auto& q : corr_j[0]) ly->x.push_back(q.get<sz_t>());
                for (auto& q : corr_j[1]) ly->id.push_back(q.get<sz_t>());
            }
        }
        if (ly) ly->init(sys);

        Assert_msg(decode.find(syn) == decode.end(), "Duplicated syndrome: " << syn);
        decode.insert({syn, ly});
        
        ++syn_it; ++corr_it;
    }
}

Decoder Decoder::parse_decoder(Sys* sys, nlohmann::json& def, sz_t n_rounds) {
    // X
    Decoder d;
    parse(sys, def["X"], d.measX, d.decodeX, true, n_rounds);

    // Z
    parse(sys, def["Z"], d.measZ, d.decodeZ, false, n_rounds);

    return d;
}

void Decoder::diagnose(std::ostream& out, const std::string& indent) {
    auto out_ = [&out, &indent](sz_t measure, const std::unordered_map<Syndrome, CorrectionLayer*>& decode) {
        out << indent << " measure: "; output_meas(measure); std:: cout << std::endl;
        out << indent << " #snydrome: " << decode.size() << std::endl;

        // sz_t correct = 0;
        // for (auto& d : decode) correct += (d.second == nullptr);

        // out << indent << " #correct syndrome: " << correct << std::endl;
    };

    out << indent << "X:" << std::endl;
    out_(measX, decodeX);
    out << indent << "Z:" << std::endl;
    out_(measZ, decodeZ);
}