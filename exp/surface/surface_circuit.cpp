#include "surface_circuit.h"

std::ostream& operator<<(std::ostream& out, const Ins& ins) {
    if (ins.type == Ins::Type::H) {
        return out << "H " << ins.target1;
    } else if (ins.type == Ins::Type::CX) {
        return out << "CX " << ins.target1 << ", " << ins.target2;
    } else return out << "None";
}

void Circuit::diagnose(std::ostream& out, const std::string& indent) {
    out << indent << "Total Gates: " << ins_list.size() << std::endl;
    for (std::size_t i = 0; i < ins_list.size(); ++i) {
        out << indent << " " << i << ": " << ins_list[i] << std::endl;
    }
}

void Circuit::run(State& s) {
    for (auto& ins : ins_list) {
        if (ins.type == Ins::Type::H) {
            h.on(ins.target1).inplace_apply(s, 0);

        } else if (ins.type == Ins::Type::CX) {
            cx.on(ins.target1, ins.target2).inplace_apply(s, 0);

        } else Error("Unknown Ins type: " << static_cast<sz_t>(ins.type));
    }
}

Circuit Circuit::parse_circuit(nlohmann::json& def) {
    Assert(def.is_array());

    Circuit cc;
    for (auto& ins_j : def) {
        std::string name = ins_j["name"].get<std::string>();
        if (name == "h") {
            cc.ins_list.push_back(
                Ins(Ins::Type::H, ins_j["qubits"][0].get<sz_t>())
            );
        }
        else if (name == "cx") {
            cc.ins_list.push_back(
                Ins(Ins::Type::CX, ins_j["qubits"][0].get<sz_t>(), ins_j["qubits"][1].get<sz_t>())
            );
        }
        else Error("Unknown Ins name: " << name);
    }

    return cc;
}
