#include <string>
#include <iostream>
#include <fstream>

#include "nlohmann/json.hpp"
#include "cxxopts.hpp"

#include "qe.h"
#include "exp.h"

using std::cout;
using std::endl;

using nlohmann::json;

class IdleOp : public Op {
public:
    void apply(State& out, const State& s, double t) {
        out = 0;
    }
};

void dump_result(const std::string& result_file, json& results) {
    std::ofstream fout{result_file};
    fout << results;
}

int main(int argc, char* argv[]) {
    cxxopts::Options options(argv[0], "Single Qubit Decay");
    options.add_options()
        ("t1", "T1", cxxopts::value<double>()->default_value("0"))
        ("t2", "T2", cxxopts::value<double>()->default_value("0"))
        ("i,init", "Init state type (0, 1, +, -)", cxxopts::value<char>()->default_value("1"))
        ("j,jump", "Jump type to detect (Amp, Ph)", cxxopts::value<std::string>()->default_value("Amp"))
        ("n,ntraj", "Total number of trajectories to run", cxxopts::value<unsigned int>()->default_value("1"))
        ("o,output", "Result output file", cxxopts::value<std::string>()->default_value("single.json"))
        ("h,help", "Print usage");

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        cout << options.help() << endl;
        return 0;
    }

    double T1 = result["t1"].as<double>();
    double T2 = result["t2"].as<double>();
    char init = result["init"].as<char>();
    std::string jump = result["jump"].as<std::string>();
    unsigned int ntraj = result["ntraj"].as<unsigned int>();
    std::string result_file = result["output"].as<std::string>();
    auto seed = std::random_device()();

    json results;
    results["config"] = {
        {"T1", T1},
        {"T2", T2},
        {"init_state", std::string(1, init)},
        {"jump", jump},
        {"seed", seed}
    };
    results["cycles"] = json::array();

    // Init
    State init_state{{2}, 0};
    switch (init) {
        case '0': break;
        case '1': init_state[0] = 0; init_state[1] = 1; break;
        case '+': init_state[0] = M_SQRT1_2; init_state[1] = M_SQRT1_2; break;
        case '-': init_state[0] = M_SQRT1_2; init_state[1] = -M_SQRT1_2; break;
        default: Error("Unknown init type: " << init);
    }

    // T1, T2
    std::vector<Op*> L;
    std::vector<Op*> Ldag;
    for (auto& l : relaxation(0, T1, T2)){
        L.push_back(l.L);
        Ldag.push_back(l.Ldag);
    }

    if (jump != "Amp" && jump != "Ph") Error("Unknown jump type: " << jump);
    if (T1 == 0 && jump == "Amp") Error("Detect Amp but T1 == 0");
    if (T2 == 0 && jump == "Ph") Error("Detect Ph but T2 == 0");

    int idx = 0;
    if (T1 != 0 && T2 != 0 && jump == "Ph") idx = 1;

    // Run
    IdleOp idle;
    StatePool pool;
    std::mt19937 eng{seed};

    Jump unr{&idle, L, Ldag, pool, eng};
    ZVODESolver solver{pool};

    for (unsigned int i = 0; i < ntraj; ++i) {
        auto s_g = pool.allocate_similar(init_state);
        auto& s = s_g.state;
        s = init_state;
        unr.new_trajectory();
        
        bool detected = false;
        for (sz_t cycle = 0; !detected; ++cycle) {
            unr.solve(&solver, s, 0, 1);
            
            for (auto& info : unr.jumps()) {
                if (info.idx == idx) {
                    double time = (double)cycle + info.time;
                    results["cycles"].push_back(time);
                    detected = true;
                    break;
                }
            }
        }
    }

    dump_result(result_file, results);

    for (auto& l : L) delete l;
    for (auto& ldag : Ldag) delete ldag;

    return 0;
}