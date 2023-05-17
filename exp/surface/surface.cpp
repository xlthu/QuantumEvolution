#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include "nlohmann/json.hpp"
#include "cxxopts.hpp"

#include "qe.h"
#include "exp.h"

#include "surface_sys.h"
#include "surface_layer.h"
#include "surface_decoder.h"
#include "surface_circuit.h"

using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::ifstream;
using std::ofstream;

using nlohmann::json;

/* ------------------------ Config ------------------------ */

static std::vector<double> load_d_vec(json& config, sz_t vec_size) {
    std::vector<double> ret;
    if (config.is_number()) {
        double s = config.get<double>();
        for (sz_t i = 0; i < vec_size; ++i) ret.push_back(s);

    } else {
        Assert(config.size() == vec_size);
        for (auto it = config.begin(); it != config.end(); ++it) ret.push_back(it->get<double>());
    }

    return ret;
}

/* ------------------------ Config. ------------------------ */

/* ------------------------ Result ------------------------ */

static void output_jumps(std::ostream& out, Sys* sys) {
    for (auto& jump : sys->unr.jumps()) {
        out << (((jump.idx & 1) == 0) ? "Amp(" : "Ph(") << jump.idx / 2 << ") ";
    }
    sys->unr.clear_jumps();
};

static void dump_result(const std::string& result_file, json& results) {
    ofstream fout{result_file};
    fout << results;
}

/* ------------------------ Result. ------------------------ */

/* ------------------------ Init ------------------------ */

template<typename OP>
static void apply_prim_op_on(State& s, const std::initializer_list<sz_t>& index) {
    OP op{0};
    for (auto i : index) {
        op.on(i).inplace_apply(s, 0);
    }
}

static void apply_logical_x(State& s) {
    SigmaX x{0};
    for (auto i : {2, 4, 6}) {
        x.on(i).inplace_apply(s, 0);
    }
}

static void apply_logical_z(State& s) {
    SigmaZ z{0};
    for (auto i : {0, 4, 8}) {
        z.on(i).inplace_apply(s, 0);
    }
}

static State init_logical_0(Sys* sys) {
    std::vector<sz_t> dims(sys->n_qubits, 2);
    State s{dims, 0};
    State tmp{s};

    auto apply_i_x = [&s, &tmp](const std::initializer_list<sz_t>& index) {
        tmp = s;
        apply_prim_op_on<SigmaX>(tmp, index);
        s += tmp;
    };

    apply_i_x({0, 1, 3, 4});
    apply_i_x({1, 2});
    apply_i_x({4, 5, 7, 8});
    apply_i_x({6, 7});

    s.normalize();

    return s;
}

static State init_logical_1(Sys* sys) {
    State s = init_logical_0(sys);
    apply_logical_x(s);
    return s;
}

static State init_logical_plus(Sys* sys) {
    std::vector<sz_t> dims(sys->n_qubits, 2);
    State plus{{2}, 0};
    plus[0] = M_SQRT1_2;
    plus[1] = M_SQRT1_2;

    State zero{{2}, 0};

    std::vector<const State*> ss(sys->n_qubits, &zero);
    for (auto i : sys->data_qubits) ss[i] = &plus;

    State s{ss};
    State tmp{s};

    auto apply_i_z = [&s, &tmp](const std::initializer_list<sz_t>& index) {
        tmp = s;
        apply_prim_op_on<SigmaZ>(tmp, index);
        s += tmp;
    };

    apply_i_z({0, 3});
    apply_i_z({1, 2, 4, 5});
    apply_i_z({3, 4, 6, 7});
    apply_i_z({5, 8});

    s.normalize();

    return s;
}

static State init_logical_minus(Sys* sys) {
    State s = init_logical_plus(sys);
    apply_logical_z(s);
    return s;
}

static void check_logical(const State& s, char type) {
    State tmp{s};
    
    auto check_z = [&s, &tmp](const std::initializer_list<sz_t>& index, bool plus = false) {
        tmp = s;
        apply_prim_op_on<SigmaZ>(tmp, index);
        if (!plus) tmp -= s; else tmp += s;
        cout << blaze::max(blaze::abs(tmp.vector())) << " ";
    };

    check_z({0, 3});
    check_z({1, 2, 4, 5});
    check_z({3, 4, 6, 7});
    check_z({5, 8});

    auto check_x = [&s, &tmp](const std::initializer_list<sz_t>& index, bool plus = false) {
        tmp = s;
        apply_prim_op_on<SigmaX>(tmp, index);
        if (!plus) tmp -= s; else tmp += s;
        cout << blaze::max(blaze::abs(tmp.vector())) << " ";
    };

    check_x({0, 1, 3, 4});
    check_x({1, 2});
    check_x({4, 5, 7, 8});
    check_x({6, 7});

    cout << "| ";

    if (type == '0') check_z({0, 4, 8});
    else if (type == '1') check_z({0, 4, 8}, true);
    else if (type == '+') check_x({2, 4, 6});
    else if (type == '-') check_x({2, 4, 6}, true);
    else Error("Unknown type: " << type);

    cout << endl;
}

/* ------------------------ Init. ------------------------ */

/* ------------------------ Extraction ------------------------ */

static void measurement_result_flip(sz_t& res, sz_t n_qubits, sz_t meas, double p_meas_flip, std::mt19937& eng) {
    std::uniform_real_distribution<double> rand;
    
    for (sz_t q = 0; q < n_qubits; ++q) {
        sz_t idx = (1 << (n_qubits - 1 - q));
        if ((meas & idx) && rand(eng) < p_meas_flip) {
            cout << q << " ";
            res ^= idx;
        }
    }
}

static void reset(State& s, sz_t n_qubits, sz_t meas) {
    SigmaX x{0};
    for (sz_t q = 0; q < n_qubits; ++q) {
        if (meas & (1 << (n_qubits - 1 - q))) {
            x.on(q).inplace_apply(s, 0);
        }
    }
}

/* ------------------------ Extraction. ------------------------ */

/* ------------------------ Correct ------------------------ */

static void correct(std::unordered_map<Syndrome, CorrectionLayer*>& decode, Syndrome& syndrome, State& s, ODESolver* solver, const string& indent) {
    bool noise_free = (solver == nullptr);

    cout << indent << "Syndrome: " << syndrome << endl;
    CorrectionLayer* cly = decode.at(syndrome);
    if (cly) {
        cout << indent << "Correction: " << endl;
        cly->diagnose(cout, indent + "  ");

        cly->apply_layer(s, solver, noise_free);

        cout << indent << "  Jumps: "; output_jumps(cout, cly->sys); cout << endl;
    } else {
        cout << indent << "Correction: None" << endl;
    }
}

/* ------------------------ Correct. ------------------------ */

/* ------------------------ Pauli ------------------------ */

static void apply_2q_pauli_error(State& s, const std::vector<Coupling>& cps, double p, std::mt19937& eng) {
    std::uniform_real_distribution<double> rand;
    std::uniform_int_distribution<int> rand_idx{0, 15};

    auto apply = [](int pauli_idx, State& s, sz_t qubit_idx) -> char {
        switch (pauli_idx) {
            case 0: return 'I';
            case 1: SigmaX{qubit_idx}.inplace_apply(s, 0); return 'X';
            case 2: SigmaY{qubit_idx}.inplace_apply(s, 0); return 'Y';
            case 3: SigmaZ{qubit_idx}.inplace_apply(s, 0); return 'Z';
            default: Error("Invalid pauli operator: " << pauli_idx); return '\0';
        }
    };

    for (auto& cp : cps) {
        if (rand(eng) < p) {
            int idx = 0;
            while (idx == 0) idx = rand_idx(eng);

            cout << apply(idx / 4, s, cp.first) << \
                    apply(idx % 4, s, cp.second) << \
                    "(" << cp.first << "," << cp.second << ") ";
        }
    }
}

/* ------------------------ Pauli. ------------------------ */


int main(int argc, char** argv) {
    cout.precision(10);

    std::string config_file_name;
    std::string result_file;
    {
        cxxopts::Options options(argv[0], "Surface Code Simulation");

        options.add_options()
            ("c,config", "Config file", cxxopts::value<std::string>()->default_value("config.json"), "filename")
            ("o,output", "Result output file", cxxopts::value<std::string>()->default_value("surface.json"))
            ("h,help", "Print usage");

        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            cout << options.help() << endl;
            return 0;
        }

        config_file_name = result["config"].as<std::string>();
        result_file = result["output"].as<std::string>();
    }

    // Load config
    json config;
    safe_open_r(config_file_name) >> config;
    cout << config << endl;

    // Init Math Libraries
    {
        // int n_threads = config["num_threads"].get<int>();
        // blas_init(n_threads);
        // blaze::setNumThreads(n_threads);
        cout << "#Threads:" << endl;
        cout << "  Blaze: " << blaze::getNumThreads() << endl;
        cout << "  MKL: " << mkl_get_max_threads() << endl;
    }

    // Random
    std::mt19937 eng;
    {
        auto seed = config["seed"].get<std::mt19937::result_type>();
        if (seed == 0) {
            seed = std::random_device()();
            std::cout << "Using seed = " << seed << std::endl;
        }
        eng.seed(seed);
    }

    // Sys
    StatePool pool;
    Sys* sys;
    {
        std::string sys_type = config["sys_type"].get<std::string>();
        sz_t T_rzx = config["T_zx"].get<sz_t>();
        sys = SysFactory::create(sys_type, T_rzx, pool, eng);
        if (!sys)
            Error("Unknown system type: " << sys_type << " (Allowed: " << SysFactory::allowed_types() << ")");
    }

    // qubits and topo
    sys->n_qubits = config["n_qubits"].get<sz_t>();
    for (auto& q_j : config["data_qubits"]) {
        sys->data_qubits.push_back({q_j.get<sz_t>()});
    }
    for (auto& q_j : config["ancilla_qubits"]) {
        sys->ancilla_qubits.push_back({q_j.get<sz_t>()});
    }
    for (auto& e_j : config["topo"]) {
        sys->topo.push_back({e_j[0].get<sz_t>(), e_j[1].get<sz_t>()});
    }

    // ZZ
    sys->reset_zz_strength(load_d_vec(config["ZZ"], sys->topo.size()));

    // Relaxation
    sys->reset_relaxation(load_d_vec(config["T1"], sys->n_qubits), load_d_vec(config["T2"], sys->n_qubits));

    // Over rotation
    {
        double or_mean = config["over_rotation"]["mean"].get<double>();
        double or_std = config["over_rotation"]["std"].get<double>();
        if (or_std) {
            cout << "Use over rotation ~ N(" << or_mean << ", " << or_std << ")" << endl;
            sys->is_over_rotation_enabled = true;
            sys->over_rotation_dis.param(
                std::normal_distribution<double>::param_type{or_mean + 1, or_std}
            );
        }
    }

    // Diagnose system
    cout << "System:" << endl;
    sys->diagnose(cout, "    ");

    // Load extraction
    std::vector<Layer> extraction;
    {
        json def;
        safe_open_r(config["extraction"].get<string>()) >> def;
        extraction = Layer::parse_layers(sys, def);
    }

    // Diagnose extraction
    for (size_t i = 0; i < extraction.size(); ++i) {
        auto& ly = extraction[i];
        cout << "Extraction " << i << ":" << endl;
        ly.diagnose(cout, "    ");
    }

    // Load extraction perfect
    Circuit extraction_perfect;
    {
        json def;
        safe_open_r(config["extraction_perfect"].get<string>()) >> def;
        extraction_perfect = Circuit::parse_circuit(def);
    }
    cout << "Extraction(Perfect):" << endl;
    extraction_perfect.diagnose(cout, "    ");

    // Load measurement
    Layer meas_layer = Layer::parse_meas(sys, config["measurement"]);
    if (meas_layer.duration == 0) {
        cout << "Use instant measurement" << endl;
    } else {
        cout << "Measure(T = " << meas_layer.duration << "):" << endl;
        meas_layer.diagnose(cout, "    ");
    }

    double p_meas_flip = config["measurement"]["prob_flip"].get<double>();
    cout << "P(measurement results flip) = " << p_meas_flip << endl;

    // Syndrome
    sz_t n_rounds = config["n_rounds"].get<sz_t>();
    Assert(n_rounds * sys->n_qubits <= Syndrome::MAX_LEN);
    Syndrome::set_syndrome_len(sys->n_qubits);

    // Load decoder
    Decoder decoder;
    {
        json def;
        safe_open_r(config["decoder"].get<string>()) >> def;
        decoder = Decoder::parse_decoder(sys, def, n_rounds);
    }
    cout << "Decoder:" << endl;
    decoder.diagnose(cout, "    ");

    // Load decoder perfect
    Decoder decoder_perfect;
    {
        json def;
        safe_open_r(config["decoder_perfect"].get<string>()) >> def;
        decoder_perfect = Decoder::parse_decoder(sys, def, 1);
    }
    cout << "Decoder(Perfect):" << endl;
    decoder_perfect.diagnose(cout, "    ");

    // Init state
    std::string _init_type = config["init_state"].get<std::string>();
    Assert(_init_type.size() == 1);
    char& init_type = _init_type[0];

    State init_state;
    cout << "Init as Logical " << init_type << ": ";
    switch (init_type) {
        case '0': init_state = init_logical_0(sys); break;
        case '1': init_state = init_logical_1(sys); break;
        case '+': init_state = init_logical_plus(sys); break;
        case '-': init_state = init_logical_minus(sys); break;
        default: Error("Unknown init_state type: " << init_type);
    }
    check_logical(init_state, init_type);

    State init_state_err{init_state};
    switch (init_type) {
        case '0': 
        case '1': apply_logical_x(init_state_err); break;
        case '+': 
        case '-': apply_logical_z(init_state_err); break;
        default: break;
    }

    // Pauli 2q error
    // double p_2q = config["p_2q"].get<double>();
    // cout << "Use 2Q Pauli error with p = " << p_2q << endl;

    // Solver
    ZVODESolver solver{pool};
    solver.set_atol(config["atol"].get<double>());
    solver.set_rtol(config["rtol"].get<double>());

    // Result output
    cout << "Output results to " << result_file << endl;
    json results;
    results["config"] = config;
    results["cycles"] = json::array();
    /* ------------------------ Configuration finished. ------------------------ */

    cout << "================================" << endl;
    
    while (true) {
        solver.set_suggested_first_step_size(0.0);
        
        // 1. Init
        auto s_g = pool.allocate_similar(init_state);
        auto& s = s_g.state;
        s = init_state;

        // 2. Repeat extraction x n_rounds + correct + detect logical error
        //    Until one logical error is detected
        for (sz_t cycle = 0; true; ++cycle) {
            cout << "Cycle " << cycle << ":" << endl;

            Syndrome syndromeX;
            Syndrome syndromeZ;

            // 2.1. Extraction x n_rounds
            for (sz_t i = 1; i <= n_rounds; ++i) {
                Timer t{"  Syndrome Extraction " + std::to_string(i)};

                cout << t.msg << ":" << endl;
                
                // 2.1.1 Extraction
                sys->unr.new_trajectory();
                cout << "    Layer(" << extraction.size() << "): "; cout.flush();
                for (sz_t layer_i = 0; layer_i < extraction.size(); ++layer_i) {
                    auto& layer = extraction[layer_i];
                    layer.apply_layer(s, &solver);
                    
                    cout << layer_i << ' ';

                    output_jumps(cout, sys);
                    // if (layer.is_2q() && p_2q) apply_2q_pauli_error(s, layer.rzx, p_2q, eng);

                    cout.flush();
                }
                cout << endl;

                // 2.1.2 Measure
                cout << "    Measure(T = " << meas_layer.duration << "): ";
                // 2.1.2.1 Idle + ID
                if (meas_layer.duration) {
                    meas_layer.apply_layer(s, &solver);
                    
                    output_jumps(cout, sys);
                }
                cout << std::endl;

                // 2.1.2.2 Perfect measure
                sz_t resX = s.measure2(decoder.measX, eng); syndromeX.shift_in(resX);
                sz_t resZ = s.measure2(decoder.measZ, eng); syndromeZ.shift_in(resZ);

                cout << "      X: "; output_meas(resX); cout << endl;
                cout << "      Z: "; output_meas(resZ); cout << endl;

                // 2.1.2.3 Flip
                if (p_meas_flip) {
                    cout << "    Flip:" << endl;

                    cout << "      X: "; measurement_result_flip(resX, sys->n_qubits, decoder.measX, p_meas_flip, eng); cout << endl;

                    cout << "         "; output_meas(resX); cout << endl;

                    cout << "      Z: "; measurement_result_flip(resZ, sys->n_qubits, decoder.measZ, p_meas_flip, eng); cout << endl;
                    cout << "         "; output_meas(resZ); cout << endl;
                }

                // 2.1.2.4 Reset
                reset(s, sys->n_qubits, resX | resZ);
            }

            // 2.2 Correct
            cout << "  Correct X:" << endl;
            correct(decoder.decodeX, syndromeX, s, &solver, "    ");

            cout << "  Correct Z:" << endl;
            correct(decoder.decodeZ, syndromeZ, s, &solver, "    ");

            // 2.3 Detect logical error
            cout << "  Detect Logical Error:" << endl;

            // 2.3.1 Clone
            auto s2_g = pool.allocate_similar(s);
            auto& s2 = s2_g.state;
            s2 = s;

            // 2.3.2 Apply circuit, measure, reset, correct
            extraction_perfect.run(s2);

            Syndrome syndromeX_perfect;
            Syndrome syndromeZ_perfect;

            sz_t resX = s2.measure2(decoder_perfect.measX, eng); syndromeX_perfect.shift_in(resX);
            sz_t resZ = s2.measure2(decoder_perfect.measZ, eng); syndromeZ_perfect.shift_in(resZ);

            reset(s2, sys->n_qubits, resX | resZ);

            cout << "    Correct X:" << endl;
            correct(decoder_perfect.decodeX, syndromeX_perfect, s2, nullptr, "      ");

            cout << "    Correct Z:" << endl;
            correct(decoder_perfect.decodeZ, syndromeZ_perfect, s2, nullptr, "      ");

            // 2.3.3 Check
            double fidelity_none = state_fidelity(s2, init_state);
            double fidelity_err = state_fidelity(s2, init_state_err);

            cout << "    Check:" << endl;
            cout << "      Fidelity (none vs error): " << fidelity_none << " " << fidelity_err << endl;

            if (fidelity_none > fidelity_err) {
                cout << "      Logical Error: None" << endl;
            } else {
                results["cycles"].push_back(cycle);

                cout << "      Logical Error Detected at Cycle " << cycle << endl;
                dump_result(result_file, results);
                break;
            }

        } // for (sz_t cycle = 0; true; ++cycle)
    } // while (true)

    delete sys;
    return 0;
}