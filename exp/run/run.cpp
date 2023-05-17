#include "cxxopts.hpp"
#include <thread>
#include <vector>
#include <regex>

#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>

template<typename T>
static std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec) {
    if (vec.empty()) return out;
    out << vec[0];
    for (std::size_t i = 1; i < vec.size(); ++i) out << "," << vec[i];
    return out;
}

static bool term_force = false;

static void sigint_handler(int signum) {
    if (term_force) {
        std::cout << "Forced termination" << std::endl;
        exit(0);
    }
    std::cout << "Terminate all processes..." << std::endl;
    term_force = true;
}

static char escape_char = '%';
static std::string id_fmt{"%ID"};

static std::string replace_id(std::string str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        if (start_pos > 0 && str[start_pos - 1] == escape_char) {
            str.erase(start_pos, 1);
            continue;
        }
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
    return str;
}

static std::vector<std::string> process_args(const std::vector<std::string>& args, const std::string& id_str) {
    std::vector<std::string> pargs; pargs.reserve(args.size());
    for (auto& arg : args) {
        pargs.push_back(replace_id(arg, id_fmt, id_str));
    }
    return pargs;
}

static void exec(const std::string& log, std::vector<std::string>& args, unsigned int i, std::vector<unsigned int> cores) {
    std::string id_str = std::to_string(i);

    if (cores.empty()) {
        std::cout << "Internal error, abort" << std::endl;
        abort();
    }

    std::stringstream cores_ss; cores_ss << cores;

    std::vector<std::string> numactl_args{
        "numactl",
        "-C",
        cores_ss.str(),
        "--localalloc"
    };

    auto pargs = process_args(args, id_str);
    auto plog = replace_id(log, id_fmt, id_str);

    std::vector<char*> all_args;

    for (auto& arg : numactl_args) all_args.push_back(const_cast<char*>(arg.c_str()));
    for (auto& arg : pargs) all_args.push_back(const_cast<char*>(arg.c_str()));
    all_args.push_back(nullptr);

    for (auto& arg: all_args) {
        if(arg) std::cout << arg << " ";
        else std::cout << std::endl;
    }

    int fd = open(plog.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);

    if (fd < 0) {
        std::cout << "Output rediction error" << std::endl;
        return;
    }

    dup2(fd, 1);
    dup2(fd, 2);

    close(fd);

    execvp(all_args[0], all_args.data());
}

static std::regex range_re{"([0-9]+)-([0-9]+)"};

static std::vector<unsigned int> expand_range(const std::vector<std::string>& args) {
    std::vector<unsigned int> ret;

    for (auto& arg : args) {
        std::smatch m;
        if (std::regex_match(arg, m, range_re)) {
            auto begin = std::stoul(m[1].str());
            auto end = std::stoul(m[2].str()) + 1;
            for (;begin != end; ++begin) ret.push_back(begin);

        } else ret.push_back(std::stoul(arg));
    }

    return ret;
}

static void help(cxxopts::Options& options) {
    std::cout << options.help() << std::endl;
    std::cout << "Notes:" << std::endl;
    std::cout << "  %ID in -o (--log), --prog, --args will be replaced by process id." << std::endl;
}

int main(int argc, char** argv) {
    cxxopts::Options options(argv[0], "Run program on every core");
    options
      .positional_help("-- prog [args]").show_positional_help();

    options.add_options()
        ("C,cores", "Available cores (Leave empty to detect all cores)", cxxopts::value<std::vector<std::string>>()->default_value(""))
        ("c,cpp", "#Cores per process (0: Divide cores to processes equally)", cxxopts::value<unsigned int>()->default_value("0"))
        ("n,np", "#Processes to run (0: Use all available cores)", cxxopts::value<unsigned int>()->default_value("1"))
        ("o,log", "Log", cxxopts::value<std::string>()->default_value("%ID.log"))
        ("prog", "Program", cxxopts::value<std::string>())
        ("args", "Args", cxxopts::value<std::vector<std::string>>());

    options.parse_positional({"prog", "args"});

    // help
    if (argc == 1) {
        help(options);
        return 0;
    }

    // Check
    bool prog_check = false;
    for (int i = 0; i < argc; ++i) {
        if (strcmp(argv[i], "--") == 0) {
            prog_check = true;
            break;
        }
    }

    if (!prog_check) {
        std::cout << "Invalid Arguments: prog and args must be after --" << std::endl;
        return -1;
    }

    cxxopts::ParseResult result;
    try {
        result = options.parse(argc, argv);
    } catch (cxxopts::OptionException& e) {
        std::cout << "Invalid Arguments: " << e.what() << std::endl;
        return -1;
    }

    if (!result.count("prog")) {
        std::cout << "Argument Missing: prog" << std::endl;
        return -1;
    }

    // Option: cores
    auto cores_s = result["cores"].as<std::vector<std::string> >();
    auto cores = expand_range(cores_s);
    if (cores.empty()) {
        unsigned int max_ncores = std::thread::hardware_concurrency();
        for (unsigned int i = 0; i < max_ncores; ++i) cores.push_back(i);
    }
    std::cout << "Available cores (" << cores.size() << "): " << cores << std::endl;

    // Option: cpp & np
    unsigned int per_core = result["cpp"].as<unsigned int>();
    unsigned int n_procs = result["np"].as<unsigned int>();

    if (per_core == 0 && n_procs == 0) {
        std::cout << "Invalid Arguments: both c,cpp and n,np are 0" << std::endl;
        return -1;
    }

    if (per_core == 0) per_core = cores.size() / n_procs;
    if (n_procs == 0) n_procs = cores.size() / per_core;

    if (n_procs * per_core > cores.size()) {
        std::cout << "Invalid Arguments: The number of cores to use exceeds the number of available cores" << std::endl;
        return -1;
    }

    // Option: prog & args
    std::string prog = result["prog"].as<std::string>();
    std::vector<std::string> args;
    if (result.count("args")) args = result["args"].as<std::vector<std::string>>();
    args.insert(args.begin(), prog);

    // Option: log
    std::string log = result["log"].as<std::string>();

    // Run
    std::cout << "Spawn " << n_procs << " process(es), each with " << per_core << " core(s), and then wait the process(es) to finish" << std::endl;

    std::map<pid_t, unsigned int> pid2id;
    for (unsigned int i = 0; i < n_procs; ++i) {
        pid_t pid = fork();

        if (pid == -1) {
            std::cout << "[fork error] " << strerror(errno) << std::endl;
            break;
        }

        if (pid == 0) {
            // Child
            std::vector<unsigned int> my_cores(cores.begin() + i * per_core, cores.begin() + (i + 1) * per_core);

            std::cout << "Process (" << i << ", " << getpid() << ") on core " << my_cores << ": ";
            
            exec(log, args, i, my_cores);
            std::cout << "[exec (" << i << ", " << getpid() << ") error] " << strerror(errno) << std::endl;
            exit(-1);
        } else {
            // Parent
            pid2id[pid] = i;
        }
    }

    struct sigaction sigint;
    sigint.sa_handler = sigint_handler;
    sigemptyset(&sigint.sa_mask);
    sigint.sa_flags = 0;

    if (sigaction(SIGINT, &sigint, nullptr) < 0) {
        std::cout << "[sigaction error] " << strerror(errno) << '\n';
        std::cout << "  Ctrl-C cannot be catched" << '\n';
        std::cout << "  And if used, it will terminate `run` and the spawned process(es) without notification" << std::endl;
    }

    while (true) {
        int wstate;
        pid_t pid = wait(&wstate);
        if (pid < 0) {
            if (errno == ECHILD) {
                std::cout << "All process(es) have exited" << std::endl;
                break;
            } else {
                std::cout << "[wait error] " << strerror(errno) << std::endl;
            }
        } else {
            if (WIFEXITED(wstate)) {
                std::cout << "Process (" << pid2id[pid] << ", " << pid << ") has finished" << std::endl;
            } else if (WIFSIGNALED(wstate)) {
                std::cout << "Process (" << pid2id[pid] << ", " << pid << ") was terminated by " << strsignal(WTERMSIG(wstate)) << std::endl;
            } else {
                std::cout << "Unknown state for process (" << pid2id[pid] << ", " << pid << ")" << std::endl;
            }
        }
    }

    return 0;
}