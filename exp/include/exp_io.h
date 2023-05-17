#ifndef _EXP_IO_H
#define _EXP_IO_H

#include <vector>
#include <fstream>

#include "exp_types.h"
#include "base/assertion.h"

inline std::ifstream safe_open_r(const std::string& file_name, std::ios::openmode mode = std::ios::in) {
    std::ifstream f{file_name, mode};
    Assert_msg(f.is_open(), file_name << " not exists.");
    return f;
}

inline std::ofstream safe_open_w(const std::string& file_name, std::ios::openmode mode = std::ios::out) {
    std::ofstream f{file_name, mode};
    Assert_msg(f.is_open(), file_name << " not exists.");
    return f;
}

dm_t load_dm(const char* fname);

std::vector<double> load_vec(const char* fname);

ket_t load_ket(const char* fname);

void dump_ket(const ket_t& ket, const char* fname);

#endif // _EXP_IO_H
