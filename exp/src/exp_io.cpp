#include "exp_io.h"

dm_t load_dm(const char* fname) {
    auto fin = safe_open_r(fname, std::ios::binary);

    unsigned m, n;
    fin.read((char*)&m, sizeof(m));
    fin.read((char*)&n, sizeof(n));

    dm_t mat(m, n);
    for (unsigned i = 0; i < m; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            double real, imag;
            fin.read((char*)&real, sizeof(real));
            fin.read((char*)&imag, sizeof(imag));
            mat(i,j) = Complex{real, imag};
        }
    }

    return mat;
}

std::vector<double> load_vec(const char* fname) {
    auto fin = safe_open_r(fname, std::ios::binary);

    unsigned n;
    fin.read((char*)&n, sizeof(n));

    std::vector<double> vec;
    vec.resize(n);
    fin.read((char*)vec.data(), sizeof(double) * n);

    return vec;
}

ket_t load_ket(const char* fname) {
    auto fin = safe_open_r(fname, std::ios::binary);

    unsigned n;
    fin.read((char*)&n, sizeof(n));

    ket_t ket;
    ket.resize(n);

    for (unsigned i = 0; i < n; ++i) {
        double real, imag;
        fin.read((char*)&real, sizeof(real));
        fin.read((char*)&imag, sizeof(imag));
        ket[i] = Complex{real, imag};
    }

    return ket;
}

void dump_ket(const ket_t& ket, const char* fname) {
    auto fout = safe_open_w(fname, std::ios::binary);

    unsigned n = ket.size();

    fout.write((char*)&n, sizeof(n));

    for (unsigned i = 0; i < n; ++i) {
        double real = ket[i].real();
        double imag = ket[i].imag();

        fout.write((char*)&real, sizeof(real));
        fout.write((char*)&imag, sizeof(imag));
    }
}
