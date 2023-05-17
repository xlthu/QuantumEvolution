#ifndef _SOP_H
#define _SOP_H

#include <vector>
#include "op.h"

#include "blaze/math/CompressedMatrix.h"
#include <mkl_spblas.h>

enum class Symmetric {
    Not = 0,
    Normal, // include diagonal
    Anti
};

enum class Triangular {
    Not = 0, // include diagonal
    Lower,
    Upper,
};

class SOpProperty {
public:
    Symmetric symmetric = Symmetric::Not;
    bool hermitian = false;
    bool diagonal = false;
    Triangular triangular = Triangular::Not;

    SOpProperty& operator+=(const SOpProperty& p);
    SOpProperty& operator-=(const SOpProperty& p);
    SOpProperty& operator*=(const SOpProperty& p);
    SOpProperty& operator*=(Complex a);
    SOpProperty& tensor(const SOpProperty& p);

    matrix_descr mkl_mat_descr() const;

private:
    sparse_matrix_type_t mkl_mat_type() const;
    sparse_fill_mode_t mkl_fill_mode() const;
};

class SOp: public Op {
public:
    using MatType = blaze::CompressedMatrix<Complex>;
    
    ~SOp();

    SOp() = default; // Undefined behaviour if used.
    // Specify mat and property
    SOp(MatType mat, SOpProperty prop) : mat(std::move(mat)), prop(std::move(prop)) {}
    // Tensor product of ops[0,n)
    explicit SOp(const std::vector<const SOp*>& ops);

    static SOp zero_like(sz_t nrows, sz_t ncols);
    static SOp zero_like(const SOp& op);
    static SOp id_like(sz_t n);
    static SOp id_like(const SOp& op);

    // Copy
    SOp(const SOp& op) { *this = op; }
    SOp& operator=(const SOp& op) {
        mat = op.mat;
        prop = op.prop;
        is_mkl_inited = false;

        return *this;
    }

    // Move
    SOp(SOp&& op) { *this = std::move(op); }
    SOp& operator=(SOp&& op) {
        mat = std::move(op.mat);
        prop = std::move(op.prop);
        is_mkl_inited = false;

        return *this;
    }

    // Unoptimized operators
    SOp& operator+=(const SOp& op);
    friend SOp operator+(const SOp& a, const SOp& b) { return SOp(a) += b; }
    SOp& operator-=(const SOp& op);
    friend SOp operator-(const SOp& a, const SOp& b) { return SOp(a) -= b; }
    SOp& operator*=(const SOp& op);
    friend SOp operator*(const SOp& a, const SOp& b) { return SOp(a) *= b; }
    SOp& operator*=(Complex a);
    friend SOp operator*(const SOp& a, Complex b) { return SOp(a) *= b; }
    friend SOp operator*(Complex a, const SOp& b) { return SOp(b) *= a; }
    SOp& tensor(const SOp& op);
    friend SOp tensor(const SOp& a, const SOp& b) { return SOp(a).tensor(b); }

    void apply(State& out, const State& s, double t) override;
    void axpy_apply(State& out, Complex a, const State& x, double t);
    void axpby_apply(State& out, Complex a, const State& x, Complex b, double t);

    // Accessor
    const MatType& matrix() const { return mat; }
    const SOpProperty& property() const { return prop; }

private:
    MatType mat;
    SOpProperty prop;

    // For MKL
    bool is_mkl_inited = false;

    sparse_matrix_t mat_mkl = nullptr;
    matrix_descr descr_mkl;
    
    void init_mkl();
    std::vector<Complex> values;
    std::vector<MKL_INT> columns;
    std::vector<MKL_INT> rows;
};

extern const SOp SSigmaX;
extern const SOp SSigmaY;
extern const SOp SSigmaZ;
extern const SOp SSigmaPlus;
extern const SOp SSigmaMinus;
extern const SOp SId2;

#endif // _SOP_H