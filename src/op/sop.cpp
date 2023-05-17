#include "op/sop.h"

#include <iostream>
#include "base/assertion.h"

#ifdef CHECK_MKL_SPARSE_CALL

#define BEFORE_ANY_MKL_SPARSE_CALL \
    sparse_status_t __mkl_sparse_status;

#define CALL_MKL_SPARSE(expr) \
    __mkl_sparse_status = (expr); \
    Assert_msg(__mkl_sparse_status == SPARSE_STATUS_SUCCESS, "MKL Sparse Error Code: " << __mkl_sparse_status)

#else // !CHECK_MKL_SPARSE_CALL

#define BEFORE_ANY_MKL_SPARSE_CALL \
    static_cast<void>(0);

#define CALL_MKL_SPARSE(expr) expr

#endif // CHECK_MKL_SPARSE_CALL

// SOpProperty

SOpProperty& SOpProperty::operator+=(const SOpProperty& p) {
    // Symmetric
    if (diagonal) symmetric = p.symmetric;
    else if (!p.diagonal) {
        if (symmetric != p.symmetric) symmetric = Symmetric::Not;
    }
    // Hermitian
    hermitian = (hermitian && p.hermitian);
    // Diagonal
    diagonal = (diagonal && p.diagonal);
    // Triangular
    if (diagonal) triangular = p.triangular;
    else if (!p.diagonal) {
        if (triangular != p.triangular) triangular = Triangular::Not;
    }

    return *this;
}

SOpProperty& SOpProperty::operator-=(const SOpProperty& p) {
    return *this += p;
}

SOpProperty& SOpProperty::operator*=(const SOpProperty& p) {
    bool commute = (diagonal && p.diagonal);
    // Symmetric
    if (!commute || symmetric == Symmetric::Not || p.symmetric == Symmetric::Not) symmetric = Symmetric::Not;
    else symmetric = (symmetric == p.symmetric) ? Symmetric::Normal : Symmetric::Anti;
    // Hermitian
    hermitian = (commute && hermitian && p.hermitian);
    // Diagonal
    diagonal = (diagonal && p.diagonal);
    // Triangular
    if (diagonal) triangular = p.triangular;
    else if (!p.diagonal) {
        if (triangular != p.triangular) triangular = Triangular::Not;
    }

    return *this;
}

SOpProperty& SOpProperty::operator*=(Complex a) {
    return *this;
}

SOpProperty& SOpProperty::tensor(const SOpProperty& p) {
    // Symmetric
    if (symmetric == Symmetric::Not || p.symmetric == Symmetric::Not) symmetric = Symmetric::Not;
    else symmetric = (symmetric == p.symmetric) ? Symmetric::Normal : Symmetric::Anti;
    // Hermitian
    hermitian = (hermitian && p.hermitian);
    // Diagonal
    diagonal = (diagonal && p.diagonal);
    // Triangular
    if (diagonal) triangular = p.triangular;
    else if (!p.diagonal) {
        if (triangular != p.triangular) triangular = Triangular::Not;
    }

    return *this;
}

matrix_descr SOpProperty::mkl_mat_descr() const {
    return matrix_descr{
        .type = mkl_mat_type(),
        .mode = mkl_fill_mode(),
        .diag = SPARSE_DIAG_NON_UNIT
    };
}

sparse_matrix_type_t SOpProperty::mkl_mat_type() const {
    if (diagonal) return SPARSE_MATRIX_TYPE_DIAGONAL;
    if (symmetric == Symmetric::Normal) return SPARSE_MATRIX_TYPE_SYMMETRIC;
    if (hermitian) return SPARSE_MATRIX_TYPE_HERMITIAN;
    if (triangular != Triangular::Not) return SPARSE_MATRIX_TYPE_TRIANGULAR;
    return SPARSE_MATRIX_TYPE_GENERAL;
}

sparse_fill_mode_t SOpProperty::mkl_fill_mode() const {
    switch (triangular) {
        case Triangular::Not: return SPARSE_FILL_MODE_FULL;
        case Triangular::Lower: return SPARSE_FILL_MODE_LOWER;
        case Triangular::Upper: return SPARSE_FILL_MODE_UPPER;
    }
    Assert_msg(false, "Unknown Triangular type :" << static_cast<std::underlying_type<Triangular>::type>(triangular));
}

// SOp

SOp::~SOp() {
    BEFORE_ANY_MKL_SPARSE_CALL

    if (mat_mkl) {
        CALL_MKL_SPARSE(mkl_sparse_destroy(mat_mkl));
    }
}

SOp::SOp(const std::vector<const SOp*>& ops) : SOp(*ops[0]) {
    for (std::size_t i = 1; i < ops.size(); ++i) {
        tensor(*ops[i]);
    }
}

SOp SOp::zero_like(sz_t nrows, sz_t ncols) {
    return SOp{
        MatType(nrows, ncols),
        SOpProperty{
            .symmetric = Symmetric::Normal,
            .hermitian = true,
            .diagonal = true,
            .triangular = Triangular::Not
        }
    };
}

SOp SOp::zero_like(const SOp& op) {
    return zero_like(op.mat.rows(), op.mat.columns());
}

SOp SOp::id_like(sz_t n) {
    return SOp{
        MatType(blaze::IdentityMatrix<Complex>(n)),
        SOpProperty{
            .symmetric = Symmetric::Normal,
            .hermitian = true,
            .diagonal = true,
            .triangular = Triangular::Not
        }
    };
}

SOp SOp::id_like(const SOp& op) {
    Assert(op.mat.rows() == op.mat.columns());

    return id_like(op.mat.rows());
}

SOp& SOp::operator+=(const SOp& op) {
    is_mkl_inited = false;
    mat += op.mat;
    prop += op.prop;

    return *this;
}

SOp& SOp::operator-=(const SOp& op) {
    is_mkl_inited = false;
    mat -= op.mat;
    prop -= op.prop;
    
    return *this;
}

SOp& SOp::operator*=(const SOp& op) {
    is_mkl_inited = false;
    mat *= op.mat;
    prop *= op.prop;

    return *this;
}

SOp& SOp::operator*=(Complex a) {
    is_mkl_inited = false;
    mat *= a;
    prop *= a;
    
    return *this;
}

SOp& SOp::tensor(const SOp& op) {
    is_mkl_inited = false;
    mat = blaze::kron(mat, op.mat);
    prop.tensor(op.prop);

    return *this;
}

// Apply Op - Optimized
void SOp::apply(State& out, const State& s, double t) {
    init_mkl();

    BEFORE_ANY_MKL_SPARSE_CALL

    if (prop.diagonal) {
        std::size_t size = values.size();
        const Complex* va = values.data();
        const Complex* vb = s.data();
        Complex* vr = out.data();

        for (std::size_t i = 0; i < size; ++i) vr[i] = va[i] * vb[i];

    } else {
        CALL_MKL_SPARSE(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, Complex{1.0}, mat_mkl, descr_mkl, s.data(), Complex{0.0}, out.data()));
    }
}

void SOp::axpy_apply(State& out, Complex a, const State& x, double t) {
    init_mkl();

    BEFORE_ANY_MKL_SPARSE_CALL

    if (prop.diagonal) {
        std::size_t size = values.size();
        const Complex* va = values.data();
        const Complex* vb = x.data();
        Complex* vr = out.data();

        for (std::size_t i = 0; i < size; ++i) vr[i] += a * va[i] * vb[i];

    } else {
        CALL_MKL_SPARSE(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, a, mat_mkl, descr_mkl, x.data(), Complex{1.0}, out.data()));
    }
}

void SOp::axpby_apply(State& out, Complex a, const State& x, Complex b, double t) {
    init_mkl();

    BEFORE_ANY_MKL_SPARSE_CALL

    if (prop.diagonal) {
        std::size_t size = values.size();
        const Complex* va = values.data();
        const Complex* vb = x.data();
        Complex* vr = out.data();

        for (std::size_t i = 0; i < size; ++i) vr[i] = b * vr[i] + a * va[i] * vb[i];

    } else {
        CALL_MKL_SPARSE(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, a, mat_mkl, descr_mkl, x.data(), b, out.data()));
    }
}

void SOp::init_mkl() {
    if (is_mkl_inited) return;

    auto nnz = mat.nonZeros();

    // CSR mat
    values.clear(); columns.clear(); rows.clear();
    values.reserve(nnz); columns.reserve(nnz); rows.reserve(mat.rows() + 1);

    MKL_INT row_cursor = 0;
    for (std::size_t i = 0; i < mat.rows(); ++i) {
        rows.push_back(row_cursor);
        for (auto it = mat.begin(i); it != mat.end(i); ++it) {
            values.push_back(it->value());
            columns.push_back(it->index());
        }
        row_cursor += mat.nonZeros(i);
    }
    rows.push_back(row_cursor);

    // MKL mat
    BEFORE_ANY_MKL_SPARSE_CALL

    if (mat_mkl) {
        CALL_MKL_SPARSE(mkl_sparse_destroy(mat_mkl));
    }

    CALL_MKL_SPARSE(mkl_sparse_z_create_csr(&mat_mkl, SPARSE_INDEX_BASE_ZERO, mat.rows(), mat.columns(), rows.data(), rows.data() + 1, columns.data(), values.data()));

    // descr_mkl = prop.mkl_mat_descr();
    descr_mkl.type = SPARSE_MATRIX_TYPE_GENERAL;

    CALL_MKL_SPARSE(mkl_sparse_set_mv_hint(mat_mkl, SPARSE_OPERATION_NON_TRANSPOSE, descr_mkl, 10000));
    CALL_MKL_SPARSE(mkl_sparse_optimize(mat_mkl));

    is_mkl_inited = true;
}

const SOp SSigmaX {
    {
        {0, 1},
        {1, 0},
    },
    SOpProperty{
        .symmetric = Symmetric::Normal,
        .hermitian = true,
        .diagonal = false,
        .triangular = Triangular::Not
    }
};

const SOp SSigmaY {
    {
        {0, _MI},
        {_I, 0},
    },
    SOpProperty{
        .symmetric = Symmetric::Anti,
        .hermitian = true,
        .diagonal = false,
        .triangular = Triangular::Not
    }
};

const SOp SSigmaZ {
    {
        {1, 0},
        {0, -1},
    },
    SOpProperty{
        .symmetric = Symmetric::Normal,
        .hermitian = true,
        .diagonal = true,
        .triangular = Triangular::Not
    }
};

const SOp SSigmaPlus {
    {
        {0, 1},
        {0, 0},
    },
    SOpProperty{
        .symmetric = Symmetric::Not,
        .hermitian = false,
        .diagonal = false,
        .triangular = Triangular::Upper
    }
};

const SOp SSigmaMinus {
    {
        {0, 0},
        {1, 0},
    },
    SOpProperty{
        .symmetric = Symmetric::Not,
        .hermitian = false,
        .diagonal = false,
        .triangular = Triangular::Lower
    }
};

const SOp SId2 {
    {
        {1, 0},
        {0, 1},
    },
    SOpProperty{
        .symmetric = Symmetric::Normal,
        .hermitian = true,
        .diagonal = true,
        .triangular = Triangular::Not
    }
};
