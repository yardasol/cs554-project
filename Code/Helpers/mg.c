#ifndef MG_C_INCLUDED
#define MG_C_INCLUDED

#include "matrix_arithmetic.c"
#include "poisson_matrix_creaters.c"
#include <assert.h>

struct MGLevelData {
    struct A_csr* A;
    float* diag;
    struct A_csr* R; /* Represents restriction from higher level to this level */
    struct A_csr* P; /* ^ */

    unsigned int n_pts_per_dim;
    unsigned int n_pts;
    struct par_multdat pmdat;
};

struct MGData {
    struct MGLevelData* levels;
    unsigned int num_levels;
    unsigned int num_pre_relax;
    unsigned int num_post_relax;
    unsigned int num_coarsest_relax;
};

/**
 * Creates an interpolation operator that coarsens every other point.
 * Requires an odd number of points on the fine grid.
 * Follows 13.3.1 in Iterative Methods for Sparse Linear Systems, Y. Saad
 */
struct A_csr* mg_1d_interpolation_odd(unsigned int n_fine, unsigned int* n_coarse_out) {
    assert(n_fine % 2 == 1);

    unsigned int n_coarse = (n_fine-1) / 2;
    unsigned int nnz = 3 * n_coarse;

    float* data = calloc(nnz, sizeof(float));
    int* col_indices = calloc(nnz, sizeof(int));
    int* row_ptrs = calloc(n_fine + 1, sizeof(int));
    row_ptrs[n_fine] = nnz;

    struct A_csr* ret = malloc(sizeof(struct A_csr));
    ret->val = data;
    ret->col_ind = col_indices;
    ret->row_ptr = row_ptrs;

    unsigned int cur_nnz = 0;
    unsigned int cur_col = 0;
    for (unsigned int row = 0; row < n_fine; row++) {
        row_ptrs[row] = cur_nnz;

        if (row == 0 ||
            row == n_fine - 1) {
            /* coarse points at end of mesh */
            data[cur_nnz] = 0.5f;
            col_indices[cur_nnz] = cur_col;
            cur_nnz++;
        } else if (row % 2 == 1) {
            /* fine points */
            data[cur_nnz] = 1.f;
            col_indices[cur_nnz] = cur_col;
            cur_nnz++;
        } else {
            /* coarse points */
            data[cur_nnz] = 0.5f;
            col_indices[cur_nnz] = cur_col;
            cur_nnz++;

            cur_col++;
            data[cur_nnz] = 0.5f;
            col_indices[cur_nnz] = cur_col;
            cur_nnz++;
        }
    }

    if (n_coarse_out != NULL) {
        *n_coarse_out = n_coarse;
    }
    return ret;
}

/**
 * Creates an interpolation operator that coarsens every two points into one.
 * Requires an even number of points on the fine grid.
 * This is completely made up and I don't know if it's a good interpolation operator.
 */
struct A_csr* mg_1d_interpolation_even(unsigned int n_fine, unsigned int* n_coarse_out) {
    assert(n_fine % 2 == 0);

    unsigned int n_coarse = n_fine / 2;
    unsigned int nnz = n_fine;

    float* data = calloc(nnz, sizeof(float));
    int* col_indices = calloc(nnz, sizeof(int));
    int* row_ptrs = calloc(n_fine + 1, sizeof(int));
    row_ptrs[n_fine] = nnz;

    struct A_csr* ret = malloc(sizeof(struct A_csr));
    ret->val = data;
    ret->col_ind = col_indices;
    ret->row_ptr = row_ptrs;

    for (unsigned int row = 0; row < n_fine; row++) {
        row_ptrs[row] = row;
        data[row] = 0.5;
        col_indices[row] = row / 2;
    }

    if (n_coarse_out != NULL) {
        *n_coarse_out = n_coarse;
    }
    return ret;
}

/**
 * Creates a general interpolation operator in 1d that roughly halves the number of fine points.
 */
struct A_csr* mg_1d_interpolation(unsigned int n_fine, unsigned int* n_coarse_out) {
    if (n_fine % 2 == 0) {
        return mg_1d_interpolation_even(n_fine, n_coarse_out);
    } else {
        return mg_1d_interpolation_odd(n_fine, n_coarse_out);
    }
}

/**
 * Creates an interpolation operator that coarsens every other point in x and y dimensions.
 * This should be passed the number of fine points in *one* dimension.
 */
struct A_csr* mg_2d_interpolation(unsigned int n_fine, unsigned int* n_coarse_out,
                                  unsigned int* n_rows_out, unsigned int* n_cols_out) {
    unsigned int n_coarse;
    struct A_csr* P_1d = mg_1d_interpolation(n_fine, &n_coarse);
    struct A_csr* P_2d = spkron(P_1d, n_fine, n_coarse,
                                P_1d, n_fine, n_coarse,
                                n_rows_out, n_cols_out);

    free(P_1d->val);
    free(P_1d->col_ind);
    free(P_1d->row_ptr);
    free(P_1d);

    if (n_coarse_out != NULL) {
        *n_coarse_out = n_coarse;
    }

    return P_2d;
}

/**
 * Forms the Galerkin coarse grid from A_H = P^T A P = R A R^T
 * The restriction operator should be passed in because it'll be more efficient for mat-mat products.
 */
struct A_csr* mg_galerkin_product(struct A_csr* R, struct A_csr* A, unsigned int R_rows, unsigned int R_cols) {
    /* Assume R=P^T */
    /* We are computing A_H = P^T A P = R A R^T */
    /* Compute RA.  Since A = A^T, RA^T = RA */
    struct A_csr* RA = spmatSpmatTransposeProductCSR(R, R_rows, R_cols,
                                                     A, R_cols, R_cols,
                                                     NULL, NULL);

    /* Compute P^T A P = R A R^T = (RA) R^T */
    struct A_csr* RARt = spmatSpmatTransposeProductCSR(RA, R_rows, R_cols,
                                                       R,  R_rows, R_cols,
                                                       NULL, NULL);

    free(RA->val);
    free(RA->col_ind);
    free(RA->row_ptr);
    free(RA);

    return RARt;
}

#endif /* MG_C_INCLUDED */
