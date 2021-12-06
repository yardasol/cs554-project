#include "matrix_arithmetic.c"
#include "poisson_matrix_creaters.c"
#include <assert.h>

/**
 * Creates an interpolation operator that coarsens every other point.
 * Requires an odd number of points on the fine grid.
 * Follows 13.3.1 in Iterative Methods for Sparse Linear Systems, Y. Saad
 */
struct A_csr* mg_1d_interpolation(unsigned int n_fine, unsigned int* n_coarse_out) {
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
 * Creates an interpolation operator that coarsens every other point in x and y dimensions.
 * This should be passed the number of fine points in *one* dimension.
 * Number of fine points should be odd.
 */
struct A_csr* mg_2d_interpolation(unsigned int n_fine, unsigned int* n_coarse_out,
                                  unsigned int* n_rows_out, unsigned int* n_cols_out) {
    assert(n_fine % 2 == 1);

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
