#ifndef MATRIX_ARITHMETIC_C
#define MATRIX_ARITHMETIC_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include "spmatrix.h" // to include A_csr and A_dia structs
#include "helpers.c"

/*Multipliers for three families of matrices*/
float* matVecProduct1(int n, float* y, float** x)
{
  float* b;
  b = malloc(sizeof(float*) * n);
  for (int i=0; i<n; i++)
  {
    for (int j=0; j<n; j++)
    {
      b[j] += x[j][i] * y[i];
    }
  }
  return b;
}

float* matVecProductCSR1(int n, float* y, struct A_csr A)
{
  float* b; int col;
  b = malloc(sizeof(float*) * n);
  for (int i=0; i<n; i++)
  {
    b[i]=0;
    for (int j=A.row_ptr[i]; j<A.row_ptr[i+1]; j++)
    {
        col = A.col_ind[j];
        b[i] += A.val[j] * y[col];
    }
  }
  return b;
}


float* matVecProductDIA1(int n, int n_diag, float* y, struct A_dia A)
{
  float* b;
  int offs, start, end;
  b = malloc(sizeof(float*) * n);
  for (int i=0; i<n; i++)
  {
    b[i]=0;
    for (int j=0; j<n_diag; j++)
    {
        start = 0; end = n;
        offs = A.off[j];
        if (offs>0){end -= offs;}
        else{start -= offs;}
        if((i>=start)&&(i<=end)){
            b[i] += A.val[j][i] * y[i+offs];
        }
    }
  }
  return b;
}

int cooEntryComparator(const void* a, const void* b) {
    const struct A_coo_entry* A_entry = a;
    const struct A_coo_entry* B_entry = b;

    if (A_entry->row == B_entry->row) {
        if (A_entry->col < B_entry->col) {
            return -1;
        } else if (A_entry->col > B_entry->col) {
            return 1;
        } else {
            return 0;
        }
    } else {
        if (A_entry->row > B_entry->row) {
            return 1;
        } else if (A_entry->row < B_entry->row) {
            return -1;
        }
    }
}

struct A_csr* coo_to_csr(struct A_coo* A) {
    /* Make a copy of the COO entries and sort by row, column ordering.
       This gives us entries that are in the same order as CSR */
    struct A_coo_entry* data_copy = calloc(A->nnz, sizeof(struct A_coo_entry));
    memcpy(data_copy, A->data, A->nnz * sizeof(struct A_coo_entry));
    qsort(data_copy, A->nnz, sizeof(struct A_coo_entry), cooEntryComparator);

    float* val = calloc(A->nnz, sizeof(float));
    int* col_ind = calloc(A->nnz, sizeof(int));
    int* row_ptr = calloc(A->rows + 1, sizeof(int));
    row_ptr[A->rows] = A->nnz;

    struct A_csr* Acsr = calloc(1, sizeof(struct A_csr));
    Acsr->val = val;
    Acsr->col_ind = col_ind;
    Acsr->row_ptr = row_ptr;

    /* Now, assemble the CSR matrix */
    int current_row = -1;
    for (unsigned int i = 0; i < A->nnz; i++) {
        struct A_coo_entry* e = data_copy + i;
        if (current_row != e->row) {
            current_row = e->row;
            row_ptr[current_row] = i;
        }

        val[i] = e->val;
        col_ind[i] = e->col;
    }

    free(data_copy);
    return Acsr;
}

struct A_csr* spmatTransposeCSR(struct A_csr* A, unsigned int A_num_rows, unsigned int A_num_cols) {
    struct A_coo* At_coo = calloc(1, sizeof(struct A_coo));
    const unsigned int nnz = A->row_ptr[A_num_rows];

    At_coo->data = calloc(nnz, sizeof(struct A_coo_entry));
    At_coo->nnz = A->row_ptr[A_num_rows];
    At_coo->rows = A_num_cols;
    At_coo->cols = A_num_rows;

    /* Convert first to a transposed COO matrix because it is easier to work with */
    unsigned int cur_nnz = 0;
    for (unsigned int row = 0; row < A_num_rows; row++) {
        for (unsigned int ptr = A->row_ptr[row]; ptr < A->row_ptr[row+1]; ptr++) {
            unsigned int col = A->col_ind[ptr];
            float val = A->val[ptr];

            At_coo->data[cur_nnz].val = val;
            At_coo->data[cur_nnz].row = col;
            At_coo->data[cur_nnz].col = row;

            cur_nnz++;
        }
    }

    /* Now, convert from COO to CSR */
    struct A_csr* At_csr = coo_to_csr(At_coo);
    free(At_coo);

    return At_csr;
}

/**
 * Computes a sparse matrix-matrix product between two CSR matrices like
 * C = AB^T
 */
struct A_csr* spmatSpmatTransposeProductCSR(
    struct A_csr* A, unsigned int A_num_rows, unsigned int A_num_cols,
    struct A_csr* B, unsigned int B_num_rows, unsigned int B_num_cols,
    unsigned int* C_rows_out, unsigned int* C_cols_out) {

    assert(A_num_cols == B_num_cols);

    struct dynamic_array coo_entries;
    dynamic_array_initialize(&coo_entries, sizeof(struct A_coo_entry));

    struct A_coo_entry entry_temp;
    for (unsigned int i = 0; i < A_num_rows; i++) {
        for (unsigned int j = 0; j < B_num_rows; j++) {
            float sum = 0.f;

            unsigned int A_start = A->row_ptr[i];
            unsigned int A_end = A->row_ptr[i+1];
            unsigned int B_start = B->row_ptr[j];
            unsigned int B_end = B->row_ptr[j+1];

            unsigned int A_current = A_start;
            unsigned int B_current = B_start;

            /* Perform the inner product */
            while (1) {
                if (A_current >= A_end ||
                    B_current >= B_end) {
                    break;
                }

                unsigned int A_col = A->col_ind[A_current];
                unsigned int B_col = B->col_ind[B_current];

                /* Add to 'sum' when the column indices match */
                if (A_col < B_col) {
                    A_current++;
                } else if (A_col > B_col) {
                    B_current++;
                } else {
                    sum += A->val[A_current] * B->val[B_current];
                    A_current++;
                    B_current++;
                }
            }

            if (sum != 0.f) {
                entry_temp.row = i;
                entry_temp.col = j;
                entry_temp.val = sum;
                dynamic_array_push_back(&coo_entries, &entry_temp);
            }
        }
    }

    struct A_coo coo_repr;
    coo_repr.data = coo_entries.data;
    coo_repr.nnz = coo_entries.elements;
    coo_repr.rows = A_num_rows;
    coo_repr.cols = B_num_rows;

    struct A_csr* csr_repr = coo_to_csr(&coo_repr);
    dynamic_array_free(&coo_entries);

    if (C_rows_out != NULL) {
        *C_rows_out = A_num_rows;
    }
    if (C_cols_out != NULL) {
        *C_cols_out = B_num_rows;
    }

    return csr_repr;
}

/**
 * Converts a sparse CSR matrix to a dense matrix.
 */
float** csr_to_dense(struct A_csr* A, unsigned int rows, unsigned int cols) {
    const unsigned int nnz = A->row_ptr[rows];
    float** A_dense = calloc(rows, sizeof(float*));

    for (unsigned int row = 0; row < rows; row++) {
        A_dense[row] = calloc(cols, sizeof(float));

        /* Populate entries of the row */
        for (unsigned int ptr = A->row_ptr[row]; ptr < A->row_ptr[row+1]; ptr++) {
            const unsigned int col = A->col_ind[ptr];
            A_dense[row][col] = A->val[ptr];
        }
    }

    return A_dense;
}

/**
 * Computes a sparse Kronecker product between two CSR matrices.
 */
struct A_csr* spkron(
    struct A_csr* A, unsigned int A_num_rows, unsigned int A_num_cols,
    struct A_csr* B, unsigned int B_num_rows, unsigned int B_num_cols,
    unsigned int* C_rows_out, unsigned int* C_cols_out) {

    const unsigned int A_nnz = A->row_ptr[A_num_rows];
    const unsigned int B_nnz = B->row_ptr[B_num_rows];
    const unsigned int C_nnz = A_nnz * B_nnz;
    const unsigned int C_num_rows = A_num_rows * B_num_rows;
    const unsigned int C_num_cols = A_num_cols * B_num_cols;

    float* C_vals = malloc(sizeof(float) * C_nnz);
    int* C_cols = malloc(sizeof(int) * C_nnz);
    int* C_rowptrs = malloc(sizeof(int) * (C_num_rows + 1));

    /* We iterate over one row at a time in C = (A kron B).
       This allows us to create C directly in CSR format */
    unsigned int cur_nnz = 0;
    for (unsigned int row = 0; row < C_num_rows; row++) {
        C_rowptrs[row] = cur_nnz;

        /* Now, loop through the respective rows in A and B */
        unsigned int A_row = row / A_num_rows;
        unsigned int B_row = row % B_num_rows;
        for (unsigned int A_ptr = A->row_ptr[A_row]; A_ptr < A->row_ptr[A_row+1]; A_ptr++) {
            unsigned int A_col = A->col_ind[A_ptr];

            for (unsigned int B_ptr = B->row_ptr[B_row]; B_ptr < B->row_ptr[B_row+1]; B_ptr++) {
                unsigned int B_col = B->col_ind[B_ptr];
                unsigned int C_col = B_num_cols * A_col + B_col;

                C_cols[cur_nnz] = C_col;
                C_vals[cur_nnz] = A->val[A_ptr] * B->val[B_ptr];
                cur_nnz++;
            }
        }
    }
    C_rowptrs[C_num_rows] = cur_nnz;

    struct A_csr* C = malloc(sizeof(struct A_csr*));
    C->val = C_vals;
    C->col_ind = C_cols;
    C->row_ptr = C_rowptrs;

    if (C_rows_out != NULL) {
        *C_rows_out = C_num_rows;
    }
    if (C_cols_out != NULL) {
        *C_cols_out = C_num_cols;
    }

    return C;
}

/*Inner product between two matrices*/
float innerProd1(int n, float* a, float* b)
{
  float c = 0.f;
  for (int i=0; i<n; i++)
    c += a[i] * b[i];
  return c;
}

/*Vector addition between matrices a and b, with a scale of k*/
float* VecAdd1(int n, float* a, float* b, float k)
{
  float* c;
  c = malloc(sizeof(float*) * n);
  for (int i=0; i<n; i++)
    c[i] = a[i] + k*b[i];
  return c;
}

#endif
