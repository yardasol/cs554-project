#ifndef MATRIX_ARITHMETIC_C
#define MATRIX_ARITHMETIC_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "poisson_matrix_creaters.c" // to include A_csr and A_dia structs
// struct A_csr {
//     float* val;
//     float* col_ind;
//     float* row_ptr;
// };

// struct A_dia {
//     float** val;
//     float* off;
// };

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

/**
 * Computes a sparse Kronecker product between two CSR matrices.
 */
struct A_csr* spkron(struct A_csr* A,
                     unsigned int A_num_rows, unsigned int A_num_cols,
                     struct A_csr* B,
                     unsigned int B_num_rows, unsigned int B_num_cols) {
    const unsigned int A_nnz = A->row_ptr[A_num_rows];
    const unsigned int B_nnz = B->row_ptr[B_num_rows];
    const unsigned int C_nnz = A_nnz * B_nnz;
    const unsigned int C_num_rows = A_num_rows * B_num_rows;
    const unsigned int C_num_cols = A_num_cols * B_num_cols;

    float* C_vals = malloc(sizeof(float) * C_nnz);
    float* C_cols = malloc(sizeof(float) * C_nnz);
    float* C_rowptrs = malloc(sizeof(float) * (C_num_rows + 1));

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

    return C;
}

/*Inner product between two matrices*/
float innerProd1(int n, float* a, float* b)
{
  float c;
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
