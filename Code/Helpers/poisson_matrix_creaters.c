#ifndef POISSON_MATRIX_CREATORS_C
#define POISSON_MATRIX_CREATORS_C

#include <time.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "spmatrix.h"

/*Create a 1D representation of matrices*/
float **create1dPoissonMat(int n){
    int i, j;
    float **A;
    A = malloc(sizeof(float*) * n);

    for(i = 0; i < n; i++) {
        A[i] = malloc(sizeof(float*) * n);
    }

    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            A[i][j] = 0;
            if (i==j) {A[i][j] = 2;}
            else {if ((i==j+1)||(i==j-1)) {A[i][j] = -1;}}
        }
    }
    return A;
}

struct A_csr create1dPoissonMatCSR(int n){
    int i,  nnz, n_diag;
    struct A_csr A;
    n_diag = 3;
    nnz = n*n_diag-2;
    A.val = malloc(sizeof(float*) * (nnz));
    A.col_ind = malloc(sizeof(int*) * (nnz));
    A.row_ptr = malloc(sizeof(int*) * (n+1));

    A.val[0]=2; A.val[1]=-1;
    A.val[nnz-2]=-1; A.val[nnz-1]=2;
    A.col_ind[0]=0; A.col_ind[1]=1;
    A.col_ind[nnz-2]=n-2; A.col_ind[nnz-1]=n-1;
    int count3 = 2;
    for(i = 2; i < nnz-2; i++){
        A.col_ind[i] = i-count3;
        if(i%3==0){A.val[i]=2;}
        else{A.val[i]=-1;}
        if((i-1)%3==0){count3+=2;}
    }

    A.row_ptr[0]=0; A.row_ptr[1]=A.row_ptr[0]+2;
    for(i = 2; i < n; i++){
        A.row_ptr[i] = A.row_ptr[i-1]+3;
    }
    A.row_ptr[n]=A.row_ptr[n-1]+2;

    return A;
}

struct A_dia create1dPoissonMatDIA(int n){
    int i, j, ndiag;
    struct A_dia A;
    ndiag = 3;

    A.val = malloc(sizeof(float*) * ndiag);

    for(i = 0; i < ndiag; i++) {
        A.val[i] = malloc(sizeof(float*) * n);
    }
    A.off = malloc(sizeof(int*) * (ndiag));

    A.off[0]=-1; A.off[1]=0; A.off[2]=1;

    for(j = 0; j < ndiag; j++){
        for(i = 0; i < n; i++){
            if(j==1){A.val[j][i]=2;}
            else{A.val[j][i]=-1;}
        }
    }

    return A;
}

/*2D representation of matrices*/
float **create2dPoissonMat(int size) {
    // Total size of matrix: size^2 * size^2
    float **A;
    int n = size*size;
    int i,j;
    A = malloc(sizeof(float*) * n);

    for (int row = 0; row < n; row++) {
        A[row] = calloc(n, sizeof(int));
    }

    for (unsigned int x = 0; x < size; x++) {
        for (unsigned int y = 0; y < size; y++) {
            unsigned int idx = y*size + x;
            A[idx][idx] = 4.f;
            if (x > 0) {
                A[idx][idx - 1] = -1.f;
            }
            if (x < size - 1) {
                A[idx][idx + 1] = -1.f;
            }
            if (y > 0) {
                A[idx][idx - size] = -1.f;
            }
            if (y < size - 1) {
                A[idx][idx + size] = -1.f;
            }
        }
    }

    return A;
}

struct A_csr create2dPoissonMatCSR(float ** matrix, int n){
    int nnz, N2D;
    float nnz_per_row;
    N2D = n*n;
    struct A_csr A;
    // To change? Variable on different sizes. We can also just do another
    // pass over matrix to see nnz total.
    nnz_per_row = 4;
    nnz = N2D*nnz_per_row-2;
    A.val = malloc(sizeof(float*) * (nnz));
    A.col_ind = malloc(sizeof(int*) * (nnz));
    A.row_ptr = malloc(sizeof(int*) * (N2D+1));
    A.row_ptr[0] = 0; // row_ptr @[0] is always 0
    // Index of where to insert elements in the .val and .col_ind
    int val_col_idx = 0;
    int nnz_encountered = 0; // Total tally for row_ptr
    int row_ptr_idx = 1; // Index position for row ptr
    // Loop over matrix
    for (int row = 0; row < N2D; row++) {
        for (int col = 0; col < N2D; col++) {
            // If non zero element found
            if (matrix[row][col] != 0) {
                // Add it to the row and col arrays created
                A.val[val_col_idx] = matrix[row][col];
                A.col_ind[val_col_idx] = col;

                nnz_encountered++;
                val_col_idx++;
            }
        }
        // Insert into row_ptr
        A.row_ptr[row_ptr_idx] = nnz_encountered;
        row_ptr_idx++;
    }
    return A;
}

struct A_dia create2dPoissonMatDIA(int n){
    int i, j, ndiag;
    struct A_dia A;
    ndiag = 5;
    int N2D = n*n;
    A.val = malloc(sizeof(float*) * ndiag);

    for(i = 0; i < ndiag; i++) {
        A.val[i] = malloc(sizeof(float*) * N2D);
    }
    A.off = malloc(sizeof(int*) * (ndiag));

    A.off[0] = -n;
    A.off[1] = -1;
    A.off[2] = 0;
    A.off[3] = 1;
    A.off[4] = n;

    // A.off[0]=-1; A.off[1]=0; A.off[2]=1;

    for(j = 0; j < ndiag; j++){
        for(i = 0; i < N2D; i++){
            if(j==2){A.val[j][i]=4;}
            else{A.val[j][i]=-1;}
        }
    }

    return A;
}

/*Zero Vector Creation for 1D and 2D*/
float *create1dZeroVec(int n){
    float *x;
    x = malloc(sizeof(float*) * n);
    for(int i = 0; i < n; i++){
        x[i] = 0;
    }
    return x;
}

float *create2dZeroVec(int size){
    float *x;
    x = malloc(sizeof(int *) * (size*size));

    for (int i = 0; i < size*size; i++) {
        x[i] = 0;
    }
    return x;
}

/*Rand Vector creation*/
float *create1dRandRHS(int n){
    float *x;
    time_t t;
    x = malloc(sizeof(float*) * n);
    /* Intializes random number generator */
    srand((unsigned) time(&t));

    for(int i = 0; i < n; i++){
        x[i] = rand()%10;
    }
    return x;
}

#endif
