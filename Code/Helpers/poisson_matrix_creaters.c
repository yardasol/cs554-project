#ifndef POISSON_MATRIX_CREATORS_C
#define POISSON_MATRIX_CREATORS_C

#include <time.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "spmatrix.h"
#include <string.h>

#include "matrix_arithmetic.c"

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

float **createPoissonILU(int n, float **A){
    int i, k, j;
    for(i = 1; i < n; i++){
        for(k = 0; k < i; k++){
            if (A[i][k] != 0){
                A[i][k] = A[i][k] / A[k][k];
                for(j = k + 1; j < n; j++){
                    if (A[i][j] != 0){
                        A[i][j] = A[i][j] - A[i][k] * A[k][j];
                    }
                }
            }
        }
    }
    return A;
}

float **getLfromPoissonILU(int n, float **ILU){
    int i, j;
    float** L = malloc(sizeof(float*) * n);
    L[0][0] = 1;
    for(i = 0; i < n; i++){
        L[i][i] = 1;
        for(j = 0; j < i; j++){
            L[i][j] = ILU[i][j];
        }
    }
    return L;
}

float **getUfromPoissonILU(int n, float **ILU){
    int i, j;
    float** U = malloc(sizeof(float*) * n);
    for(i = 0; i < n; i++){
        for(j = i; j >= i; j++){
            U[i][j] = ILU[i][j];
        }
    }
    return U;
}



int find_l(int row, int col, int nnz, struct A_csr *A_ptr){
	int l,i;
	int next_row = A_ptr->row_ptr[row+1];
	for (i = (int) A_ptr->row_ptr[row]; i < next_row; i++){
		if ((int) A_ptr->col_ind[i] == col){
			l = i;
			break;
		}
	}
	if (i == next_row){
		l = nnz;
	}
	return l;
}

struct A_csr createIdentityCSR(int n) {
    struct A_csr I;
    const unsigned int nnz = n;
    I.val = malloc(sizeof(float) * nnz);
    for (unsigned int i = 0; i < nnz; i++) {
        I.val[i] = 1.f;
    }

    I.row_ptr = malloc(sizeof(int) * (nnz + 1));
    I.col_ind = malloc(sizeof(int) * n);

    for (unsigned int i = 0; i < nnz; i++) {
        I.row_ptr[i] = i;
        I.col_ind[i] = i;
    }

    I.row_ptr[n] = nnz;

    return I;
}

struct A_csr createPoissonILUCSR(int n, int n_diag, int offset, struct A_csr *A_ptr){
	struct A_csr ILU_csr;
	int nnz = n * n_diag - offset;

	//copy A
	ILU_csr.val = malloc(sizeof(float*) * (nnz));
	for (int a = 0; a < nnz; a++){ILU_csr.val[a] = A_ptr->val[a];}

        ILU_csr.row_ptr = malloc(sizeof(float*) * (n + 1));
        memcpy(ILU_csr.row_ptr, A_ptr->row_ptr, sizeof(float) * (n + 1));
        ILU_csr.col_ind = malloc(sizeof(float*) * nnz);
        memcpy(ILU_csr.col_ind, A_ptr->col_ind, sizeof(float) * nnz);

	int l;
    int i, k, j;
	int ik, kk, ij, kj;
    for(i = 1; i < n; i++){
        for(k = 0; k < i; k++){
			//find index correspoinding to index ik
			ik = find_l(i, k, nnz, A_ptr);

			//find index corresponding to index kk
			kk = find_l(k, k, nnz, A_ptr);
			if ((ik != nnz) && (kk != nnz)){
				ILU_csr.val[ik] = ILU_csr.val[ik] / ILU_csr.val[kk];

				for (j = k + 1; j < n; j++){
					//find index corresponding to index ij
					ij = find_l(i, j, nnz, A_ptr);

					//find index corresponding to index kj
					kj = find_l(k, j, nnz, A_ptr);

					//value at index kj can be zero, resulting in ILU[ij] == A[ij]
					if ((ij != nnz) && (kj == nnz)){
						ILU_csr.val[ij] = ILU_csr.val[ij];
					}
					if ((ij != nnz) && (kj != nnz)){
						ILU_csr.val[ij] = ILU_csr.val[ij] - ILU_csr.val[ik] * ILU_csr.val[kj];
					}
				}
            }
        }
    }
    return ILU_csr;
}

struct A_csr getLfromPoissonILUCSR(int n, int n_diag, int offset, struct A_csr *ILU_ptr){
	struct A_csr L_csr;
	int nnz = n * ((n_diag + 1) / 2) - offset / 2; //adjust for triangular matrix

	L_csr.val = malloc(sizeof(float*) * (nnz));
	L_csr.row_ptr = malloc(sizeof(int*) * (n + 1));
	L_csr.col_ind = malloc(sizeof(int*) * (nnz));

	int row, col;
	int row_ptr, next_row_ptr;
	int l,i;
	l = 0;
	for (row = 0; row < n; row++){
		row_ptr = ILU_ptr->row_ptr[row];
		next_row_ptr = ILU_ptr->row_ptr[row+1];
		for (i = row_ptr; i < next_row_ptr; i++){
			col = ILU_ptr->col_ind[i];
			if (row == col){
				L_csr.val[l] = 1;
				L_csr.col_ind[l] = col;
				L_csr.row_ptr[row+1] = l+1;
				l += 1;
			}
			if (row > col){

				L_csr.val[l] = ILU_ptr->val[i];
				L_csr.col_ind[l] = col;
				l += 1;
			}

		}
	}
	return L_csr;
}

struct A_csr getUfromPoissonILUCSR(int n, int n_diag, int offset, struct A_csr *ILU_ptr){
	struct A_csr U_csr;
	int nnz = n * ((n_diag + 1) / 2) - offset / 2; //adjust for triangular matrix

	U_csr.val = malloc(sizeof(float*) * (nnz));
	U_csr.row_ptr = malloc(sizeof(int*) * (n + 1));
	U_csr.col_ind = malloc(sizeof(int*) * (nnz));

	int row, col;
	int row_ptr, next_row_ptr;
	int l,i;
	l = 0;
	for (row = 0; row < n; row++){
		row_ptr = ILU_ptr->row_ptr[row];
		next_row_ptr = ILU_ptr->row_ptr[row+1];
		for (i = row_ptr; i < next_row_ptr; i++){
			col = ILU_ptr->col_ind[i];
			if (row <= col){
				U_csr.val[l] = ILU_ptr->val[i];
				U_csr.col_ind[l] = col;
				l += 1;
			}

		}

		U_csr.row_ptr[row+1] = l;
	}
	return U_csr;
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

struct A_csr create2dPoissonMatCSR(int n) {
    struct A_csr A_1d = create1dPoissonMatCSR(n);
    struct A_csr I = createIdentityCSR(n);

    struct A_csr* AI = spkron(&A_1d, n, n, &I, n, n, NULL, NULL);
    struct A_csr* IA = spkron(&I, n, n, &A_1d, n, n, NULL, NULL);

    /* Add AI + IA */
    struct dynamic_array coo_entries;
    dynamic_array_initialize(&coo_entries, sizeof(struct A_coo_entry));

    for (unsigned int row = 0; row < n*n; row++) {
        unsigned int ptr_A = AI->row_ptr[row];
        unsigned int ptr_B = IA->row_ptr[row];
        unsigned int end_A = AI->row_ptr[row+1];
        unsigned int end_B = IA->row_ptr[row+1];

        struct A_coo_entry entry;

        while (1) {
            if (ptr_A >= end_A || ptr_B >= end_B) {
                break;
            }

            int col_A = AI->col_ind[ptr_A];
            int col_B = IA->col_ind[ptr_B];

            if (col_A == col_B) {
                entry.val = AI->val[ptr_A] + IA->val[ptr_B];
                entry.row = row;
                entry.col = col_A;

                ptr_A++;
                ptr_B++;
            } else if (col_A < col_B) {
                entry.val = AI->val[ptr_A];
                entry.row = row;
                entry.col = col_A;

                ptr_A++;
            } else {
                entry.val = IA->val[ptr_B];
                entry.row = row;
                entry.col = col_B;

                ptr_B++;
            }

            dynamic_array_push_back(&coo_entries, &entry);
        }

        /* Add leftover entries in row */
        while (ptr_A < end_A) {
            entry.val = AI->val[ptr_A];
            entry.row = row;
            entry.col = AI->col_ind[ptr_A];
            dynamic_array_push_back(&coo_entries, &entry);
            ptr_A++;
        }
        while (ptr_B < end_B) {
            entry.val = IA->val[ptr_B];
            entry.row = row;
            entry.col = IA->col_ind[ptr_B];
            dynamic_array_push_back(&coo_entries, &entry);
            ptr_B++;
        }
    }

    struct A_coo coo_repr;
    coo_repr.data = coo_entries.data;
    coo_repr.nnz = coo_entries.elements;
    coo_repr.rows = n*n;
    coo_repr.cols = n*n;

    struct A_csr csr_repr = *coo_to_csr(&coo_repr);
    dynamic_array_free(&coo_entries);

    return csr_repr;
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
