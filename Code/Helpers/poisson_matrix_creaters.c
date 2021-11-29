#ifndef POISSON_MATRIX_CREATORS_C
#define POISSON_MATRIX_CREATORS_C

#include <time.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

struct A_csr {
    float* val;
    int* col_ind;
    int* row_ptr;
};

struct A_dia {
    float** val;
    float* off;
};

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

struct A_csr createPoissonILUCSR(int n, int n_diag, int offset, struct A_csr *A_ptr){
	struct A_csr ILU_csr;
	int nnz = n * n_diag - offset;

	//copy A
	ILU_csr.val = malloc(sizeof(float*) * (nnz));
	for (int a = 0; a < nnz; a++){ILU_csr.val[a] = A_ptr->val[a];}
	ILU_csr.row_ptr = A_ptr->row_ptr;
	ILU_csr.col_ind = A_ptr->col_ind;

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
				L_csr.row_ptr[row+1] = i+1;
				l += 1;
			}
			if (row > col){
				//printf("row, col: %d, %d\n", row, col);
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
	
		U_csr.row_ptr[row+1] = i;
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
    // Total size of matrix: size^2 * size^2
    float **A;
    int n = size*size;
    int i,j;
    A = malloc(sizeof(float*) * n);

    for (int row = 0; row < n; row++) {
        A[row] = malloc(sizeof(int *) * n);
    }
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            A[i][j] = 0;
            if (i==j) {A[i][j] = 4;}
            else{
                if ((i==j+1)||(i==j-1)){A[i][j] = -1;}
                else{if ((i==j+size)||(i==j-size)) {A[i][j] = -1;}}    
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
