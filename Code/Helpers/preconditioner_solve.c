#ifndef PRECONDITIONER_SOLVE_C
#define PRECONDITIONER_SOLVE_C

#include <string.h>

#include "poisson_matrix_creaters.c"

struct PCret {
    float* diag;
    float* sol;
};

//conssuder using pointers for the setup and update helper functions
int setupL (int nnz, struct A_csr *L_ptr){
    struct L;
    L_ptr->val = malloc(sizeof(float*) * (nnz));
    L_ptr->col_ind = malloc(sizeof(int*) * (nnz));
    L_ptr->row_ptr = malloc(sizeof(int*) * (nnz + 1));
    L_ptr->val[0] = 1; L_ptr->col_ind[0] = 0; L_ptr->row_ptr[0] = 0;
    return 0;
}

int updateL(){
}

int setupU(int nnz, int n_diag, struct A_csr *U_ptr, struct A_csr *A_ptr){
    U_ptr->val = malloc(sizeof(float*) * (nnz));
    U_ptr->col_ind = malloc(sizeof(int*) * (nnz));
    U_ptr->row_ptr = malloc(sizeof(int*) * (nnz + 1));
    for (int l = 0; l < n_diag; l++){
        U_ptr->val[l] = A_ptr->val[l];
	U_ptr->col_ind[l] = A_ptr->col_ind[l];
	U_ptr->row_ptr[l] = A_ptr->row_ptr[l];
    }
    return 0;
}

int updateU(){
}

// Get ILU in-place factorization
// hardcoded for 1D poisson matrix structure right now
// WE MAY NEED TO USE POINTERS HERE to return the funcs
int getILU(int n, struct A_csr A, struct A_csr arr[]){
    int i, nnz, n_diag;
    struct A_csr L, U;
    n_diag = 2; //hardcoded
    nnz = n*n_diag - 1; //-1 hardcoded

    n_diag_A = 3;
    nnz_A = n*n_diag_A - 2;

    setupL(nnz, &L);
    setupU(nnz, n_diag, &U, &A);

    int cnt_L = 1;
    int cnt_U = n_diag; 
    int cnt_A = n_diag;

    int diag, row, col;
    int index_L_ik, index_U_kj;
    while (cnt_A < nnz_A){
      row = A.row_ptr[cnt_A];
      col = A.col_ind[cnt_A];
      if (row > col){
	//update L
	diag = cnt_A;
	while (A.row_ptr[diag] > A.col_ind[diag]) {diag += 1;}
	L.val[cnt_L] = A.val[cnt_A]/A.val[diag];
	L.col_ind[cnt_L] = col;
	L.row_ptr[cnt_L] = row;
	cnt_L += 1;
      }
      else {
	index_L_ik = cnt_L-1;
	//add 1 to main diag of L
	if (row == col){
	  L.val[cnt_L] = 1;
	  L.col_ind[cnt_L] = col;
	  L.row_ptr[cnt_L] = row;
	  cnt_L += 1;
	}
	index_U_kj = ...;
	//update U
	U.val[cnt_U] = A.val[cnt_A] - L.val[index_L_ik]*U.val[index_U_kj];
	U.row_ptr[cnt_U] = row;
	U.col_ind[cnt_U] = col;
	cnt_U += 1
      }	      
      cnt_A += 1;
    }
    arr[0] = L;
    arr[1] = U;
    return 0;
}

/* Create an ILU-preconditioned A */
//float **create1dILUPoissonMat(int n, float **A){
//    int i, k, j;

    /* need to add some lines that multiply the L and U components */

//    for(i = 1; i < n; i++){
//   		
//      for(k = 0; k < i-1; k++){
//        if (A[i][k] != 0){
//          L[i][k] = A[i][k] / A[k][k];
//            for(j = k + 1; j < n; j++){
//              if (A[i][j] != 0){
//                U[i][j] = A[i][j] - L[i][k] * U[k][j]
//              }
//            } 
//        }
//      }
//    }
//    return A;
//}

// Extract diagonal from CSR matrix
float* getCSRdiagonal(int n, struct A_csr A){
    float* d = create1dZeroVec(n);
    int col=0;
    for (int i = 0; i<n; i++){   
        for (int j=A.row_ptr[i]; j<A.row_ptr[i+1]; j++){
            col = A.col_ind[j];
            if(col==i){
                d[col] = A.val[j];
                col++; 
                if(col>n){break;}   
                }
            }
        }     
    return d; 
}

struct PCret Jacobi_PC_Solve(int n, struct A_csr A, float* r){
    struct PCret pcret;
    float* z = create1dZeroVec(n);
    pcret.diag = getCSRdiagonal(n,A);
    for(int i = 0; i < n; i++){
        z[i] = r[i]/pcret.diag[i];
    }
    pcret.sol = z;
    return pcret;
}

struct PCret ILU_PC_Solve(int n, struct A_csr A, float* r){
    struct PCret pcret;
    struct A_csr  = ...;
    float* z = create1dZeroVec(n);
    pcret.diag = getILU(n,A);
    for(int i = 0; i < n; i++){
        z[i] = r[i]/pcret.diag[i];
    }
    pcret.sol = z;
    return pcret;
}

/*Preconditioner solve*/
struct PCret PC_Solve(int n, float* r, struct A_csr A, char* pctype){

    struct PCret pcret;

    if (strcmp(pctype,"Jacobi")==0){
    	pcret = Jacobi_PC_Solve(n,A,r);
    }

    return pcret;
}


#endif
